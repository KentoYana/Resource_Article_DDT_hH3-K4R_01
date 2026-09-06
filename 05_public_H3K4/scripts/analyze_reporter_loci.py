#!/usr/bin/env python3
"""Quantify and plot H3K4 ChIP-seq signal at candidate reporter loci.

The five selected runs remain separate.  CPM values are never pooled across
studies or compared across histone marks as if they were absolute quantities.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import platform
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig


DEFAULT_WORK_ROOT = Path("/Volumes/Garage/Re_analysis/260906_issue69_H3K4")
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR.parent / "output" / "reporter_loci"
WINDOW_FLANK_BP = 2_000
PROMOTER_UPSTREAM_BP = 1_000
PROMOTER_DOWNSTREAM_BP = 200
PROFILE_BIN_BP = 10


@dataclass(frozen=True)
class Target:
    display_name: str
    locus_tag: str
    role: str


@dataclass(frozen=True)
class Gene:
    locus_tag: str
    symbol: str
    chrom: str
    start0: int
    end0: int
    strand: str

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class Run:
    sample_id: str
    study: str
    mark: str


TARGETS = (
    Target("pan-2", "NCU10048", "focal"),
    Target("ad-3A", "NCU03166", "alternative"),
    Target("ad-3B", "NCU03194", "alternative"),
    Target("ad-8", "NCU09789", "alternative"),
    Target("mtr", "NCU06619", "alternative"),
    Target("his-3", "NCU03139", "alternative"),
    Target("csr-1", "NCU00726", "exploratory"),
)

RUNS = (
    Run("Ferraro2021_WT_H3K4me1", "Ferraro et al. 2021", "H3K4me1"),
    Run("Ferraro2021_WT_H3K4me2", "Ferraro et al. 2021", "H3K4me2"),
    Run("Sasaki2014_WT_H3K4me2", "Sasaki et al. 2014", "H3K4me2"),
    Run("Ferraro2021_WT_H3K4me3", "Ferraro et al. 2021", "H3K4me3"),
    Run("Storck2020_WT_H3K4me3", "Storck et al. 2020", "H3K4me3"),
)

VARIANTS = {
    "nonduplicate": "q20.nonduplicate.cpm.bw",
    "all_mapped": "q20.all_mapped.cpm.bw",
}

METRICS = (
    "promoter",
    "gene_body",
    "gene_body_plus_minus_2kb",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-root",
        type=Path,
        default=DEFAULT_WORK_ROOT,
        help=f"Raw reprocessing directory (default: {DEFAULT_WORK_ROOT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for tables and figures (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def load_genes(gff_path: Path) -> dict[str, Gene]:
    genes: dict[str, Gene] = {}
    with gff_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            if attributes.get("gene_biotype") != "protein_coding":
                continue
            locus_tag = attributes.get("locus_tag")
            if not locus_tag:
                continue
            chrom = fields[0]
            if not re.fullmatch(r"CM002(?:23[6-9]|24[0-2])\.1", chrom):
                continue
            gene = Gene(
                locus_tag=locus_tag,
                symbol=attributes.get("gene", attributes.get("Name", locus_tag)),
                chrom=chrom,
                start0=int(fields[3]) - 1,
                end0=int(fields[4]),
                strand=fields[6],
            )
            if locus_tag in genes:
                raise ValueError(f"Duplicate gene feature for {locus_tag}")
            genes[locus_tag] = gene
    return genes


def interval_for_metric(
    gene: Gene, metric: str, chrom_length: int
) -> tuple[int, int]:
    if metric == "gene_body":
        start, end = gene.start0, gene.end0
    elif metric == "gene_body_plus_minus_2kb":
        start = gene.start0 - WINDOW_FLANK_BP
        end = gene.end0 + WINDOW_FLANK_BP
    elif metric == "promoter":
        if gene.strand == "+":
            start = gene.start0 - PROMOTER_UPSTREAM_BP
            end = gene.start0 + PROMOTER_DOWNSTREAM_BP
        elif gene.strand == "-":
            start = gene.end0 - PROMOTER_DOWNSTREAM_BP
            end = gene.end0 + PROMOTER_UPSTREAM_BP
        else:
            raise ValueError(f"Unsupported strand for {gene.locus_tag}: {gene.strand}")
    else:
        raise ValueError(f"Unknown metric: {metric}")
    return max(0, start), min(chrom_length, end)


def exact_mean(
    bigwig: pyBigWig.pyBigWig, gene: Gene, metric: str, chrom_length: int
) -> float:
    start, end = interval_for_metric(gene, metric, chrom_length)
    result = bigwig.stats(gene.chrom, start, end, type="mean", exact=True)[0]
    return math.nan if result is None else float(result)


def empirical_midrank(value: float, population: list[float]) -> tuple[float, float, int]:
    finite = sorted(item for item in population if math.isfinite(item))
    if not math.isfinite(value) or not finite:
        return math.nan, math.nan, len(finite)
    left = bisect.bisect_left(finite, value)
    right = bisect.bisect_right(finite, value)
    ties = right - left
    percentile = 100.0 * (left + 0.5 * ties) / len(finite)
    rank_desc = (len(finite) - right) + (ties + 1) / 2.0
    return percentile, rank_desc, len(finite)


def region_details(
    bigwig: pyBigWig.pyBigWig, gene: Gene, metric: str, chrom_length: int
) -> dict[str, float | int]:
    start, end = interval_for_metric(gene, metric, chrom_length)
    values = np.asarray(bigwig.values(gene.chrom, start, end, numpy=True), dtype=float)
    valid = np.isfinite(values)
    covered = int(valid.sum())
    region_length = int(values.size)
    if covered:
        valid_values = values[valid]
        mean_cpm = float(valid_values.mean())
        max_cpm = float(valid_values.max())
        nonzero_fraction = float(np.count_nonzero(valid_values > 0) / covered)
        exact_zero_bases = int(np.count_nonzero(valid_values == 0))
    else:
        mean_cpm = math.nan
        max_cpm = math.nan
        nonzero_fraction = math.nan
        exact_zero_bases = 0
    return {
        "region_start_1based": start + 1,
        "region_end_1based": end,
        "region_length_bp": region_length,
        "covered_bases": covered,
        "missing_bases": region_length - covered,
        "exact_zero_bases": exact_zero_bases,
        "mean_cpm": mean_cpm,
        "max_cpm": max_cpm,
        "nonzero_fraction": nonzero_fraction,
    }


def format_number(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.6f}"
    return str(value)


def write_locus_table(path: Path, genes: dict[str, Gene]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "display_name",
                "locus_tag",
                "role",
                "chrom",
                "start_1based",
                "end_1based",
                "strand",
                "length_bp",
                "gff_gene_name",
            ]
        )
        for target in TARGETS:
            gene = genes[target.locus_tag]
            writer.writerow(
                [
                    target.display_name,
                    target.locus_tag,
                    target.role,
                    gene.chrom,
                    gene.start0 + 1,
                    gene.end0,
                    gene.strand,
                    gene.length,
                    gene.symbol,
                ]
            )


def collect_summary(
    genes: dict[str, Gene], bigwig_root: Path
) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    genome_genes = list(genes.values())
    target_lookup = {target.locus_tag: target for target in TARGETS}

    for variant, suffix in VARIANTS.items():
        for run in RUNS:
            bigwig_path = bigwig_root / f"{run.sample_id}.{suffix}"
            if not bigwig_path.is_file():
                raise FileNotFoundError(f"Missing bigWig: {bigwig_path}")
            with pyBigWig.open(str(bigwig_path)) as bigwig:
                chrom_lengths = bigwig.chroms()
                for gene in genome_genes:
                    if gene.chrom not in chrom_lengths:
                        raise ValueError(f"{gene.chrom} is absent from {bigwig_path.name}")

                genome_values: dict[str, list[float]] = {metric: [] for metric in METRICS}
                gene_values: dict[tuple[str, str], float] = {}
                for gene in genome_genes:
                    for metric in METRICS:
                        value = exact_mean(
                            bigwig, gene, metric, chrom_lengths[gene.chrom]
                        )
                        genome_values[metric].append(value)
                        gene_values[(gene.locus_tag, metric)] = value

                for locus_tag, target in target_lookup.items():
                    gene = genes[locus_tag]
                    for metric in METRICS:
                        details = region_details(
                            bigwig, gene, metric, chrom_lengths[gene.chrom]
                        )
                        value = gene_values[(locus_tag, metric)]
                        percentile, rank_desc, population_n = empirical_midrank(
                            value, genome_values[metric]
                        )
                        reporter_values = [
                            gene_values[(item.locus_tag, metric)] for item in TARGETS
                        ]
                        _, reporter_rank, reporter_n = empirical_midrank(
                            value, reporter_values
                        )
                        summary.append(
                            {
                                "sample_id": run.sample_id,
                                "study": run.study,
                                "mark": run.mark,
                                "track_variant": variant,
                                "bigwig_file": bigwig_path.name,
                                "display_name": target.display_name,
                                "locus_tag": target.locus_tag,
                                "role": target.role,
                                "chrom": gene.chrom,
                                "gene_start_1based": gene.start0 + 1,
                                "gene_end_1based": gene.end0,
                                "strand": gene.strand,
                                "metric": metric,
                                **details,
                                "genome_percentile_midrank": percentile,
                                "genome_rank_desc": rank_desc,
                                "genome_population_n": population_n,
                                "reporter_set_rank_desc": reporter_rank,
                                "reporter_set_n": reporter_n,
                            }
                        )
    return summary


def write_summary_table(path: Path, summary: list[dict[str, object]]) -> None:
    fieldnames = list(summary[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in summary:
            writer.writerow({key: format_number(value) for key, value in row.items()})


def summary_index(
    summary: list[dict[str, object]],
) -> dict[tuple[str, str, str, str], dict[str, object]]:
    return {
        (
            str(row["sample_id"]),
            str(row["track_variant"]),
            str(row["locus_tag"]),
            str(row["metric"]),
        ): row
        for row in summary
    }


def profile_for_gene(
    bigwig: pyBigWig.pyBigWig, gene: Gene, chrom_length: int
) -> tuple[np.ndarray, np.ndarray]:
    start, end = interval_for_metric(
        gene, "gene_body_plus_minus_2kb", chrom_length
    )
    n_bins = max(1, math.ceil((end - start) / PROFILE_BIN_BP))
    raw = bigwig.stats(gene.chrom, start, end, nBins=n_bins, type="mean")
    values = np.asarray(
        [math.nan if value is None else float(value) for value in raw], dtype=float
    )
    if gene.strand == "-":
        values = values[::-1]
    x = np.linspace(
        -WINDOW_FLANK_BP + (end - start) / n_bins / 2,
        gene.length + WINDOW_FLANK_BP - (end - start) / n_bins / 2,
        n_bins,
    )
    return x, values


def save_figure(fig: plt.Figure, base_path: Path) -> None:
    pdf_path = base_path.parent / f"{base_path.name}.pdf"
    png_path = base_path.parent / f"{base_path.name}.png"
    fig.savefig(
        pdf_path,
        metadata={"Creator": "analyze_reporter_loci.py", "CreationDate": None},
    )
    fig.savefig(png_path, dpi=220)
    plt.close(fig)


def plot_profiles(
    variant: str,
    genes: dict[str, Gene],
    bigwig_root: Path,
    output_dir: Path,
    lookup: dict[tuple[str, str, str, str], dict[str, object]],
) -> None:
    profiles: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
    row_maxima: dict[str, float] = {}
    suffix = VARIANTS[variant]
    for run in RUNS:
        bigwig_path = bigwig_root / f"{run.sample_id}.{suffix}"
        row_max = 0.0
        with pyBigWig.open(str(bigwig_path)) as bigwig:
            chrom_lengths = bigwig.chroms()
            for target in TARGETS:
                gene = genes[target.locus_tag]
                x, values = profile_for_gene(
                    bigwig, gene, chrom_lengths[gene.chrom]
                )
                profiles[(run.sample_id, target.locus_tag)] = (x, values)
                finite = values[np.isfinite(values)]
                if finite.size:
                    row_max = max(row_max, float(finite.max()))
        row_maxima[run.sample_id] = max(1.0, row_max * 1.03)

    fig, axes = plt.subplots(
        len(RUNS), len(TARGETS), figsize=(22, 12.5), squeeze=False
    )
    for row_index, run in enumerate(RUNS):
        y_max = row_maxima[run.sample_id]
        for column_index, target in enumerate(TARGETS):
            gene = genes[target.locus_tag]
            ax = axes[row_index][column_index]
            x, values = profiles[(run.sample_id, target.locus_tag)]
            ax.axvspan(
                -PROMOTER_UPSTREAM_BP,
                PROMOTER_DOWNSTREAM_BP,
                color="#f2c14e",
                alpha=0.18,
                linewidth=0,
            )
            ax.axvspan(0, gene.length, color="#7a7a7a", alpha=0.12, linewidth=0)
            ax.axvline(0, color="#555555", linestyle="--", linewidth=0.6)
            ax.plot(x, values, color="#1546a0", linewidth=0.8)
            ax.fill_between(x, 0, values, color="#2f66c5", alpha=0.55)
            ax.set_xlim(-WINDOW_FLANK_BP, gene.length + WINDOW_FLANK_BP)
            ax.set_ylim(0, y_max)
            ax.set_yticks([0, y_max])
            ax.tick_params(axis="both", labelsize=6, length=2)
            ax.set_xticks(
                [-WINDOW_FLANK_BP, 0, gene.length, gene.length + WINDOW_FLANK_BP],
                ["-2 kb", "TSS", "TES", "+2 kb"],
            )
            body = lookup[
                (run.sample_id, variant, target.locus_tag, "gene_body")
            ]
            promoter = lookup[
                (run.sample_id, variant, target.locus_tag, "promoter")
            ]
            ax.text(
                0.98,
                0.94,
                "body pct {:.0f} · promoter pct {:.0f}".format(
                    float(body["genome_percentile_midrank"]),
                    float(promoter["genome_percentile_midrank"]),
                ),
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=5.4,
                color="#333333",
            )
            if row_index == 0:
                ax.set_title(
                    f"{target.display_name}\n{target.locus_tag} ({gene.strand})",
                    fontsize=9,
                )
            if column_index == 0:
                ax.set_ylabel(
                    f"{run.study}\n{run.mark}\nCPM",
                    fontsize=7.5,
                )
            else:
                ax.set_yticklabels([])
    fig.suptitle(
        f"Candidate reporter loci: gene-oriented H3K4me profiles ({variant})",
        fontsize=14,
    )
    fig.text(
        0.5,
        0.012,
        "Each row has one shared y-axis across all seven loci; y-axes are not shared between runs. "
        "Gold: promoter (-1 kb to +200 bp); gray: annotated gene body. Annotations are within-run genome-wide midrank percentiles.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0.02, 0.035, 1, 0.965), h_pad=1.0, w_pad=0.5)
    save_figure(fig, output_dir / f"reporter_locus_profiles.{variant}")


def plot_percentile_heatmaps(
    variant: str,
    output_dir: Path,
    lookup: dict[tuple[str, str, str, str], dict[str, object]],
) -> None:
    fig = plt.figure(figsize=(17, 6.2))
    grid = fig.add_gridspec(
        1, len(METRICS) + 1, width_ratios=[1, 1, 1, 0.045], wspace=0.42
    )
    axes = [fig.add_subplot(grid[0, index]) for index in range(len(METRICS))]
    metric_titles = {
        "promoter": "Promoter (-1 kb/+200 bp)",
        "gene_body": "Gene body",
        "gene_body_plus_minus_2kb": "Gene body ±2 kb",
    }
    image = None
    for metric_index, metric in enumerate(METRICS):
        ax = axes[metric_index]
        matrix = np.asarray(
            [
                [
                    float(
                        lookup[
                            (run.sample_id, variant, target.locus_tag, metric)
                        ]["genome_percentile_midrank"]
                    )
                    for target in TARGETS
                ]
                for run in RUNS
            ]
        )
        image = ax.imshow(matrix, vmin=0, vmax=100, cmap="viridis", aspect="auto")
        ax.set_title(metric_titles[metric], fontsize=10)
        ax.set_xticks(
            range(len(TARGETS)), [target.display_name for target in TARGETS], rotation=45, ha="right"
        )
        ax.set_yticks(
            range(len(RUNS)),
            [f"{run.study}\n{run.mark}" for run in RUNS],
        )
        ax.tick_params(labelsize=7)
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                value = matrix[row_index, column_index]
                color = "white" if value < 45 or value > 80 else "black"
                ax.text(
                    column_index,
                    row_index,
                    f"{value:.0f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=color,
                )
    assert image is not None
    colorbar_axis = fig.add_subplot(grid[0, -1])
    colorbar = fig.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("Within-run genome-wide percentile (midrank)", fontsize=8)
    fig.suptitle(
        f"Candidate reporter H3K4me signal percentiles ({variant})", fontsize=13
    )
    fig.subplots_adjust(left=0.12, right=0.95, top=0.84, bottom=0.21)
    save_figure(fig, output_dir / f"reporter_locus_percentiles.{variant}")


def write_analysis_summary(
    path: Path,
    lookup: dict[tuple[str, str, str, str], dict[str, object]],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Reporter-locus H3K4me descriptive summary\n\n")
        handle.write(
            "This table reports the nonduplicate tracks for `pan-2`. Percentiles are "
            "calculated separately for every run and region definition across all annotated "
            "protein-coding genes on the seven nuclear chromosomes.\n\n"
        )
        handle.write(
            "| Study | Mark | Region | Mean CPM | Nonzero fraction | Genome percentile | Reporter-set rank |\n"
        )
        handle.write("|---|---|---|---:|---:|---:|---:|\n")
        for run in RUNS:
            for metric in METRICS:
                row = lookup[(run.sample_id, "nonduplicate", "NCU10048", metric)]
                handle.write(
                    "| {} | {} | {} | {:.4f} | {:.3f} | {:.1f} | {:.1f}/7 |\n".format(
                        run.study,
                        run.mark,
                        metric,
                        float(row["mean_cpm"]),
                        float(row["nonzero_fraction"]),
                        float(row["genome_percentile_midrank"]),
                        float(row["reporter_set_rank_desc"]),
                    )
                )
        handle.write("\n## Interpretation boundaries\n\n")
        handle.write(
            "- These are descriptive CPM summaries, not statistical tests; no selected mark has biological replication within every study.\n"
            "- Nonzero coverage is not equivalent to enrichment because no matched input is used.\n"
            "- Studies and marks remain separate; absolute CPM values are not pooled or treated as directly exchangeable.\n"
            "- Exact zero bases are counted separately from missing bases in the long-form TSV.\n"
            "- `csr-1` is exploratory and is not treated as an established forward-mutation reporter.\n"
        )


def write_parameters(path: Path, gff_path: Path) -> None:
    payload = {
        "analysis_script": Path(__file__).name,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "pyBigWig": pyBigWig.__version__,
        "reference_gff": gff_path.name,
        "reference_accession": "GCA_000182925.2",
        "input_bigwigs": [
            f"{run.sample_id}.{suffix}"
            for variant, suffix in VARIANTS.items()
            for run in RUNS
        ],
        "runs": [asdict(run) for run in RUNS],
        "targets": [asdict(target) for target in TARGETS],
        "normalization": "CPM; inherited from the uniformly generated bigWig inputs",
        "minimum_mapping_quality": 20,
        "profile_bin_bp": PROFILE_BIN_BP,
        "window_flank_bp": WINDOW_FLANK_BP,
        "promoter_definition": "strand-aware TSS -1000 bp through +200 bp",
        "gene_body_definition": "full NCBI GFF gene feature",
        "ranking_population": "protein-coding genes on the seven NC12 nuclear chromosomes",
        "percentile_method": "empirical midrank, calculated independently for each run, track variant, and region metric",
        "exact_zero_treatment": "included as zero; counted separately in reporter_locus_signal_summary.tsv",
        "missing_treatment": "excluded from regional means and counted separately; NA retained when a region has no covered bases",
        "random_seed": None,
        "pooling": "none",
    }
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_checksums(output_dir: Path) -> None:
    checksum_path = output_dir / "output_sha256.txt"
    paths = sorted(
        path
        for path in output_dir.iterdir()
        if path.is_file() and path.name != checksum_path.name
    )
    with checksum_path.open("w", encoding="utf-8") as handle:
        for path in paths:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{digest}  {path.name}\n")


def main() -> int:
    args = parse_args()
    work_root = args.work_root.resolve()
    output_dir = args.output_dir.resolve()
    gff_path = work_root / "reference" / "GCA_000182925.2_NC12_genomic.gff"
    bigwig_root = work_root / "bigwig"
    if not gff_path.is_file():
        raise FileNotFoundError(f"Missing reference GFF: {gff_path}")
    if not bigwig_root.is_dir():
        raise FileNotFoundError(f"Missing bigWig directory: {bigwig_root}")
    output_dir.mkdir(parents=True, exist_ok=True)

    genes = load_genes(gff_path)
    missing_targets = [target.locus_tag for target in TARGETS if target.locus_tag not in genes]
    if missing_targets:
        raise ValueError(f"Target loci absent from GFF: {', '.join(missing_targets)}")

    write_locus_table(output_dir / "reporter_loci.tsv", genes)
    summary = collect_summary(genes, bigwig_root)
    write_summary_table(output_dir / "reporter_locus_signal_summary.tsv", summary)
    lookup = summary_index(summary)
    for variant in VARIANTS:
        plot_profiles(variant, genes, bigwig_root, output_dir, lookup)
        plot_percentile_heatmaps(variant, output_dir, lookup)
    write_analysis_summary(output_dir / "analysis_summary.md", lookup)
    write_parameters(output_dir / "analysis_parameters.json", gff_path)
    write_checksums(output_dir)
    print(f"Analyzed {len(TARGETS)} reporter loci across {len(RUNS)} runs and {len(VARIANTS)} track variants")
    print(f"Genome-wide ranking population: {len(genes)} protein-coding genes")
    print(f"Outputs: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
