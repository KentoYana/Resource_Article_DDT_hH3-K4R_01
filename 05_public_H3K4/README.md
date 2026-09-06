# Public H3K4 ChIP-seq dataset

This directory records the public wild-type *Neurospora crassa* H3K4
ChIP-seq datasets selected for locus-level inspection. It includes a raw-read
reprocessing workflow and the resulting reporter-locus summaries.

## Dataset identification

The selected series is
[GSE154497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154497),
"Chromatin accessibility profiling in *Neurospora crassa* reveals features
associated with accessible and inaccessible chromatin." The series contains
wild-type mycelial ChIP-seq samples for H3K4me1, H3K4me2, and H3K4me3.

The H3K4me2 sample relevant to inspection of the `pan-2` locus is:

- GEO sample: [GSM4672246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672246)
- SRA experiment: [SRX8737830](https://www.ncbi.nlm.nih.gov/sra/SRX8737830)
- SRA run: [SRR12229307](https://www.ncbi.nlm.nih.gov/sra/SRR12229307)

The series includes a public CPM-normalized bigWig for each of the three H3K4
marks. These deposited tracks are the preferred starting point for reproducing
a browser view because they retain the processing selected by the original
authors.

The Ferraro et al. accession and sample crosswalk is recorded in
[`dataset/accessions.tsv`](dataset/accessions.tsv).

## Raw-read reprocessing workflow

[`scripts/reprocess_selected_h3k4_chipseq.sh`](scripts/reprocess_selected_h3k4_chipseq.sh)
downloads and uniformly processes five selected wild-type H3K4 ChIP-seq runs:

- Ferraro et al. (2021): SRR12229306, SRR12229307, and SRR12229308;
- Sasaki et al. (2014): SRR1295547; and
- Storck et al. (2020): SRR12202381.

By default, all downloaded data and generated outputs are kept outside this
repository in `/Volumes/Garage/Re_analysis/260906_issue69_H3K4`. The workflow
can be inspected without downloading data by running:

```sh
./05_public_H3K4/scripts/reprocess_selected_h3k4_chipseq.sh check
```

Run `--help` to list the resumable stages. The workflow obtains raw SRA data
and the NC12 reference anew, performs raw and trimmed FastQC, maps all runs to
the same reference, marks duplicate reads without destroying them, and creates
both duplicate-retaining and nonduplicate CPM bigWig tracks. It intentionally
does not pool studies or call peaks without a matched input.

After bigWig generation, the `reporters` stage runs
[`scripts/analyze_reporter_loci.R`](scripts/analyze_reporter_loci.R). It
compares `pan-2` with `ad-3A`, `ad-3B`, `ad-8`, `mtr`, `his-3`, and the
exploratory `csr-1` locus. Coordinates and strand are read from the downloaded
NC12 GFF. The analysis reports strand-aware promoter (-1 kb to +200 bp), gene
body, and gene-body-plus-or-minus-2-kb CPM values, while retaining exact zero
and missing bases as distinct states. Percentiles are calculated independently
for each run and track variant across all protein-coding genes on the seven
nuclear chromosomes.

The reporter analysis requires R with `tidyverse`, `tikzDevice`, `patchwork`,
`here`, `IRanges`, `rtracklayer`, `digest`, and `jsonlite`. The workflow's
`check` stage verifies these packages in the active R library before starting
an analysis.

As with the other Resource analyses, the script can also be opened in RStudio
and run with **Source**. It uses `here::i_am()` to locate the repository root,
so open this repository as the RStudio working project before sourcing
`scripts/analyze_reporter_loci.R`. The figures are first drawn on the standard
R graphics device for inspection in RStudio's Plots pane, then drawn again with
`tikzDevice` for final output. The default bigWig and reference location is the
external SSD path
`/Volumes/Garage/Re_analysis/260906_issue69_H3K4`.
Set `H3K4_WORK_ROOT` before sourcing to use a different analysis directory:

```r
Sys.setenv(H3K4_WORK_ROOT = "/path/to/issue69_H3K4")
source("05_public_H3K4/scripts/analyze_reporter_loci.R")
```

`H3K4_OUTPUT_DIR` may similarly override the output directory. With neither
variable set, output is written to `05_public_H3K4/output/reporter_loci`.

```sh
./05_public_H3K4/scripts/reprocess_selected_h3k4_chipseq.sh reporters
```

Tables and fixed-scale TikZ figures are written to
[`output/reporter_loci`](output/reporter_loci). Each profile-plot row shares
one y-axis across all seven loci, but y-axes are not shared across studies or
histone marks. Profile figures are separated into H3K4me1, H3K4me2, and
H3K4me3 outputs so that labels and annotations remain legible at the target
width. Both duplicate-retaining and nonduplicate results are kept; the studies
are never pooled. The figure files are LaTeX fragments containing
`tikzpicture` environments generated from `ggplot2`/`patchwork` objects with
`tikzDevice`; they require TikZ when included in a document. TikZ figures use
the same 7.5-inch width and `lwdUnit = 72.27 / 96` as the reversion-assay
figure. The heatmap uses a continuous `#0068b7`-to-`#f39800` scale with its
color bar below the panels.

## Experimental metadata

The GEO records describe the samples as wild-type mycelial cultures grown in
Vogel's medium with 1.5% sucrose. The associated paper specifies shaking at
32 degrees C for 18 hours. Libraries were sequenced on an Illumina NextSeq 500
as single-end reads with a nominal length of 75 bp.

The deposited processing description reports:

- Trim Galore with `--length 20 --fastqc`;
- BWA-MEM 0.7.15 with `-M` or Bowtie2 2.4.1 with `--very-sensitive`;
- SAMtools 1.3.1;
- deepTools 3.3.1 with `--normalizeUsing CPM`; and
- genome build `GCA_000182925.2 (NC12)`.

GEO records the Trim Galore version as "version 4.0". This value is reproduced
as deposited and has not been silently corrected. GEO also does not identify
which of the two listed aligners was used for each H3K4 sample.

## Publication

Ferraro AR, Ameri AJ, Lu Z, Kamei M, Schmitz RJ, and Lewis ZA (2021).
"Chromatin accessibility profiling in *Neurospora crassa* reveals molecular
features associated with accessible and inaccessible chromatin."
*BMC Genomics* 22:459.

- DOI: [10.1186/s12864-021-07774-0](https://doi.org/10.1186/s12864-021-07774-0)
- PubMed: [34147068](https://pubmed.ncbi.nlm.nih.gov/34147068/)
- BioProject: [PRJNA646493](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA646493)
- SRA study: [SRP272099](https://www.ncbi.nlm.nih.gov/sra/?term=SRP272099)

## Interpretation limits

- Each H3K4 mark is represented by one GEO sample and one SRA run; biological
  replication is not available for these marks.
- The H3K4me2 and H3K4me3 runs are shallow (649,309 and 170,871 spots,
  respectively). A weak or absent browser signal is not a general proof that a
  locus can never be methylated.
- The series contains an input sample, GSM4672252/SRR12229313, but its genotype
  is `hH3-3xFLAG`, whereas the three H3K4-mark samples are recorded as wild
  type. It is therefore not treated here as a matched input for the H3K4-mark
  ChIP-seq samples.
- CPM values from different marks should not be compared as absolute
  methylation amounts because the samples use different antibodies and have
  different sequencing depths.
- The samples are basal, untreated mycelial cultures. They cannot determine
  whether H3K4 methylation changes after UV irradiation or replication stress.
- A locus-level pattern can support a condition-specific descriptive statement,
  but it cannot by itself establish a direct or indirect causal mechanism.

The RNA-seq runs SRR5177529 and SRR5177530 are not part of GSE154497. They come
from a separate study and are not matched to these ChIP-seq cultures; they are
therefore not included in the ChIP-seq accession table.

## Selection boundary

GSE154497 is selected as the best-supported public source for the wild-type
H3K4 browser tracks. A browser screenshot without track metadata is not, by
itself, sufficient to prove that the displayed track is byte-for-byte identical
to the deposited bigWig. This identification therefore records the public
dataset and the evidence-supported correspondence without treating an
unlabelled screenshot as the primary data source.

## Primary metadata sources

- [GSE154497 series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154497)
- [GSM4672245, H3K4me1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672245)
- [GSM4672246, H3K4me2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672246)
- [GSM4672247, H3K4me3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672247)
- [SRA Run Selector for SRP272099](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP272099)
- [Ferraro et al. (2021)](https://doi.org/10.1186/s12864-021-07774-0)
