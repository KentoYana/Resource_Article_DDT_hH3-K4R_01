# Public H3K4 ChIP-seq dataset

This directory contains the public wild-type *Neurospora crassa* H3K4
ChIP-seq datasets selected for locus-level inspection, a raw-read reprocessing
workflow, and the resulting reporter-locus and repair/DDT-gene summaries.

## Dataset identification

The selected series,
[GSE154497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154497), is
"Chromatin accessibility profiling in *Neurospora crassa* reveals features
associated with accessible and inaccessible chromatin." It contains wild-type
mycelial ChIP-seq samples for H3K4me1, H3K4me2, and H3K4me3.

The H3K4me2 sample relevant to inspection of the `pan-2` locus is:

- GEO sample: [GSM4672246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672246)
- SRA experiment: [SRX8737830](https://www.ncbi.nlm.nih.gov/sra/SRX8737830)
- SRA run: [SRR12229307](https://www.ncbi.nlm.nih.gov/sra/SRR12229307)

The series includes a public CPM-normalized bigWig for each H3K4 mark. These
deposited tracks retain the original authors' processing and are the preferred
starting point for reproducing a browser view.

The Ferraro et al. accession and sample crosswalk is recorded in
[`dataset/accessions.tsv`](dataset/accessions.tsv).

## Raw-read reprocessing workflow

[`scripts/reprocess_selected_h3k4_chipseq.sh`](scripts/reprocess_selected_h3k4_chipseq.sh)
downloads and uniformly processes five selected wild-type H3K4 ChIP-seq runs:

- Ferraro et al. (2021): SRR12229306, SRR12229307, and SRR12229308;
- Sasaki et al. (2014): SRR1295547; and
- Storck et al. (2020): SRR12202381.

By default, downloaded data and generated outputs remain outside this
repository in `/Volumes/Garage/Re_analysis/260906_issue69_H3K4`. To inspect the
workflow without downloading data, run:

```sh
./05_public_H3K4/scripts/reprocess_selected_h3k4_chipseq.sh check
```

Run `--help` to list the resumable stages. The workflow obtains raw SRA data
and the NC12 reference, performs raw and trimmed FastQC, and maps all runs to
the same reference. It marks duplicate reads without discarding them and
creates duplicate-retaining and nonduplicate CPM bigWig tracks. Studies are
not pooled, and peaks are not called without a matched input.

After bigWig generation, the `reporters` stage runs
[`scripts/analyze_reporter_loci.R`](scripts/analyze_reporter_loci.R). It
compares `pan-2` with `ad-3A`, `ad-3B`, `ad-8`, `his-3`, `mtr`, and the
exploratory `csr-1` locus. Coordinates and strand come from the downloaded NC12
GFF. The analysis reports strand-aware promoter (-1 kb to +200 bp), gene-body,
and gene-body-plus-or-minus-2-kb CPM values. Exact zero and missing bases remain
distinct. Percentiles are calculated separately for each run and track variant
across all protein-coding genes on the seven nuclear chromosomes.

The reporter analysis requires R with `tidyverse`, `tikzDevice`, `patchwork`,
`here`, `IRanges`, `rtracklayer`, `digest`, and `jsonlite`. The workflow's
`check` stage verifies these packages in the active R library before starting
an analysis.

The script can also be sourced in RStudio. It uses `here::i_am()` to locate the
repository root, so open this repository as the RStudio project before sourcing
`scripts/analyze_reporter_loci.R`. The script first draws figures on the
standard R graphics device for inspection in the Plots pane. It then redraws
them with `tikzDevice` for final output. By default, bigWig and reference files
are read from the external SSD path
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
[`output/reporter_loci`](output/reporter_loci). Each profile-plot row shares a
y-axis across all seven loci; y-axes are not shared across studies or histone
marks. Separate H3K4me1, H3K4me2, and H3K4me3 files keep labels and annotations
legible at the target width. Duplicate-retaining and nonduplicate results are
both retained, and studies are not pooled.

The figure files are LaTeX fragments containing `tikzpicture` environments
generated from `ggplot2`/`patchwork` objects with `tikzDevice`. They require
TikZ when included in a document. The files use the same 7.5-inch width and
`lwdUnit = 72.27 / 96` as the reversion-assay figure. Following the color-scale
example at
<https://okumuralab.org/~okumura/stat/colors.html>, the heatmap maps 0 to
`#0068b7`, 50 to white, and 100 to `#f39800`, with its color bar below the
panels.

### Ferraro-only manuscript display

The revision manuscript uses the Ferraro et al. (2021) series as its primary
display because all three H3K4 marks were generated within that study. Generate
the dedicated reporter-locus output without overwriting the cross-study
sensitivity analysis:

```sh
H3K4_STUDY_SELECTION=Ferraro2021 \
  Rscript 05_public_H3K4/scripts/analyze_reporter_loci.R
```

This writes to `output/reporter_loci_ferraro`. For the repair/DDT target set,
run:

```sh
H3K4_TARGET_SET=repair H3K4_STUDY_SELECTION=Ferraro2021 \
  Rscript 05_public_H3K4/scripts/analyze_reporter_loci.R
```

This writes to `output/repair_genes_ferraro`. The existing `reporter_loci` and
`repair_genes` directories retain the separately displayed Sasaki and Storck
runs as sensitivity analyses. Reads and CPM values are not pooled across
studies or histone marks.
Percentile-heatmap LaTeX filenames record both the track variant and selected
dataset, for example
`reporter_locus_percentiles_nonduplicate_Ferraro2021.tex`.

## Repair/DDT-related gene inspection

For Reviewer 1's question about a possible indirect effect through repair-gene
regulation, the same uniformly processed H3K4me tracks can be inspected at the
eight genes tested directly in the genetic-interaction experiments: `mus-9`,
`uvs-2`, `mus-26`, `polh`, `qde-3`, `recQ2`, `mei-3`, and `mus-11`. No new
download, alignment, or normalization is performed for this target set.

Run only this analysis stage from the shell with:

```sh
./05_public_H3K4/scripts/reprocess_selected_h3k4_chipseq.sh repairs
```

To run the same analysis in RStudio, set the target collection before sourcing
the script:

```r
Sys.setenv(H3K4_TARGET_SET = "repair")
source("05_public_H3K4/scripts/analyze_reporter_loci.R")
```

Results are written to [`output/repair_genes`](output/repair_genes). The output
uses the same promoter, gene-body, and plus-or-minus-2-kb definitions,
within-run genome-wide percentiles, duplicate-retaining/nonduplicate variants,
7.5-inch TikZ width, and blue-white-orange scale as the reporter-locus output.
The studies are displayed separately and are not pooled. These basal wild-type
profiles cannot determine transcriptional effects of `hH3-K4R` or changes after
UV irradiation or replication stress.

The nonduplicate promoter profiles are heterogeneous rather than uniformly
H3K4me-rich. `uvs-2`, `mus-26`, and `mus-11` have high H3K4me1 and H3K4me2
promoter percentiles; `mus-26` and `mus-11` also have high H3K4me3 promoter
percentiles. By contrast, `recQ2` has low H3K4me2 and H3K4me3 promoter
percentiles, while `qde-3` and `mei-3` are not consistently high for H3K4me3.
These descriptive patterns make a gene-specific indirect effect plausible but
do not establish altered expression or causality.

## Nucleosome occupancy and accessibility workflow

Low H3K4me ChIP-seq signal can reflect a low fraction of methylated H3 or low
local nucleosome occupancy. To examine these alternatives with public data,
[`scripts/reprocess_selected_nucleosome_occupancy.sh`](scripts/reprocess_selected_nucleosome_occupancy.sh)
uniformly reprocesses the paired-end raw reads listed in
[`dataset/nucleosome_accessions.tsv`](dataset/nucleosome_accessions.tsv):

- three wild-type MNase-seq biological replicates and two wild-type total-H3
  ChIP-seq biological replicates from Kamei et al. (2021), GSE150758; and
- two wild-type ATAC-seq biological replicates from Ferraro et al. (2021),
  GSE154497.

By default, the workflow stores raw and intermediate data on the external SSD
at `/Volumes/Garage/Re_analysis/260907_issue69_nucleosome`. It downloads NC12
and all seven SRA runs, then performs paired-end trimming and alignment. The
workflow retains proper primary pairs with MAPQ of at least 20 on the seven
nuclear chromosomes and marks duplicates without discarding them. It generates
both duplicate-retaining and nonduplicate tracks.

MNase tracks count the central three bases of 130-200-bp fragments as
nucleosome dyads. Total-H3 tracks represent paired-fragment coverage. ATAC
tracks represent Tn5-shifted cut sites. The R analysis uses 10-bp ATAC tracks,
while separate 1-bp tracks are retained for IGV. ATAC is interpreted as
accessibility rather than a direct measurement of nucleosome occupancy. Each
run is normalized to CPM independently; assays and biological replicates are
not pooled.

Run the complete resumable workflow with:

```sh
THREADS=6 caffeinate -dimsu \
  ./05_public_H3K4/scripts/reprocess_selected_nucleosome_occupancy.sh all
```

The final stage runs
[`scripts/analyze_reporter_loci_nucleosome.R`](scripts/analyze_reporter_loci_nucleosome.R).
Like the H3K4 script, it can be opened and sourced from RStudio. Required
packages and functions are declared before data processing, tidyverse is used
for tabular processing, figures are drawn first on the standard R graphics
device, and final plots are written as `tikzpicture` fragments with
`width = 7.5` and `lwdUnit = 72.27 / 96`.

```r
Sys.setenv(
  NUCLEOSOME_WORK_ROOT = "/path/to/issue69_nucleosome",
  NUCLEOSOME_OUTPUT_DIR = "/path/to/output"
)
source("05_public_H3K4/scripts/analyze_reporter_loci_nucleosome.R")
```

The tracked default output directory is
[`output/reporter_loci_nucleosome`](output/reporter_loci_nucleosome). Reporter
loci are ordered `pan-2`, `ad-3A`, `ad-3B`, `ad-8`, `his-3`, `mtr`, and
`csr-1`. Profile plots are split by assay, while heatmaps use within-run
genome-wide midrank percentiles with the same blue-white-orange scale and
bottom color bar as the H3K4 figures. The output also preserves the exact
command-line software versions recorded during preprocessing.

## Experimental metadata

The GEO records describe wild-type mycelial cultures grown in Vogel's medium
with 1.5% sucrose. The associated paper specifies 18 hours of shaking at 32
degrees C. Libraries were sequenced as nominal 75-bp single-end reads on an
Illumina NextSeq 500.

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

## Current author interpretation and decision (2026-09-07)

- No additional wet-lab experiment will be performed for this issue. The
  response will use this public-data reanalysis while retaining the limitations
  stated here.
- No prior ChIP-seq study of a `pan-2` mutant has been identified in the sources
  examined for this analysis. The chromatin state of the strain background used
  for the mutation assay is unknown. Differences between mutant and wild-type
  strains and between mycelia and conidia remain untested. Literature support
  for developmental-state differences has not yet been identified and verified
  for citation.
- In the selected wild-type mycelial datasets, `pan-2` does not show the
  consistently high genome-relative H3K4me1, H3K4me2, or H3K4me3 signal expected
  of an H3K4me-rich locus. Signal is not uniformly zero, and the
  gene-body-plus-or-minus-2-kb H3K4me3 percentiles are intermediate. This is
  therefore a bounded descriptive conclusion rather than evidence of complete
  absence.
- This pattern argues against the simplest model in which the observed
  mutation-frequency and indel-size changes depend directly on abundant
  H3K4me1/2/3 at the assayed `pan-2` locus. It does not exclude a locus-local
  effect. It is also compatible with indirect effects mediated through gene
  expression or with direct effects of larger-scale chromatin organization.
  The present data do not distinguish these possibilities.
- The complementary Kamei MNase/total-H3 and Ferraro ATAC reanalysis does not
  support describing the entire `pan-2` locus as nucleosome-free or extremely
  open. In both total-H3 replicates, the annotated gene body is high relative to
  other protein-coding genes (81.3-94.0th percentile), whereas the promoter is
  low (6.7-18.9th percentile). ATAC signal is intermediate at the promoter
  (55.0-59.3th percentile), low in the gene body (18.4-24.8th percentile), and
  higher when the plus-or-minus-2-kb flanks are included (66.9-75.0th
  percentile). MNase estimates vary among replicates. These data support
  substantial nucleosome occupancy across the gene body while leaving
  promoter-local depletion possible. They do not establish the state of the
  mutation-assay `pan-2` mutant background or dormant conidia. Detailed results
  and QC are in [`output/reporter_loci_nucleosome`](output/reporter_loci_nucleosome).

The RNA-seq runs SRR5177529 and SRR5177530 are not part of GSE154497. They come
from a separate study and are not matched to these ChIP-seq cultures, so the
ChIP-seq accession table does not include them.

## Selection boundary

GSE154497 is the best-supported public source identified for the wild-type
H3K4 browser tracks. A browser screenshot without track metadata cannot show
that the displayed track is byte-for-byte identical to the deposited bigWig.
The identification therefore records the public dataset and the
evidence-supported correspondence without treating an unlabelled screenshot
as the primary data source.

## Primary metadata sources

- [GSE154497 series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154497)
- [GSM4672245, H3K4me1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672245)
- [GSM4672246, H3K4me2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672246)
- [GSM4672247, H3K4me3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672247)
- [SRA Run Selector for SRP272099](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP272099)
- [Ferraro et al. (2021)](https://doi.org/10.1186/s12864-021-07774-0)
