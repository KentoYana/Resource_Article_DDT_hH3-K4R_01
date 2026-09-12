# Reporter-locus H3K4me descriptive summary

This table reports the nonduplicate tracks for `pan-2`. Percentiles are calculated separately for every run and region definition across all annotated protein-coding genes on the seven nuclear chromosomes.

| Study | Mark | Region | Mean CPM | Nonzero fraction | Genome percentile | Reporter-set rank |
|---|---|---|---:|---:|---:|---:|
| Ferraro et al. 2021 | H3K4me1 | promoter | 0.6654 | 0.733 | 30.1 | 7.0/7 |
| Ferraro et al. 2021 | H3K4me1 | gene_body | 0.9658 | 0.789 | 31.0 | 7.0/7 |
| Ferraro et al. 2021 | H3K4me1 | gene_body_plus_minus_2kb | 2.1893 | 0.888 | 45.4 | 7.0/7 |
| Ferraro et al. 2021 | H3K4me2 | promoter | 0.0000 | 0.000 | 4.8 | 7.0/7 |
| Ferraro et al. 2021 | H3K4me2 | gene_body | 0.3962 | 0.240 | 39.8 | 7.0/7 |
| Ferraro et al. 2021 | H3K4me2 | gene_body_plus_minus_2kb | 1.1399 | 0.278 | 38.5 | 7.0/7 |
| Sasaki et al. 2014 | H3K4me2 | promoter | 0.4212 | 0.742 | 26.0 | 7.0/7 |
| Sasaki et al. 2014 | H3K4me2 | gene_body | 0.4823 | 0.819 | 34.9 | 7.0/7 |
| Sasaki et al. 2014 | H3K4me2 | gene_body_plus_minus_2kb | 1.1134 | 0.850 | 40.8 | 6.0/7 |
| Ferraro et al. 2021 | H3K4me3 | promoter | 0.8735 | 0.133 | 46.5 | 5.0/7 |
| Ferraro et al. 2021 | H3K4me3 | gene_body | 0.2923 | 0.045 | 16.9 | 7.0/7 |
| Ferraro et al. 2021 | H3K4me3 | gene_body_plus_minus_2kb | 2.5649 | 0.191 | 63.2 | 5.0/7 |
| Storck et al. 2020 | H3K4me3 | promoter | 1.5552 | 1.000 | 37.2 | 4.0/7 |
| Storck et al. 2020 | H3K4me3 | gene_body | 1.3137 | 1.000 | 10.7 | 7.0/7 |
| Storck et al. 2020 | H3K4me3 | gene_body_plus_minus_2kb | 1.9490 | 1.000 | 52.5 | 5.0/7 |

## Interpretation boundaries

- These are descriptive CPM summaries, not statistical tests; no selected mark has biological replication within every study.
- Nonzero coverage is not equivalent to enrichment because no matched input is used.
- Studies and marks remain separate; absolute CPM values are not pooled or treated as directly exchangeable.
- Exact zero bases are counted separately from missing bases in the long-form TSV.
- `csr-1` is exploratory and is not treated as an established forward-mutation reporter.

## Current author interpretation and decision (2026-09-07)

- No additional wet-lab experiment will be performed for this issue; the response will use the available public-data reanalysis with the limitations below.
- No prior ChIP-seq study of a `pan-2` mutant has been identified in the sources examined for this analysis. The chromatin state of the strain background used for the mutation assay is therefore unknown. Possible differences between mutant and wild-type strains and between mycelia and conidia remain untested here; literature support for developmental-state differences has not yet been identified and verified for citation.
- In the selected wild-type mycelial datasets, `pan-2` does not show the consistently high genome-relative H3K4me1, H3K4me2, or H3K4me3 signal expected of an H3K4me-rich locus. Signal is not uniformly zero, and the gene-body-plus-or-minus-2-kb H3K4me3 percentiles are intermediate, so this is a bounded descriptive conclusion rather than evidence of complete absence.
- This pattern argues against the simplest model in which the observed mutation-frequency and indel-size changes depend directly on abundant H3K4me1/2/3 at the assayed `pan-2` locus. It does not exclude a locus-local effect and is also compatible with indirect effects mediated through gene expression or with direct effects of larger-scale chromatin organization; the present data do not distinguish these possibilities.
- The complementary Kamei MNase/total-H3 and Ferraro ATAC reanalysis does not support describing the entire `pan-2` locus as nucleosome-free or extremely open. In both total-H3 replicates, the annotated gene body is high relative to other protein-coding genes (81.3-94.0th percentile), whereas the promoter is low (6.7-18.9th percentile). ATAC signal is intermediate at the promoter (55.0-59.3th percentile), low in the gene body (18.4-24.8th percentile), and higher when the plus-or-minus-2-kb flanks are included (66.9-75.0th percentile); MNase estimates vary among replicates. These data therefore support substantial nucleosome occupancy across the gene body while leaving promoter-local depletion possible. They do not establish the state of the mutation-assay `pan-2` mutant background or dormant conidia; detailed results and QC are in `output/reporter_loci_nucleosome`.
