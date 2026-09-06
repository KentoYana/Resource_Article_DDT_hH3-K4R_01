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
