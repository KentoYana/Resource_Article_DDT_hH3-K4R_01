# Reporter-locus nucleosome and accessibility summary

This table reports nonduplicate tracks for `pan-2`. Percentiles are calculated separately for every run, assay, and region definition across all annotated protein-coding genes on the seven nuclear chromosomes.

| Study | Assay | Region | Mean CPM | Nonzero fraction | Genome percentile | Reporter-set rank |
|---|---|---|---:|---:|---:|---:|
| Kamei et al. 2021 rep1 | MNase | promoter | 0.0346 | 0.025 | 49.7 | 5.0/7 |
| Kamei et al. 2021 rep1 | MNase | gene_body | 0.0607 | 0.040 | 76.7 | 2.0/7 |
| Kamei et al. 2021 rep1 | MNase | gene_body_plus_minus_2kb | 0.0430 | 0.030 | 54.1 | 6.0/7 |
| Kamei et al. 2021 rep2 | MNase | promoter | 0.0388 | 0.040 | 54.9 | 6.0/7 |
| Kamei et al. 2021 rep2 | MNase | gene_body | 0.0539 | 0.052 | 67.7 | 4.0/7 |
| Kamei et al. 2021 rep2 | MNase | gene_body_plus_minus_2kb | 0.0420 | 0.043 | 50.3 | 4.0/7 |
| Kamei et al. 2021 rep3 | MNase | promoter | 0.0202 | 0.028 | 25.1 | 6.5/7 |
| Kamei et al. 2021 rep3 | MNase | gene_body | 0.0462 | 0.052 | 40.8 | 4.0/7 |
| Kamei et al. 2021 rep3 | MNase | gene_body_plus_minus_2kb | 0.0336 | 0.042 | 19.9 | 6.0/7 |
| Kamei et al. 2021 rep1 | total H3 | promoter | 2.7335 | 0.892 | 6.7 | 6.0/7 |
| Kamei et al. 2021 rep1 | total H3 | gene_body | 8.8289 | 1.000 | 81.3 | 1.0/7 |
| Kamei et al. 2021 rep1 | total H3 | gene_body_plus_minus_2kb | 5.7860 | 0.955 | 17.2 | 5.0/7 |
| Kamei et al. 2021 rep2 | total H3 | promoter | 4.1250 | 1.000 | 18.9 | 4.0/7 |
| Kamei et al. 2021 rep2 | total H3 | gene_body | 9.1005 | 1.000 | 94.0 | 1.0/7 |
| Kamei et al. 2021 rep2 | total H3 | gene_body_plus_minus_2kb | 6.1353 | 1.000 | 23.4 | 3.0/7 |
| Ferraro et al. 2021 rep1 | ATAC | promoter | 0.4288 | 0.868 | 59.3 | 6.0/7 |
| Ferraro et al. 2021 rep1 | ATAC | gene_body | 0.0962 | 0.716 | 18.4 | 7.0/7 |
| Ferraro et al. 2021 rep1 | ATAC | gene_body_plus_minus_2kb | 0.3118 | 0.846 | 75.0 | 4.0/7 |
| Ferraro et al. 2021 rep2 | ATAC | promoter | 0.3083 | 0.959 | 55.0 | 6.0/7 |
| Ferraro et al. 2021 rep2 | ATAC | gene_body | 0.1684 | 0.894 | 24.8 | 7.0/7 |
| Ferraro et al. 2021 rep2 | ATAC | gene_body_plus_minus_2kb | 0.2606 | 0.941 | 66.9 | 4.0/7 |

## Interpretation boundaries

- These are descriptive CPM summaries, not statistical tests. Biological replicates remain visible and are not pooled.
- MNase dyad density and total-H3 coverage are complementary proxies for nucleosome occupancy; neither alone uniquely measures absolute occupancy.
- ATAC cut-site density measures accessibility and is not treated as a direct nucleosome-occupancy measurement.
- Assays remain separate; absolute CPM values are not pooled or treated as directly exchangeable.
- Exact zero bases are counted separately from missing bases in the long-form TSV.
- `csr-1` is exploratory and is not treated as an established forward-mutation reporter.
- Genome-wide nonduplicate replicate Spearman correlations range from 0.89-0.97 for ATAC, 0.64-0.78 for MNase, and 0.50-0.65 for total H3 across the three region definitions. The moderate total-H3 agreement and shallow rep1 reduce the precision of locus-level inference.
- MNase pairs in the 130-200-bp analysis window account for 12.0%, 15.0%, 22.5% of retained pairs in replicates 1-3, respectively.

## Interpretation status

- In the two total-H3 replicates, the `pan-2` gene body is high relative to other protein-coding genes (81.3-94.0th percentile), whereas its promoter is low (6.7-18.9th percentile).
- `pan-2` ATAC accessibility is intermediate at the promoter (55.0-59.3th percentile), low in the gene body (18.4-24.8th percentile), and higher only when the plus-or-minus-2-kb flanks are included (66.9-75.0th percentile).
- MNase dyad percentiles vary among the three replicates: promoter 25.1-54.9, gene body 40.8-76.7, and gene body plus or minus 2 kb 19.9-54.1. This variation is retained rather than averaged away.
- Taken together, these public wild-type data do not support describing the entire `pan-2` locus as nucleosome-free or extremely open. They are compatible with relatively low total-H3 signal at the promoter but substantial nucleosome occupancy across the annotated gene body.
- Wild-type public data do not establish the chromatin state of the mutation-assay `pan-2` mutant background or dormant conidia.
- No additional wet-lab experiment is planned for this issue; any conclusion will remain limited to the public-data conditions analyzed here.
