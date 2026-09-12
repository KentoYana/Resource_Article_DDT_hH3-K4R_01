# Repair/DDT-gene H3K4me descriptive summary

This table reports within-run genome-wide midrank percentiles from the nonduplicate tracks. Values separated by a hyphen give the range across the selected studies; studies are not treated as biological replicates and are not pooled.

| Gene | H3K4me1 promoter | H3K4me1 gene body | H3K4me2 promoter | H3K4me2 gene body | H3K4me3 promoter | H3K4me3 gene body |
|---|---:|---:|---:|---:|---:|---:|
| `mus-9` | 84.8 | 45.8 | 79.6-84.5 | 56.9-58.9 | 16.0-72.6 | 53.2-54.9 |
| `uvs-2` | 97.3 | 92.9 | 91.2-91.9 | 78.2-90.0 | 48.3-65.7 | 8.4-34.4 |
| `mus-26` | 92.2 | 37.3 | 94.5-98.9 | 27.5-46.8 | 83.0-88.0 | 43.2-47.6 |
| `polh` | 78.2 | 67.9 | 71.2-78.6 | 62.8-66.0 | 48.4-55.9 | 41.5-54.8 |
| `qde-3` | 60.4 | 48.4 | 57.2-60.8 | 54.4-57.7 | 34.1-60.5 | 13.5-29.1 |
| `recQ2` | 63.6 | 53.5 | 35.6-40.9 | 40.0-50.1 | 7.7-23.1 | 17.2-26.3 |
| `mei-3` | 68.9 | 64.9 | 15.4-55.9 | 44.8-54.0 | 16.2-38.6 | 49.9-52.6 |
| `mus-11` | 98.6 | 72.6 | 87.4-94.9 | 58.7-60.5 | 86.5-91.0 | 53.5-56.6 |

## Descriptive result

- The promoter profiles are heterogeneous rather than uniformly H3K4me-rich.
- `uvs-2`, `mus-26`, and `mus-11` have high H3K4me1 and H3K4me2 promoter percentiles; `mus-26` and `mus-11` also have high H3K4me3 promoter percentiles.
- `recQ2` has low H3K4me2 and H3K4me3 promoter percentiles, while `qde-3` and `mei-3` are not consistently high for H3K4me3.
- These basal profiles make a gene-specific indirect effect through altered transcription plausible, but do not demonstrate an expression change or causality.

## Interpretation boundaries

- These are descriptive CPM summaries and genome-relative percentiles, not statistical tests.
- No selected mark has biological replication within every study, and the studies remain separate.
- Nonzero coverage is not equivalent to enrichment because no matched input is used.
- Basal wild-type H3K4me profiles can show whether these genes lie in H3K4me-rich chromatin under the public-data conditions, but cannot determine transcriptional effects of `hH3-K4R` or changes after UV irradiation or replication stress.
- Exact zero bases are retained and counted separately from missing bases in the long-form TSV.
