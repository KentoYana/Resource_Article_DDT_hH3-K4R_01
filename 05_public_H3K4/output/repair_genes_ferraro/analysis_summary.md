# Repair/DDT-gene H3K4me descriptive summary

This table reports within-run genome-wide midrank percentiles from the nonduplicate tracks. When more than one study is selected for a mark, values separated by a hyphen give the cross-study range; studies are not treated as biological replicates and are not pooled.

| Gene | H3K4me1 promoter | H3K4me1 gene body | H3K4me2 promoter | H3K4me2 gene body | H3K4me3 promoter | H3K4me3 gene body |
|---|---:|---:|---:|---:|---:|---:|
| `mus-9` | 84.8 | 45.8 | 79.6 | 58.9 | 16.0 | 53.2 |
| `uvs-2` | 97.3 | 92.9 | 91.9 | 78.2 | 48.3 | 34.4 |
| `mus-26` | 92.2 | 37.3 | 98.9 | 27.5 | 83.0 | 47.6 |
| `polh` | 78.2 | 67.9 | 78.6 | 62.8 | 48.4 | 41.5 |
| `qde-3` | 60.4 | 48.4 | 60.8 | 54.4 | 60.5 | 13.5 |
| `recQ2` | 63.6 | 53.5 | 35.6 | 40.0 | 7.7 | 17.2 |
| `mei-3` | 68.9 | 64.9 | 15.4 | 44.8 | 16.2 | 52.6 |
| `mus-11` | 98.6 | 72.6 | 87.4 | 58.7 | 86.5 | 53.5 |

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
