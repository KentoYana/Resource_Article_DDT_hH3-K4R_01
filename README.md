# Resource_Article_DDT_hH3-K4R_01

This repository contains datasets, analysis scripts, and generated output files associated with a manuscript on the role of histone H3 lysine 4 methylation in DNA damage tolerance in *Neurospora crassa*.

The repository is intended to support transparency and reproducibility of the analyses reported in the manuscript. Each analysis directory contains the input data, analysis script, and output files used for the corresponding part of the study.

## Repository structure

```text
.
├── 01_killing_test
├── 02_reversion_assay
├── 03_seq_analysis
├── 04_quantitative_spot_test
├── 91_cross_information
└── 99_pan-2_identification
```

## Contents

### `01_killing_test`

This directory contains the dataset and R script for the UV killing-test analysis.

Main files:

- `dataset/killing_data.csv`: raw killing-test data.
- `dataset/strain_list.csv`: strain and genotype annotation table.
- `script_killing_test.R`: R script for data processing, beta-regression analysis, post hoc comparisons, and figure generation.
- `output/`: processed data, statistical summaries, comparison results, and TikZ figure output.

Run from the repository root:

```sh
Rscript 01_killing_test/script_killing_test.R
```

### `02_reversion_assay`

This directory contains the dataset and R script for the reversion-frequency analysis.

Main files:

- `dataset/reversion_data.csv`: raw reversion-assay data.
- `dataset/strain_list.csv`: strain annotation table.
- `dataset/genotype_list.csv`: genotype annotation table.
- `script_reversion_frequency_AUC.R`: R script for reversion-frequency calculation, AUC-based analysis, pairwise comparisons, and figure generation.
- `output/`: processed data, AUC data, pairwise comparison results, Welch's t-test summaries, and TikZ figure output.

Run from the repository root:

```sh
Rscript 02_reversion_assay/script_reversion_frequency_AUC.R
```

### `03_seq_analysis`

This directory contains the dataset and scripts for sequence-based analysis of revertants.

Main files:

- `dataset/reversion_sequence.csv`: revertant sequence information.
- `dataset/pan-2_genomic.fasta`: reference genomic sequence for `pan-2`.
- `dataset/palindrome.txt`: sequence information used for palindrome-related inspection.
- `dataset/ir_pairs.csv`: inverted-repeat pair information used for indel classification.
- `dataset/strain_list.csv`: strain annotation table.
- `script_rev_seq_analysis.R`: R script for indel extraction, positional summaries, indel-size analysis, negative-binomial GLM analysis, and figure generation.
- `run_clustalo_B36.zsh`: shell script for generating a Clustal-format alignment using Clustal Omega.
- `output/`: processed sequence data, indel summaries, statistical results, FASTA files, alignment files, and TikZ figure output.

Run the R analysis from the repository root:

```sh
Rscript 03_seq_analysis/script_rev_seq_analysis.R
```

Run the Clustal Omega alignment script from the repository root:

```sh
zsh 03_seq_analysis/run_clustalo_B36.zsh
```

### `04_quantitative_spot_test`

This directory contains the datasets and R script for quantitative spot-test analysis.

Main files:

- `dataset/exp_list.csv`: list of quantitative spot-test experiments.
- `dataset/qSpot-*`: SpotScanner output files and strain annotation files for each target background.
- `script_qSpot_test.R`: R script for importing SpotScanner results, calculating normalized spot-coverage ratios, beta-regression analysis, likelihood-ratio tests, slope comparisons, dose-wise comparisons, and figure generation.
- `output/`: processed data, statistical summaries, target-specific output directories, and TikZ figure output.

The target-specific dataset directories include:

- `qSpot-mei-3`
- `qSpot-mus-11`
- `qSpot-mus-26-polh`
- `qSpot-mus-9`
- `qSpot-recQ`
- `qSpot-uvs-2`

Run from the repository root:

```sh
Rscript 04_quantitative_spot_test/script_qSpot_test.R
```

### `91_cross_information`

This directory contains the dataset and R script for cross-information or segregation analysis related to `hH3-K4R`.

Main files:

- `dataset/cross_data.csv`: cross-information dataset.
- `script_segregation_hH3-K4R.R`: R script for processing cross data and generating summary output.
- `output/`: processed data and summary tables.

Run from the repository root:

```sh
Rscript 91_cross_information/script_segregation_hH3-K4R.R
```

### `99_pan-2_identification`

This directory contains sequence files and a shell script used for `pan-2` identification and related sequence alignment.

Main files:

- `dataset/pan-2_genomic.fasta`: genomic sequence file for `pan-2`.
- `dataset/PTHR20881.fasta`: protein-family sequence file used for comparative alignment.
- `run_clustalo_pan-2_analysis.zsh`: shell script for generating Clustal-format alignments using Clustal Omega.
- `output/`: Clustal-format alignment files.

Run from the repository root:

```sh
zsh 99_pan-2_identification/run_clustalo_pan-2_analysis.zsh
```

## Requirements

The analyses were performed using R scripts and shell scripts.

### R

The following R packages are required:

```r
install.packages(c(
  "tidyverse",
  "ggfortify",
  "betareg",
  "emmeans",
  "patchwork",
  "RColorBrewer",
  "ggh4x",
  "MASS",
  "seqinr",
  "tikzDevice",
  "here"
))
```

The scripts use the `here` package to locate files relative to the repository root. For reproducibility, run the scripts from the repository root.

### Clustal Omega

The shell scripts require Clustal Omega.

On macOS with Homebrew:

```sh
brew install clustal-omega
```

On Debian/Ubuntu-based systems:

```sh
sudo apt-get install clustalo
```

Check that the command is available:

```sh
clustalo --version
```

## Reproducing the analyses

After cloning the repository and installing the required software, the main analyses can be rerun from the repository root as follows:

```sh
Rscript 01_killing_test/script_killing_test.R
Rscript 02_reversion_assay/script_reversion_frequency_AUC.R
Rscript 03_seq_analysis/script_rev_seq_analysis.R
zsh 03_seq_analysis/run_clustalo_B36.zsh
Rscript 04_quantitative_spot_test/script_qSpot_test.R
Rscript 91_cross_information/script_segregation_hH3-K4R.R
zsh 99_pan-2_identification/run_clustalo_pan-2_analysis.zsh
```

Each script reads input files from the corresponding `dataset/` directory and writes generated files to the corresponding `output/` directory.

## Output files

The `output/` directories contain processed data tables, statistical summaries, model outputs, and figure source files.

Common output file types include:

- `.csv`: processed data, summary tables, model outputs, and pairwise comparison results.
- `.txt`: human-readable summaries of model fits and statistical tests.
- `.tex`: TikZ/PGFPlots figure source files.
- `.fasta`: sequence files generated during sequence analysis.
- `.aln`: Clustal-format alignment files.

## Notes on TikZ figure files

Some plots are exported as TikZ/PGFPlots `.tex` files using the R package `tikzDevice`.

The TikZ files are intended as editable figure source files. In particular, axis labels, legend labels, panel labels, text placement, and other label-level formatting may be adjusted by directly editing the generated TeX source files before inclusion in the manuscript or figure assembly workflow.

Rerunning the analysis scripts may regenerate the TikZ files and overwrite manual TeX-level edits. If manual label edits are made, keep a separate copy of the edited TeX files or reapply the edits after regeneration.

## Repository hygiene

The repository includes a `.gitignore` file to exclude system-generated, session-specific, and user-specific files that are not required for reproducibility.

Examples include macOS metadata files, R history files, R session files, RStudio project-user files, R Markdown cache files, and local environment files.

## License

This repository uses separate licenses for code and research data.

- Analysis scripts and shell scripts are licensed under the MIT License.
- Datasets, processed output files, alignment files, figure source files, and documentation are licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0), unless otherwise noted.

See the license files in this repository for the full license terms.

## Citation

If you use the datasets, analysis scripts, or generated output files from this repository, please cite the associated manuscript.

Full citation information should be added here after publication or preprint release.

## Contact

For questions about the datasets, analysis scripts, or generated output files, please contact the corresponding author of the associated manuscript.

Questions and comments may also be submitted through GitHub Issues in this repository. Kento Momma may also be contacted through the GitHub profile associated with this repository or through the ORCID record below.

- GitHub: [KentoYana](https://github.com/KentoYana)
- ORCID: [Kento Momma](https://orcid.org/0009-0004-7481-1364)
