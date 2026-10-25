# Race-tube hyphal growth (Figure S1e)

This workflow provides the measured distances and descriptive growth curves for
wild type, hH3-K4R, and set-1 KO. It reproduces the numerical quantities plotted
in Figure S1e. It performs no hypothesis tests or fitted growth models.

This workflow is an addition after v1.1.6; it is not contained in that release.

## Input and provenance

`dataset/racetube210111.csv` is the author-supplied CSV, preserved byte for byte,
including its empty trailing column. The empty column is not used in calculations.
It contains 234 records: three groups, three tube identifiers per group, and 26
time points including the zero-distance starting time (2021-01-11 21:00:00 Japan
time). Measurements end on 2021-01-20 at 20:17:00 (215.283333333333 hours).

| CSV identifier prefix | Plot label | Tube suffixes |
| --- | --- | --- |
| C1-T10-28a | wild-type | -1, -2, -3 |
| KNT-Op3-87a | hH3-K4R | -1, -2, -3 |
| FGSC15828 | set-1 KO | -1, -2, -3 |

The mapping is retained from the author's `analysis/script230525.R`.
Suffixes identify tube trajectories; repeated time points are not independent
replicates. No additional strain identity or independence is inferred from the
file names.

The two supplied missing values are retained:

- KNT-Op3-87a-3 at 2021-01-14 09:00:00.
- FGSC15828-2 at 2021-01-13 20:00:00.

For those two group/time combinations, n = 2; otherwise n = 3. No interpolation,
imputation, additional exclusion, or conversion of missing values to zero is used.

The supplied CSV matches the earlier per-group WT, K4R, and set1 text data at
every time point, including missing values. In the WT text file, Tube1 and Tube2
correspond to CSV suffixes -2 and -1, respectively; Tube3 corresponds to -3.
The K4R and set1 suffixes match the earlier tube order. This ordering difference
does not change group means or standard errors. The CSV identifiers are retained.

Input SHA-256: `98675ebfa3c7e18534a57bbc407b76b8162cdec59746a1d8a688acc895227413`.

## Calculation and adaptation

`script_race_tube.R` is adapted from the supplied `script230525.R` (whose internal
header says `script221018a.R`). The original script SHA-256 is `1a855b413f090b2e82badb7b35fc26dbb39242fc3cfd53227d51f3d2c50974ab`.
The local source's absolute paths and interactive-session files are not included.

The existing calculation is preserved:

- Parse timestamps in the Japan time zone and subtract the earliest timestamp.
- Use measured cumulative distance directly, without subtracting the first
  nonzero measurement, normalization, or transformation of zero/one values.
- At each genotype/time point, calculate the arithmetic mean and sample standard
  deviation divided by the square root of the number of nonmissing observations.
- Retain full precision during computation, with no early significant-digit rounding.
- Retain the original endpoint calculation: each final distance divided by the
  elapsed time, followed by its group mean and SE. This is a descriptive average
  rate from the origin, not a fitted slope or an inferential test.

Only execution and output handling have been adapted: paths are resolved relative
to the script, unused dependencies and interactive plotting are removed, numerical
outputs and sample counts are saved, and PDF/TikZ export is enabled. The unit label
is rendered as `mm` in both formats. Group order, shapes, colors, error-bar widths,
and the plotting aspect ratio are retained. No random operation or seed is used.

The old 2021 `RaceTube.R` is not this workflow: it skipped the commented origin,
subtracted the first nonzero measurement, and rounded summaries. It is not used to
reproduce Figure S1e. The supplied 2023/2024 TeX exports are also not adopted as the
current manuscript figure.

## Run

Requires R, tidyverse, tikzDevice, and a working LaTeX installation for tikzDevice's
font metrics. No saved R workspace, user profile, or network access is required.

From the repository root:

```sh
Rscript --vanilla 07_race_tube/script_race_tube.R
```

The script can also run from another working directory using its absolute path.
An optional first argument selects an output directory for an isolated rerun:

```sh
Rscript --vanilla 07_race_tube/script_race_tube.R /tmp/race-tube-check
```

The standalone PDF is written directly by R. The TikZ file is an editable fragment
for inclusion within a LaTeX document using `tikz`. It is not a standalone document.

## Outputs

- `observations.csv`: all tube/time records, measured distances, elapsed hours and
  genotype mapping; missing values remain `NA`.
- `growth_summary.csv`: 78 group/time means, SEs, nonmissing n and missing counts.
- `endpoint_velocity.csv`: the original descriptive endpoint rates in mm/hour.
- `race_tube.pdf`: a standalone preview of the regenerated resource plot.
- `race_tube.tex`: an editable TikZ plot with a deterministic generator header.
- `R-session.txt`: R and package versions for this regeneration, not a claim about
  the software versions used for the historical manuscript export.

The original manuscript has a curated TikZ presentation copy (generated in July
2025). Its bounding region and typography can differ from a fresh resource export.
This workflow preserves its numerical content and does not replace that copy.

## Verification (2026-10-26)

- The supplied CSV and the repository copy have identical SHA-256 hashes.
- The original script was run in an isolated directory with only its working path
  redirected and its interactive plot display suppressed. All observations,
  group means/SEs and endpoint rates agree exactly with the adapted script's
  numerical CSV fields. The adaptation adds n and missing-count columns only.
- Independent arithmetic from the CSV agrees with all 78 mean/SE pairs within
  3.5e-13 mm (floating-point/CSV serialization precision).
- All 78 mean positions and SE bars agree with the existing Figure S1e TikZ
  coordinates within 0.009 pt, below its 0.01-pt coordinate precision. This check
  uses the existing plot's rounded axis ticks, so it is a coordinate-precision
  comparison, not a claim of byte-identical figure files.
- A separate copy of this workflow was run with `Rscript --vanilla` from an
  unrelated working directory, without `.RData` or history files. Its three CSV
  outputs and TikZ source are byte-identical. The PDF page content is identical;
  PDF creation-time metadata can vary between runs.
- The standalone PDF and a separately compiled TikZ preview passed `qpdf --check`.
  Both rendered views were inspected for axes, units, group mapping, shapes,
  colors, error bars, legend readability and clipping. No manuscript figure was
  replaced. The original point sizes and overlapping trajectories are retained.

The regeneration used R 4.6.1, tidyverse 2.0.0, dplyr 1.2.1, ggplot2 4.0.3,
stringr 1.6.0 and tikzDevice 0.12.6. See `output/R-session.txt` for the complete
session. These are reproduction-environment versions, not historical versions.

## License

The script uses the repository MIT License. Original data, outputs and this
README use CC BY 4.0, as specified in the repository root.
