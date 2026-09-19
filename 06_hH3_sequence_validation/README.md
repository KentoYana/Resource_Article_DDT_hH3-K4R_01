# hH3 K4 sequence validation

This directory contains the raw Sanger chromatograms and base-call exports used
to validate the H3 K4 codon shown in Supplementary Figure S1.

## Contents

- `dataset/`: four original AB1 chromatograms and their matching PHD exports
  for C1-T10-28a (wild type) and KNT-Op3-87a (`hH3-K4R`), sequenced in both
  directions.
- `script_hH3_sequence_validation.R`: parses the ABIF tags directly, verifies
  that each PHD export matches the corresponding AB1 sequence, quality scores,
  and peak positions, checks the AAG/AGA K4 codons, and writes the trace panel
  as editable TikZ.
- `output/`: verification tables and metadata, the R session record, and the
  generated TikZ trace panel.

Run from the repository root:

```sh
Rscript 06_hH3_sequence_validation/script_hH3_sequence_validation.R
```

The script can also be opened and sourced from RStudio with the repository root
as the project working directory.

## Interpretation limits

The four saved AB1/PHD pairs confirm AAG in C1-T10-28a and AGA in
KNT-Op3-87a in both sequencing directions. The displayed traces are from the
forward reads; signal heights are scaled independently and are not allele
fractions. These files do not independently document strain 84A, the reported
full-length checks, or nuclear purity. The manuscript interprets homokaryotic
status using the combined molecular validation, stable vegetative propagation,
and transmission through crosses.

## Raw metadata exception

The original AB1 instrument metadata include the plate identifier
`Plate_20211216_120326_yanagisawa`. The author explicitly approved publication
of the unmodified raw files as an exception; the embedded surname has therefore
not been removed or rewritten.
