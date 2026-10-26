# v1.1.7

Resource package for the revised manuscript by Yanagisawa and colleagues,
"Histone H3K4 methylation-related chromatin regulation modulates UV-induced
mutagenic outcomes and RecQ-associated DNA damage tolerance in Neurospora crassa".

## Additions since v1.1.6

- `07_race_tube`: original distance measurements and a portable workflow for
  Figure S1e, with observations, growth summaries, endpoint velocities, PDF/TikZ
  outputs and the recorded R environment. It reproduces the previously used
  descriptive calculations without introducing hypothesis tests.
- `08_source_images`: 20 original images underlying Figures 1a/b and S1b/d,
  image checksums, and concise figure-to-file and sample correspondence.
- `08_source_images/qspot_representatives.md`: correspondence for Figures 4a
  and 5a, linking 24 original qSpot images and their SpotScanner CSVs to the
  individual strain rows in Table S4. These image/CSV files already existed in
  v1.1.6; the guide does not change their contents or measurements.
- Updated manuscript title and data documentation. The conidial-development
  photograph uses the descriptive filename `conidiation_slants.JPG`.

## Preservation and reproducibility

The earlier analysis directories retain their inputs, scripts and outputs.
The 79 statistical CSV/TXT files reproduced in Supplementary File S2 retain
exactly the same bytes as in v1.1.6, including nonsignificant results. No new
inferential analysis or reinterpretation is introduced by this version update.

Reproduction instructions and limitations are documented in the repository
README and individual workflow READMEs. Analysis scripts use the MIT License;
original data, outputs and documentation use CC BY 4.0. Third-party materials
retain their original terms.
