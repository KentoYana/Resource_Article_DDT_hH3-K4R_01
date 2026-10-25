# Original images for Figures 1 and S1

This directory contains 20 original image files corresponding to Figure 1a/b
and Figure S1b/d. Files are preserved byte for byte as supplied by the authors.
`SHA256SUMS` records their checksums. These files were added after v1.1.6 and
are not included in that release.

## Figure-to-file correspondence

| Figure | Files | Capture date (EXIF) | Subjects shown in the manuscript |
| --- | --- | --- | --- |
| 1a | `figure1/Spot230119/*.JPG` (9 files) | 2023-01-19 | wild type and hH3-K4R; upper two spot rows |
| 1b | `figure1/Spot230126/*.JPG` (9 files) | 2023-01-26 | wild type and delta set-1; upper two spot rows |
| S1b | `figureS1/Band211214a.jpg` | 2021-12-14 | marker, C1-T10-28a and KNT-Op3-87a |
| S1d | `figureS1/P4173247.JPG` | 2024-04-17 | wild type, hH3-K4R and delta set-1; left to right |

Capture dates are the dates recorded by the cameras, not inferred inoculation
or experiment-start dates. Figure labels refer to the manuscript's final-revision
figure numbering.

## Figure 1a/b: qualitative spot assays

The two experiment sets have different row orders. Read rows from top to bottom
in the stored photographs; conidial concentration decreases from left to right.
The author-provided experiment records identify the incubation time as two days.

| Image set | Row 1 | Row 2 | Row 3 | Row 4 | Rows used in Figure 1 |
| --- | --- | --- | --- | --- | --- |
| Spot230119 | wild type | hH3-K4R | hH3-K4A | delta set-1 | 1 and 2 (panel a) |
| Spot230126 | wild type | delta set-1 | hH3-K4R | hH3-K4A | 1 and 2 (panel b) |

The full photographs retain all four rows. The following filenames occur in both
sets; treatment conditions are those stated in the Figure 1 legend.

| Filename | Treatment | Concentration or UV dose |
| --- | --- | --- |
| `Control.JPG` | untreated control | none |
| `BLM.JPG` | bleomycin | 1.0 micromolar |
| `CPT.JPG` | camptothecin | 0.57 micromolar |
| `UV.JPG` | UV irradiation | 300 J/m^2 |
| `4NQO.JPG` | 4-nitroquinoline 1-oxide | 0.13 micromolar |
| `CDDP.JPG` | cisplatin | 0.10 millimolar |
| `MMS.JPG` | methyl methanesulfonate | 2.4 millimolar |
| `MNNG.JPG` | N-methyl-N'-nitro-N-nitrosoguanidine | 7.1 micromolar |
| `HU.JPG` | hydroxyurea | 15 millimolar |

Only the conditions displayed in Figure 1a/b are included in this image set.
These are the representative images underlying the panels, not a collection of
all experimental repeats. Interpretation follows the manuscript legend.

## Figure S1b: genomic PCR gel

`Band211214a.jpg` is the full gel photograph. To read it in the manuscript's
orientation, turn the photograph through 180 degrees so that the wells are at
the top. The lanes then run from left to right as follows:

1. lambda/HindIII marker.
2. C1-T10-28a (wild type).
3. KNT-Op3-87a (hH3-K4R).
4. Additional samples (lanes 4 and 5), not used in Figure S1b.

Figure S1b shows the first three contiguous lanes in this orientation.
The PCR primers were hH3-For-TA and hH3-Rev-TA. The matching sequence-validation
files are in `../06_hH3_sequence_validation/`.

## Figure S1d: conidial development

`P4173247.JPG` shows all three slants in one photograph. From left to right:
wild type, hH3-K4R and delta set-1. The manuscript conditions are 10^5 conidia
inoculated onto glycerol complete medium slants, incubated at 25 degrees C for
seven days. The three slants correspond to the three images in Figure S1d.

## Integrity and license

From this directory, verify the image files with:

```sh
shasum -a 256 -c SHA256SUMS
```

Original images and this documentation use the repository's CC BY 4.0 license.
