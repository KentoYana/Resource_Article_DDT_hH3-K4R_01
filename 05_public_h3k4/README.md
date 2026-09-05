# Public H3K4 ChIP-seq dataset

This directory records the public wild-type *Neurospora crassa* H3K4
ChIP-seq dataset selected for locus-level inspection. This metadata record does
not itself contain a reanalysis workflow or newly generated results.

## Dataset identification

The selected series is
[GSE154497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154497),
"Chromatin accessibility profiling in *Neurospora crassa* reveals features
associated with accessible and inaccessible chromatin." The series contains
wild-type mycelial ChIP-seq samples for H3K4me1, H3K4me2, and H3K4me3.

The H3K4me2 sample relevant to inspection of the `pan-2` locus is:

- GEO sample: [GSM4672246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4672246)
- SRA experiment: [SRX8737830](https://www.ncbi.nlm.nih.gov/sra/SRX8737830)
- SRA run: [SRR12229307](https://www.ncbi.nlm.nih.gov/sra/SRR12229307)

The series includes a public CPM-normalized bigWig for each of the three H3K4
marks. These deposited tracks are the preferred starting point for reproducing
a browser view because they retain the processing selected by the original
authors.

## Selection boundary

GSE154497 is selected as the best-supported public source for the wild-type
H3K4 browser tracks. A browser screenshot without track metadata is not, by
itself, sufficient to prove that the displayed track is byte-for-byte identical
to the deposited bigWig. This identification therefore records the public
dataset and the evidence-supported correspondence without treating an
unlabelled screenshot as the primary data source.
