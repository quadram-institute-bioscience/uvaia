## Reference database for testing uvaia

The file [03.unique_acgt.aln.xz](03.unique_acgt.aln.xz) contains a small database with 9185 aligned sequences from the
UK sequenced at the [QIB](https://quadram.ac.uk/).
This file can be used as the reference to test `uvaia` against your aligned query sequences.
The sequences were aligned with `uvaialign` against the reference genome [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) without trimming. 

Their PANGO lineages as estimated by [pangolin v4.0.6](https://github.com/cov-lineages/pangolin) can be found in file 
[03.unique_acgt.lineage_report.csv](03.unique_acgt.lineage_report.csv).
As the name suggests these sequences are a random sample of the full database such that every sequence has a unique
frequency of ACGT.

For a [COGUK alignment](https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/), or if
you have access, to the unaligned sequences at [GISAID](https://www.gisaid.org/) (we suggest `uvaialign` to align the sequences :wink:).

## Sequences used in manuscript

The files [04.sample_1_1k.names](04.sample_1_1k.names) and [04.sample_3_5k.names](04.sample_3_5k.names) contain the
sequences names of the samples used in the manuscript, and can all be found in the reference alignment above. 

For the timing analyses, we used [the COGUK data set](https://www.cogconsortium.uk/) ("Unmasked alignment" available from [the archived version of 7th May 2023](https://webarchive.nationalarchives.gov.uk/ukgwa/20230507102210/https://www.cogconsortium.uk/priority-areas/data-linkage-analysis/public-data-analysis/))
