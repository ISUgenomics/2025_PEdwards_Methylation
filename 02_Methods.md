# Methods

Written summary of methods performed in this repo. This is the methods write up for the paper.

## Sample Collection

## Reference Genome and Annotation

The northern red-backed vole (*Clethrionomys rutilus*) reference genome assembly mMyoRut1_p1.0 (GenBank accession GCA_040207285.2) was obtained from Ensembl. The genome consists of 28 chromosomes plus the mitochondrial genome and 55 unplaced scaffolds. Gene annotations were retrieved from the Ensembl Rapid Release database (version 2025_03), comprising 20,257 protein-coding genes, 3,100 non-coding RNA genes, 969 long non-coding RNAs, and 425 pseudogenes. Genome annotation statistics were extracted using gffutils (v0.12) in Python.

## Genomic Feature Annotation

DMLs and DMRs were annotated to genomic features using a custom strand-aware annotation pipeline implemented in R using the rtracklayer, GenomicRanges, and GenomicFeatures Bioconductor packages. For each differentially methylated site or region, the nearest gene was identified using the `distanceToNearest()` function from GenomicRanges. Strand-aware annotation was performed by computing the transcription start site (TSS) and transcription end site (TES) based on gene strand orientation.

Genomic features were classified as follows:
- **Promoter**: Within 2 kb upstream of the TSS
- **Exon**: Overlapping an annotated exon
- **Intron**: Within the gene body but not overlapping an exon
- **Downstream**: Within 2 kb downstream of the TES
- **Intergenic**: More than 2 kb from any annotated gene

Exon/intron discrimination was performed by overlapping DML/DMR coordinates with exon annotations from the GFF3 file using `findOverlaps()` from GenomicRanges. The direction of each site relative to its nearest gene (upstream, downstream, or overlapping) was determined in a strand-aware manner.

## Software and Packages

The following software and R packages were used in this analysis:
- **Python**: gffutils for GFF3 parsing
- **R/Bioconductor**: rtracklayer (GFF3 import), GenomicRanges and GenomicFeatures (genomic interval operations), ChIPseeker (annotation), txdbmaker (TxDb construction)
- **R/CRAN**: dplyr and tidyr (data manipulation), ggplot2 (visualization). 