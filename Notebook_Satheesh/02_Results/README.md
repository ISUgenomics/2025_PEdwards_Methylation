# 02_Results

This directory contains the strand-aware annotation results for differentially methylated loci (DMLs) and differentially methylated regions (DMRs) from the DSS analysis.

## Directory Structure

```
02_Results/
└── results/
    ├── DML_nearest_genes_*_strandAware.csv   # Annotated DML files
    ├── DMR_nearest_genes_*_strandAware.csv   # Annotated DMR files
    ├── DML_feature_summary_*.png             # DML feature distribution plots
    ├── DMR_feature_summary_*.png             # DMR feature distribution plots
```

## Output Files

### CSV Files (`*_strandAware.csv`)

Each CSV contains the original DML/DMR data plus strand-aware annotations:

| Column | Description |
|--------|-------------|
| `chr`, `pos`/`start`, `end` | Genomic coordinates |
| `gene_id`, `gene_name`, `gene_desc` | Nearest gene information |
| `distance` | Distance to nearest gene (bp) |
| `gene_strand` | Strand of nearest gene (+/-) |
| `TSS`, `TES` | Transcription start/end sites (strand-aware) |
| `direction` | Position relative to gene: upstream, downstream, or overlapping |
| `refined_feature` | Genomic feature: promoter, exon, intron, downstream, or intergenic |

### PNG Files (`*_feature_summary_*.png`)

Bar plots showing the distribution of DMLs/DMRs across genomic features.

## Feature Definitions

| Feature | Definition |
|---------|------------|
| **promoter** | Within 2kb upstream of TSS |
| **exon** | Overlaps an exon |
| **intron** | Within gene body, not overlapping exon |
| **downstream** | Within 2kb downstream of TES |
| **intergenic** | More than 2kb from any gene |

## Script

These results were generated using `00_Scripts/02_batch_strand_aware_annotation.R`.
