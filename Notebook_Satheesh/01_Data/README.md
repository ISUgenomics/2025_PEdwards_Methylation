Date: `2025-10-31`

workdir: `/work/gif4/satheesh/2025/12_PEdwards_Methylation/`

## Setting up data and scripts directories

```bash
# Create the scripts directory
mkdir -p 00_Scripts

# Create data directory
mkdir -p 01_Data
```

## Downloading the genome and gff3 file

workdir: `/work/gif4/satheesh/2025/12_PEdwards_Methylation/01_Data`

```bash
ssh novadtn

# Download the genome
wget -O crutilus_genome_softmasked.fa.gz https://ftp.ebi.ac.uk/pub/ensemblorganisms/Clethrionomys_rutilus/GCA_040207285.2/genome/softmasked.fa.gz

# Download the gff3 file
wget -O crutilus_genes.gff3.gz ftp://ftp.ebi.ac.uk/pub/ensemblorganisms/Clethrionomys_rutilus/GCA_040207285.2/ensembl/geneset/2025_03/genes.gff3.gz
```

Let us now get the stats of the genome gff file

```bash
# activate the methylation conda environment
conda activate methylation

# Install gffutils
mamba install bioconda::gffutils

# Unzip the gff3 file
unpigz crutilus_genes.gff3.gz

# Run the gffutils stats script
time python ../00_Scripts/01_run_gffutils.py
```

Output:

```bash
real    0m54.334s
user    0m43.519s
sys     0m0.927s
```

<details>
<summary>01_run_gffutils.py</summary>

```python
#!/usr/bin/env python

import gffutils

# Create a database
db = gffutils.create_db("crutilus_genes.gff3", dbfn="vole.db", force=True, merge_strategy="create_unique")

# Extract gene-level info
genes = []
for gene in db.features_of_type("gene"):
    genes.append({
        "gene_id": gene.id,
        "chrom": gene.chrom,
        "start": gene.start,
        "end": gene.end,
        "strand": gene.strand,
        "length": gene.end - gene.start + 1,
        "attributes": gene.attributes
    })

import pandas as pd
gene_df = pd.DataFrame(genes)
gene_df.to_csv("gene_annotation_summary.tsv", sep="\t", index=False)
```
</details>

**What does the above script do?**

The above script creates a database from the gff3 file and extracts gene-level info. The script then converts the gene-level info to a pandas DataFrame and saves it to a TSV file.

**Contents of the TSV file:**

```bash
head -n 11 gene_annotation_summary.tsv
```

*Output*

```bash
gene_id chrom   start   end     strand  length  attributes
gene:ENSBOLG00010002046 1       5075    6010    -       936     "ID: ['gene:ENSBOLG00010002046']
Name: ['OR4S2']
biotype: ['protein_coding']
description: ['olfactory receptor family 4 subfamily S member 2 [Ensembl NN prediction with score 100%]']
gene_id: ['ENSBOLG00010002046']
version: ['1']"
gene:ENSBOLG00010002049 1       210717  211190  -       474     "ID: ['gene:ENSBOLG00010002049']
Name: ['BTF3L4']
biotype: ['protein_coding']
description: ['basic transcription factor 3 like 4 [Ensembl NN prediction with score 93.37%]']
gene_id: ['ENSBOLG00010002049']
version: ['1']"
```

**Description of the output:**

This TSV file contains a simplified summary of all genes in the vole genome. Each row represents one gene with the following information:

- **gene_id**: Unique Ensembl gene identifier
- **chrom**: Chromosome number where the gene is located
- **start/end**: Genomic coordinates (base pair positions) marking the gene boundaries
- **strand**: DNA strand orientation (+ or -)
- **length**: Gene size in base pairs (calculated as end - start + 1)
- **attributes**: Additional gene metadata including:
  - Gene name (e.g., OR4S2, BTF3L4)
  - Biotype (e.g., protein_coding)
  - Functional description
  - Version number

This format should probably make it easy to filter, sort, and analyze gene features for downstream methylation analysis.

## Getting more inforamtion about the genes

workdir: `/work/gif4/satheesh/2025/12_PEdwards_Methylation/01_Data`

```bash 
conda activate methylation 

# activating ipython 
ipython
```

```ipython
import gffutils

# Load the database 
db = gffutils.FeatureDB("vole.db", keep_order=True)

# List all feature types 
for ftype in db.featuretypes():
    print(ftype)
```

*Output:*

```
CDS
C_gene_segment
J_gene_segment
V_gene_segment
Y_RNA
exon
five_prime_UTR
gene
lnc_RNA
mRNA
miRNA
ncRNA_gene
pseudogene
pseudogenic_transcript
rRNA
region
scRNA
snRNA
snoRNA
tRNA
three_prime_UTR
transcript
```

**Get different gene types (e.g. biotype or gene_class)**

Extracting biotype, gene_biotype, or gene_type attribute, if available in the database. 

```ipython
from collections import Counter

# Define what counts as a "gene-like" feature
gene_like_types = [
    "gene", "mRNA", "lnc_RNA", "miRNA", "ncRNA_gene",
    "pseudogene", "pseudogenic_transcript", "rRNA",
    "snRNA", "snoRNA", "scRNA", "tRNA",
    "C_gene_segment", "J_gene_segment", "V_gene_segment", "Y_RNA"
]

counts = Counter()

for ftype in db.featuretypes():
    if ftype in gene_like_types:
        counts[ftype] = db.count_features_of_type(ftype)

# Print summary
for k, v in counts.items():
    print(f"{k}\t{v}")
```

*Output:*

| Type | Count |
| ---- | ----- |
| C_gene_segment | 11 |
| J_gene_segment | 11 |
| V_gene_segment | 164 |
| Y_RNA | 36 |
| gene | 20257 |
| lnc_RNA | 969 |
| mRNA | 27212 |
| miRNA | 70 |
| ncRNA_gene | 3100 |
| pseudogene | 425 |
| pseudogenic_transcript | 425 |
| rRNA | 203 |
| scRNA | 50 |
| snRNA | 1021 |
| snoRNA | 857 |
| tRNA | 22 |