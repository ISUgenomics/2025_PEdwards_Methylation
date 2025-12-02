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
head -n 13 gene_annotation_summary.tsv
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

## Getting more information about the genes

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

## DMLs and DMRs

```bash
 ltree 01_Data/dmls_dmrs_dss/
```

```
├── DMLs_LowDecline.csv
├── DMLs_PeakDecline.csv
├── DMLs_PeakLow.csv
├── DMRs_LowDecline.csv
├── DMRs_PeakDecline.csv
└── DMRs_PeakLow.csv
```

## Convert DMLs and DMRs to BED format

```bash
# DMLs
awk -F',' 'NR>1 {print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5}' 01_Data/dmls_dmrs_dss/DMLs_LowDecline.csv > 01_Data/dmls_dmrs_dss/DMLs_LowDecline.bed

# DMRs
awk -F',' 'NR>1 {print $1"\t"$2"\t"$3"\t"$6}' 01_Data/dmls_dmrs_dss/DMRs_LowDecline.csv > 01_Data/dmls_dmrs_dss/DMRs_LowDecline.bed
```

## Prepare annotation (genes/TSS)

```bash
# Extract genes from GFF3
awk -F '\t' '$3=="gene"{print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' 01_Data/crutilus_genes.gff3 > 01_Data/crutilus_genes.bed

# Extract TSS from GFF3
awk -F '\t' '$3=="gene"{

  if($7=="+")
    print $1"\t"$4-1"\t"$4"\t"$9"\t.\t"$7;
  else
    print $1"\t"$5-1"\t"$5"\t"$9"\t.\t"$7
}' 01_Data/crutilus_genes.gff3 > 01_Data/crutilus_tss.bed
```

## R

```R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.20") # with R 4.4

BiocManager::install(c("Biostrings", "GenomeInfoDb", "GenomicRanges"))
BiocManager::install("genomation")
```

All packages are installed.

```R
library(genomation)
library(GenomicRanges)

# Read genes (or TSSs)
genes <- readTranscriptFeatures("genes.bed")

# Read DML/DMR as GRanges
dml <- readGeneric("DMLs.bed", header=FALSE)
dmr <- readGeneric("DMRs.bed", header=FALSE)

# Compute densities around gene body or TSS
mat_dml <- ScoreMatrixBin(target=dml, windows=genes, bin.num=100)
mat_dmr <- ScoreMatrixBin(target=dmr, windows=genes, bin.num=100)

# Plot average profiles
plotMeta(mat_dml, main="DML density across genes")
plotMeta(mat_dmr, main="DMR density across genes")
```

### Download the genome from NCBI

```bash
 /work/gif3/satheesh/programs/datasets download genome accession GCA_040207285.2 --include genome,gbff
```

*Converting gbff to gff3*

workdir: `/work/gif4/satheesh/2025/12_PEdwards_Methylation/01_Data/ncbi_genome/ncbi_dataset/data/GCA_040207285.2`

```bash
module load per-bioperl
time bp_genbank2gff3 genomic.gbff
```

*Output:*

```bash
real    2m47.305s
user    2m39.465s
sys     0m5.156s
```

```bash
cp ../../../../../00_Scripts/01_run_gffutils.py . // changed the input filename in the script to `genomic.gbff.gff`

# run the script
conda activate methylation
python 01_run_gffutils.py
```

workdir: `/work/gif4/satheesh/2025/12_PEdwards_Methylation/01_Data/ncbi_genome/ncbi_dataset/data/GCA_040207285.2`

Getting the chromosome list: 

```bash
grep ">" GCA_040207285.2_mMyoRut1_p1.0_genomic.fna | cut -f1 -d ' ' | sort > chm_list.txt
cat chm_list.txt
```

<details>
<summary>chm_list.txt</summary>

```bash
>CM105764.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 1, whole genome shotgun sequence
>CM105765.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 2, whole genome shotgun sequence
>CM105766.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 3, whole genome shotgun sequence
>CM105767.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 4, whole genome shotgun sequence
>CM105768.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 5, whole genome shotgun sequence
>CM105769.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 6, whole genome shotgun sequence
>CM105770.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 7, whole genome shotgun sequence
>CM105771.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 8, whole genome shotgun sequence
>CM105772.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 9, whole genome shotgun sequence
>CM105773.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 10, whole genome shotgun sequence
>CM105774.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 11, whole genome shotgun sequence
>CM105775.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 12, whole genome shotgun sequence
>CM105776.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 13, whole genome shotgun sequence
>CM105777.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 14, whole genome shotgun sequence
>CM105778.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 15, whole genome shotgun sequence
>CM105779.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 16, whole genome shotgun sequence
>CM105780.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 17, whole genome shotgun sequence
>CM105781.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 18, whole genome shotgun sequence
>CM105782.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 19, whole genome shotgun sequence
>CM105783.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 20, whole genome shotgun sequence
>CM105784.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 21, whole genome shotgun sequence
>CM105785.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 22, whole genome shotgun sequence
>CM105786.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 23, whole genome shotgun sequence
>CM105787.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 24, whole genome shotgun sequence
>CM105788.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 25, whole genome shotgun sequence
>CM105789.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 26, whole genome shotgun sequence
>CM105790.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 27, whole genome shotgun sequence
>CM105791.1 Clethrionomys rutilus isolate NRedBackedVole01 chromosome 28, whole genome shotgun sequence
>JBCDMM020000029.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_29, whole genome shotgun sequence
>JBCDMM020000030.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_30, whole genome shotgun sequence
>JBCDMM020000031.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_31, whole genome shotgun sequence
>JBCDMM020000032.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_32, whole genome shotgun sequence
>JBCDMM020000033.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_33, whole genome shotgun sequence
>JBCDMM020000034.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_34, whole genome shotgun sequence
>JBCDMM020000035.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_35, whole genome shotgun sequence
>JBCDMM020000036.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_36, whole genome shotgun sequence
>JBCDMM020000037.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_37, whole genome shotgun sequence
>JBCDMM020000038.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_38, whole genome shotgun sequence
>JBCDMM020000039.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_39, whole genome shotgun sequence
>JBCDMM020000040.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_40, whole genome shotgun sequence
>JBCDMM020000041.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_41, whole genome shotgun sequence
>JBCDMM020000042.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_42, whole genome shotgun sequence
>JBCDMM020000043.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_43, whole genome shotgun sequence
>JBCDMM020000044.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_80, whole genome shotgun sequence
>JBCDMM020000045.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_44, whole genome shotgun sequence
>JBCDMM020000046.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_45, whole genome shotgun sequence
>JBCDMM020000047.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_46, whole genome shotgun sequence
>JBCDMM020000048.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_86, whole genome shotgun sequence
>JBCDMM020000049.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_47, whole genome shotgun sequence
>JBCDMM020000050.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_48, whole genome shotgun sequence
>JBCDMM020000051.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_90, whole genome shotgun sequence
>JBCDMM020000052.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_87, whole genome shotgun sequence
>JBCDMM020000053.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_88, whole genome shotgun sequence
>JBCDMM020000054.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_81, whole genome shotgun sequence
>JBCDMM020000055.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_92, whole genome shotgun sequence
>JBCDMM020000056.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_94, whole genome shotgun sequence
>JBCDMM020000057.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_95, whole genome shotgun sequence
>JBCDMM020000058.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_49, whole genome shotgun sequence
>JBCDMM020000059.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_151, whole genome shotgun sequence
>JBCDMM020000060.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_97, whole genome shotgun sequence
>JBCDMM020000061.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_98, whole genome shotgun sequence
>JBCDMM020000062.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_154, whole genome shotgun sequence
>JBCDMM020000063.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_156, whole genome shotgun sequence
>JBCDMM020000064.1 Clethrionomys rutilus isolate NRedBackedVole01 SUPER_50, whole genome shotgun sequence
>JBCDMM020000065.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_158, whole genome shotgun sequence
>JBCDMM020000066.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_99, whole genome shotgun sequence
>JBCDMM020000067.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_160, whole genome shotgun sequence
>JBCDMM020000068.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_161, whole genome shotgun sequence
>JBCDMM020000069.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_163, whole genome shotgun sequence
>JBCDMM020000070.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_164, whole genome shotgun sequence
>JBCDMM020000071.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_166, whole genome shotgun sequence
>JBCDMM020000072.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_100, whole genome shotgun sequence
>JBCDMM020000073.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_168, whole genome shotgun sequence
>JBCDMM020000074.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_169, whole genome shotgun sequence
>JBCDMM020000075.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_170, whole genome shotgun sequence
>JBCDMM020000076.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_101, whole genome shotgun sequence
>JBCDMM020000077.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_172, whole genome shotgun sequence
>JBCDMM020000078.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_173, whole genome shotgun sequence
>JBCDMM020000079.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_102, whole genome shotgun sequence
>JBCDMM020000080.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_175, whole genome shotgun sequence
>JBCDMM020000081.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_176, whole genome shotgun sequence
>JBCDMM020000082.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_177, whole genome shotgun sequence
>JBCDMM020000083.1 Clethrionomys rutilus isolate NRedBackedVole01 Scaffold_103, whole genome shotgun sequence
>CM079637.2 Clethrionomys rutilus isolate NRedBackedVole01 mitochondrion, complete sequence, whole genome shotgun sequence
```
</details>

**Getting the chromosomes list from the genome assembly from Ensembl database**

workdir: `/work/gif4/satheesh/2025/12_PEdwards_Methylation/01_Data`

```bash
zcat crutilus_genome_softmasked.fa.gz |grep ">" > ensembl_genome_chm_list.txt
```

```bash
cat ensembl_genome_chm_list.txt
```

<details>
<summary>ensembl_genome_chm_list.txt</summary>

```bash
>1 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:1:1:145147427:1
>2 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:2:1:140625002:1
>3 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:3:1:130821871:1
>4 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:4:1:122679291:1
>5 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:5:1:122277089:1
>6 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:6:1:121003038:1
>7 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:7:1:109483915:1
>8 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:8:1:108756654:1
>9 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:9:1:103215971:1
>10 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:10:1:93108836:1
>11 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:11:1:87441830:1
>12 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:12:1:83336943:1
>13 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:13:1:81387531:1
>14 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:14:1:80906502:1
>15 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:15:1:77386962:1
>16 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:16:1:77136395:1
>17 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:17:1:76608580:1
>18 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:18:1:68917561:1
>19 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:19:1:68389746:1
>20 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:20:1:67601279:1
>21 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:21:1:66203084:1
>22 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:22:1:61476354:1
>23 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:23:1:61000435:1
>24 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:24:1:59988077:1
>25 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:25:1:50293231:1
>26 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:26:1:46425774:1
>27 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:27:1:43130027:1
>28 softmasked:chromosome primary_assembly:mMyoRut1_p1.0:28:1:5370437:1
>MT softmasked:chromosome primary_assembly:mMyoRut1_p1.0:MT:1:16350:1
>JBCDMM020000029.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000029.1:1:3667573:1
>JBCDMM020000030.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000030.1:1:3255780:1
>JBCDMM020000031.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000031.1:1:3242292:1
>JBCDMM020000032.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000032.1:1:2865281:1
>JBCDMM020000033.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000033.1:1:2788517:1
>JBCDMM020000034.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000034.1:1:2639075:1
>JBCDMM020000035.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000035.1:1:1965298:1
>JBCDMM020000036.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000036.1:1:1762251:1
>JBCDMM020000037.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000037.1:1:1761906:1
>JBCDMM020000038.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000038.1:1:1688118:1
>JBCDMM020000039.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000039.1:1:1673564:1
>JBCDMM020000040.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000040.1:1:1558053:1
>JBCDMM020000041.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000041.1:1:1554092:1
>JBCDMM020000042.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000042.1:1:1433687:1
>JBCDMM020000043.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000043.1:1:763005:1
>JBCDMM020000044.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000044.1:1:683041:1
>JBCDMM020000045.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000045.1:1:548625:1
>JBCDMM020000046.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000046.1:1:544334:1
>JBCDMM020000047.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000047.1:1:453655:1
>JBCDMM020000048.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000048.1:1:452412:1
>JBCDMM020000049.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000049.1:1:398832:1
>JBCDMM020000050.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000050.1:1:368720:1
>JBCDMM020000051.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000051.1:1:365695:1
>JBCDMM020000052.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000052.1:1:322718:1
>JBCDMM020000053.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000053.1:1:237209:1
>JBCDMM020000054.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000054.1:1:226206:1
>JBCDMM020000055.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000055.1:1:215521:1
>JBCDMM020000056.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000056.1:1:136054:1
>JBCDMM020000057.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000057.1:1:58489:1
>JBCDMM020000058.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000058.1:1:57681:1
>JBCDMM020000059.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000059.1:1:53663:1
>JBCDMM020000060.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000060.1:1:50667:1
>JBCDMM020000061.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000061.1:1:38136:1
>JBCDMM020000062.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000062.1:1:37372:1
>JBCDMM020000063.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000063.1:1:34573:1
>JBCDMM020000064.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000064.1:1:34507:1
>JBCDMM020000065.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000065.1:1:33904:1
>JBCDMM020000066.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000066.1:1:32977:1
>JBCDMM020000067.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000067.1:1:30911:1
>JBCDMM020000068.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000068.1:1:29704:1
>JBCDMM020000069.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000069.1:1:28918:1
>JBCDMM020000070.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000070.1:1:28852:1
>JBCDMM020000071.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000071.1:1:27017:1
>JBCDMM020000072.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000072.1:1:24652:1
>JBCDMM020000073.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000073.1:1:24133:1
>JBCDMM020000074.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000074.1:1:23955:1
>JBCDMM020000075.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000075.1:1:23372:1
>JBCDMM020000076.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000076.1:1:22452:1
>JBCDMM020000077.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000077.1:1:21500:1
>JBCDMM020000078.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000078.1:1:20404:1
>JBCDMM020000079.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000079.1:1:19709:1
>JBCDMM020000080.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000080.1:1:18420:1
>JBCDMM020000081.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000081.1:1:16996:1
>JBCDMM020000082.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000082.1:1:16497:1
>JBCDMM020000083.1 softmasked:primary_assembly mMyoRut1_p1.0:JBCDMM020000083.1:1:14107:1
```
</details>

Both genomes are the same biological assembly *(mMyoRut1_p1.0 / GCA_040207285.2)* but have different sequence IDs. 
Thus, coordinates are identical, but chromosome names differ, blocking overlap. 

## Initial results (not complete, analysis steps to be added)

Test DMRs Filename: `DMRs_LowDecline.bed`

// TODO: Add figure 1 and 2

Figure 1. Genomic feature distribution of Low–Decline DMRs.

Pie chart summarizing the locations of differentially methylated regions (DMRs) across major gene-associated and intergenic features. The majority of DMRs fall within distal intergenic regions (55.6%), followed by introns (1st intron: 11.0%; other introns: 24.7%). Only a small proportion map to promoter-proximal regions (≤1 kb upstream: 3.66%; 1–2 kb upstream: 2.75%) or exonic/UTR regions. This distribution indicates that methylation differences between Low and Decline phases predominantly affect distal regulatory elements and intragenic regulatory domains.

Figure 2. Stacked-bar representation of DMR feature distribution.

Horizontal bar plot illustrating the proportional contribution of each genomic feature category to the full Low–Decline DMR set. Colors correspond to the same categories as in Figure 1. This representation highlights the strong enrichment of DMRs within intergenic and intronic compartments and the relative scarcity of promoter-region methylation changes.

// Script: `annotate_DMRs_strand_aware.r`

```R
library(rtracklayer)
library(dplyr)

# Load nearest-gene table
df <- read.csv("DMR_nearest_genes_LowDecline.csv")

# Load GFF3 to get gene strand
gff <- import("crutilus_genes_NCBI.gff3")
genes <- gff[gff$type == "gene"]

gene_meta <- data.frame(
  gene_id = genes$ID,
  gene_strand = as.character(strand(genes)),
  stringsAsFactors = FALSE
)

# Merge strand info
df2 <- df %>% left_join(gene_meta, by = "gene_id")

# Compute TSS/TES
df2 <- df2 %>% mutate(
  gene_strand = ifelse(is.na(gene_strand), "+", gene_strand),
  TSS = ifelse(gene_strand == "+", gene_start, gene_end),
  TES = ifelse(gene_strand == "+", gene_end, gene_start)
)

# Overlap check
df2 <- df2 %>% mutate(
  overlaps_gene = (DMR_end >= gene_start & DMR_start <= gene_end)
)

# Distance to TSS/TES (absolute)
df2 <- df2 %>% mutate(
  distance_bp = pmin(abs(DMR_start - TSS), abs(DMR_end - TSS)),
  dist_to_TSS = pmin(abs(DMR_start - TSS), abs(DMR_end - TSS)),
  dist_to_TES = pmin(abs(DMR_start - TES), abs(DMR_end - TES))
)

# Biological direction (strand-aware)
df2 <- df2 %>% mutate(
  direction = case_when(
    # Plus strand
    DMR_end < gene_start & gene_strand == "+" ~ "upstream",
    DMR_start > gene_end & gene_strand == "+" ~ "downstream",
    
    # Minus strand
    DMR_end < gene_start & gene_strand == "-" ~ "downstream",
    DMR_start > gene_end & gene_strand == "-" ~ "upstream",
    
    TRUE ~ "overlapping"
  )
)

# Feature assignment (≥2 kb threshold)
df2 <- df2 %>% mutate(
  feature_strand = case_when(
    overlaps_gene ~ "gene_body",
    !overlaps_gene & dist_to_TSS <= 2000 ~ "promoter",
    !overlaps_gene & dist_to_TES <= 2000 ~ "downstream",
    TRUE ~ "intergenic"
  )
)

# Use GFF exon coordinates to annotate exon/intron
exons <- gff[gff$type == "exon"]
exon_df <- data.frame(
  exon_chr = as.character(seqnames(exons)),
  exon_start = start(exons),
  exon_end = end(exons),
  gene_id = exons$Parent,
  stringsAsFactors = FALSE
)

# Convert Parent field to match gene_id naming
exon_df$gene_id <- sub("transcript:|gene:", "", exon_df$gene_id.value)

# Join exon info
library(GenomicRanges)

dmr_gr <- GRanges(df2$DMR_chr, IRanges(df2$DMR_start, df2$DMR_end))
exon_gr <- GRanges(exon_df$exon_chr, IRanges(exon_df$exon_start, exon_df$exon_end), gene_id = exon_df$gene_id)

hits_exon <- findOverlaps(dmr_gr, exon_gr)

# Mark exonic DMRs
df2$exon_overlap <- FALSE
df2$exon_overlap[queryHits(hits_exon)] <- TRUE

# Define intron if gene_body but not exon
df2 <- df2 %>% mutate(
  refined_feature = case_when(
    exon_overlap ~ "exon",
    feature_strand == "gene_body" & !exon_overlap ~ "intron",
    TRUE ~ feature_strand
  )
)

# Summary plot for refined_feature
library(ggplot2)
png("DMR_feature_summary.png", width=800, height=600)
ggplot(df2, aes(x=refined_feature)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  ggtitle("DMR Feature Distribution (strand-aware, exon/intron refined)") +
  xlab("Feature Type") + ylab("Count")
dev.off()

# Save output
write.csv(df2, "DMR_nearest_genes_LowDecline_strandAware.csv", row.names = FALSE)
```

## Final complete analysis 

Running the Rscript `00_Scripts/02_batch_strand_aware_annotation.R`. This will analyse all DML and DMR results files obtained from the DSS analysis.

```bash 
Rscript 00_Scripts/02_batch_strand_aware_annotation.R
```

```
=== Summary Statistics ===


 DML_nearest_genes_LowDecline_strandAware.csv
  Feature distribution:

downstream       exon intergenic     intron   promoter 
        12         18        554        265          3 
  Direction distribution:

 downstream overlapping    upstream 
        256         281         315 

 DML_nearest_genes_PeakDecline_strandAware.csv
  Feature distribution:

downstream       exon intergenic     intron   promoter 
         8          7        479        243          9 
  Direction distribution:

 downstream overlapping    upstream 
        205         251         290 

 DML_nearest_genes_PeakLow_strandAware.csv
  Feature distribution:

downstream       exon intergenic     intron   promoter 
        13         19        841        479         15 
  Direction distribution:

 downstream overlapping    upstream 
        387         500         480 

 DMR_nearest_genes_LowDecline_strandAware.csv
  Feature distribution:

downstream       exon intergenic     intron   promoter 
         4         14        243        170          7 
  Direction distribution:

 downstream overlapping    upstream 
        107         186         145 

 DMR_nearest_genes_PeakDecline_strandAware.csv
  Feature distribution:

downstream       exon intergenic     intron   promoter 
         6          9        207        132          3 
  Direction distribution:

 downstream overlapping    upstream 
         89         141         127 

 DMR_nearest_genes_PeakLow_strandAware.csv
  Feature distribution:

downstream       exon intergenic     intron   promoter 
         4         14        279        166          3 
  Direction distribution:

 downstream overlapping    upstream 
        134         182         150 
```

### Understanding the Output

The script produces two types of classifications for each DML/DMR:

**Feature Distribution** — Classifies where the DML/DMR falls relative to gene structure:

| Feature | Definition |
|---------|------------|
| **promoter** | Within 2kb upstream of the transcription start site (TSS) |
| **exon** | Overlaps an exon |
| **intron** | Within gene body but not overlapping an exon |
| **downstream** | Within 2kb downstream of the transcription end site (TES) |
| **intergenic** | More than 2kb from any gene |

**Direction Distribution** — Classifies the DML/DMR position relative to the nearest gene (strand-aware):

| Direction | Definition |
|-----------|------------|
| **upstream** | Located before the gene (5' side, accounting for strand) |
| **downstream** | Located after the gene (3' side, accounting for strand) |
| **overlapping** | Overlaps the gene body |

**Note:** "Upstream" in Direction and "promoter" in Feature are related but not identical. A DML/DMR classified as "upstream" in direction will be labeled "promoter" in feature only if it is ≤2kb from the TSS; otherwise it is classified as "intergenic".
