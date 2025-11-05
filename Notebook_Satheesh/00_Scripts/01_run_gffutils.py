#!/usr/bin/env python

import gffutils

# Create or load a database
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