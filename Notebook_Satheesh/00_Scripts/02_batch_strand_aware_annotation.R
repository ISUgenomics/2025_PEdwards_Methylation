#!/usr/bin/env Rscript
# =============================================================================
# Batch Strand-Aware DMR/DML Annotation Script
# =============================================================================
# Description: Processes all DML and DMR files, finds nearest genes, and 
#              performs strand-aware annotation with exon/intron refinement.
# 
# Input:  CSV files in dmls_dmrs_dss/ folder (DMLs_*.csv, DMRs_*.csv)
# Output: Annotated CSV files and summary plots in results/ folder
# =============================================================================

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(txdbmaker)
library(dplyr)
library(ggplot2)

# =============================================================================
# Configuration - Update these paths as needed
# =============================================================================
base_dir    <- "/Users/vsatheesh/Documents/ISU/GIF/gif_projects/2025/04_PEdwards_Methylation"
input_dir   <- file.path(base_dir, "dmls_dmrs_dss")
output_dir  <- file.path(base_dir, "results")
gff_file    <- file.path(base_dir, "crutilus_genes_NCBI.gff3")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load GFF3 and prepare gene/exon annotations (done once)
# =============================================================================
message("Loading GFF3 annotation...")
gff <- import(gff_file)

# Extract genes
genes <- gff[gff$type == "gene"]
gene_meta <- data.frame(
  gene_id = genes$ID,
  gene_strand = as.character(strand(genes)),
  stringsAsFactors = FALSE
)

# Extract exons for exon/intron annotation
exons <- gff[gff$type == "exon"]
exon_df <- data.frame(
  exon_chr = as.character(seqnames(exons)),
  exon_start = start(exons),
  exon_end = end(exons),
  gene_id = sub("transcript:|gene:", "", exons$Parent),
  stringsAsFactors = FALSE
)
exon_gr <- GRanges(exon_df$exon_chr, 
                   IRanges(exon_df$exon_start, exon_df$exon_end), 
                   gene_id = exon_df$gene_id)

# Build TxDb for ChIPseeker (optional, for additional annotation)
message("Building TxDb...")
txdb <- makeTxDbFromGFF(gff_file, format = "gff3")

# =============================================================================
# Function: Find nearest genes and create initial annotation
# =============================================================================
find_nearest_genes <- function(input_file, file_type) {
  message(paste("  Finding nearest genes for:", basename(input_file)))
  
  # Read input file
  data <- read.csv(input_file)
  

  # Determine column names based on file type
  if (file_type == "DML") {
    chr_col <- "chr"
    pos_col <- "pos"
    # For DMLs, start and end are the same position
    data_gr <- GRanges(
      seqnames = data[[chr_col]],
      ranges = IRanges(start = data[[pos_col]], end = data[[pos_col]])
    )
  } else {
    # DMR files have start and end columns
    chr_col <- "chr"
    start_col <- "start"
    end_col <- "end"
    data_gr <- GRanges(
      seqnames = data[[chr_col]],
      ranges = IRanges(start = data[[start_col]], end = data[[end_col]])
    )
  }
  
  # Find nearest gene
  # Some DMLs/DMRs may be on chromosomes with no genes - handle missing hits
  hit <- distanceToNearest(data_gr, genes)
  
  # Build output dataframe
  prefix <- ifelse(file_type == "DML", "DML", "DMR")
  
  # Initialize df with original data
  df <- data.frame(data, stringsAsFactors = FALSE)
  
  # Add standardized columns for downstream processing
  if (file_type == "DML") {
    df[[paste0(prefix, "_chr")]] <- data[[chr_col]]
    df[[paste0(prefix, "_start")]] <- data[[pos_col]]
    df[[paste0(prefix, "_end")]] <- data[[pos_col]]
  } else {
    df[[paste0(prefix, "_chr")]] <- data[[chr_col]]
    df[[paste0(prefix, "_start")]] <- data[[start_col]]
    df[[paste0(prefix, "_end")]] <- data[[end_col]]
  }
  
  # Initialize gene columns with NA
  df$distance <- NA_integer_
  df$gene_chr <- NA_character_
  df$gene_start <- NA_integer_
  df$gene_end <- NA_integer_
  df$gene_id <- NA_character_
  df$gene_name <- NA_character_
  df$gene_desc <- NA_character_
  
  # Fill in values only for rows with hits
  matched_idx <- queryHits(hit)
  nearestGenes <- genes[subjectHits(hit)]
  
  df$distance[matched_idx] <- mcols(hit)$distance
  df$gene_chr[matched_idx] <- as.character(seqnames(nearestGenes))
  df$gene_start[matched_idx] <- start(nearestGenes)
  df$gene_end[matched_idx] <- end(nearestGenes)
  df$gene_id[matched_idx] <- mcols(nearestGenes)$ID
  df$gene_name[matched_idx] <- mcols(nearestGenes)$Name
  df$gene_desc[matched_idx] <- mcols(nearestGenes)$description
  
  # Report unmatched rows
  n_unmatched <- nrow(df) - length(matched_idx)
  if (n_unmatched > 0) {
    message(paste("  Warning:", n_unmatched, "rows had no nearest gene (likely on scaffolds without genes)"))
  }
  
  return(list(df = df, prefix = prefix))
}

# =============================================================================
# Function: Perform strand-aware annotation
# =============================================================================
strand_aware_annotation <- function(df, prefix) {
  message("  Performing strand-aware annotation...")
  
  # Define column names based on prefix
  chr_col <- paste0(prefix, "_chr")
  start_col <- paste0(prefix, "_start")
  end_col <- paste0(prefix, "_end")
  
  # Merge strand info
  df2 <- df %>% left_join(gene_meta, by = "gene_id")
  
  # Compute TSS/TES based on strand
  df2 <- df2 %>% mutate(
    gene_strand = ifelse(is.na(gene_strand), "+", gene_strand),
    TSS = ifelse(gene_strand == "+", gene_start, gene_end),
    TES = ifelse(gene_strand == "+", gene_end, gene_start)
  )
  
  # Check for overlap with gene body
  df2 <- df2 %>% mutate(
    overlaps_gene = (.data[[end_col]] >= gene_start & .data[[start_col]] <= gene_end)
  )
  
  # Distance to TSS/TES (absolute)
  df2 <- df2 %>% mutate(
    distance_bp = pmin(abs(.data[[start_col]] - TSS), abs(.data[[end_col]] - TSS)),
    dist_to_TSS = pmin(abs(.data[[start_col]] - TSS), abs(.data[[end_col]] - TSS)),
    dist_to_TES = pmin(abs(.data[[start_col]] - TES), abs(.data[[end_col]] - TES))
  )
  
  # Biological direction (strand-aware)
  df2 <- df2 %>% mutate(
    direction = case_when(
      # Plus strand
      .data[[end_col]] < gene_start & gene_strand == "+" ~ "upstream",
      .data[[start_col]] > gene_end & gene_strand == "+" ~ "downstream",
      # Minus strand
      .data[[end_col]] < gene_start & gene_strand == "-" ~ "downstream",
      .data[[start_col]] > gene_end & gene_strand == "-" ~ "upstream",
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
  
  # Exon/intron refinement using GenomicRanges overlap
  data_gr <- GRanges(df2[[chr_col]], IRanges(df2[[start_col]], df2[[end_col]]))
  hits_exon <- findOverlaps(data_gr, exon_gr)
  
  df2$exon_overlap <- FALSE
  df2$exon_overlap[queryHits(hits_exon)] <- TRUE
  
  # Define refined feature
  df2 <- df2 %>% mutate(
    refined_feature = case_when(
      exon_overlap ~ "exon",
      feature_strand == "gene_body" & !exon_overlap ~ "intron",
      TRUE ~ feature_strand
    )
  )
  
  return(df2)
}

# =============================================================================
# Function: Generate summary plot
# =============================================================================
generate_plot <- function(df, output_file, title) {
  message("  Generating summary plot...")
  
  png(output_file, width = 800, height = 600)
  p <- ggplot(df, aes(x = refined_feature)) +
    geom_bar(fill = "steelblue") +
    theme_minimal() +
    ggtitle(title) +
    xlab("Feature Type") + 
    ylab("Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()
}

# =============================================================================
# Main Processing Loop
# =============================================================================
message("\n=== Starting Batch Processing ===\n")

# Get all CSV files in input directory
input_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

for (input_file in input_files) {
  filename <- basename(input_file)
  message(paste("\nProcessing:", filename))
  
  # Determine file type (DML or DMR)
  if (grepl("^DMLs_", filename)) {
    file_type <- "DML"
  } else if (grepl("^DMRs_", filename)) {
    file_type <- "DMR"
  } else {
    message(paste("  Skipping unknown file type:", filename))
    next
  }
  
  # Extract comparison name (e.g., "LowDecline" from "DMRs_LowDecline.csv")
  comparison <- sub("^(DMLs_|DMRs_)", "", sub("\\.csv$", "", filename))
  
  # Step 1: Find nearest genes
  result <- find_nearest_genes(input_file, file_type)
  df <- result$df
  prefix <- result$prefix
  
  # Step 2: Strand-aware annotation
  df_annotated <- strand_aware_annotation(df, prefix)
  
  # Step 3: Save annotated CSV
  output_csv <- file.path(output_dir, paste0(file_type, "_nearest_genes_", comparison, "_strandAware.csv"))
  write.csv(df_annotated, output_csv, row.names = FALSE)
  message(paste("  Saved:", basename(output_csv)))
  
  # Step 4: Generate summary plot
  plot_file <- file.path(output_dir, paste0(file_type, "_feature_summary_", comparison, ".png"))
  plot_title <- paste0(file_type, " Feature Distribution - ", comparison, " (strand-aware)")
  generate_plot(df_annotated, plot_file, plot_title)
  message(paste("  Saved:", basename(plot_file)))
}

message("\n=== Batch Processing Complete ===")
message(paste("Results saved to:", output_dir))

# =============================================================================
# Summary Statistics
# =============================================================================
message("\n=== Summary Statistics ===\n")

# Read all output files and summarize
output_files <- list.files(output_dir, pattern = "_strandAware\\.csv$", full.names = TRUE)

for (out_file in output_files) {
  df <- read.csv(out_file)
  message(paste("\n", basename(out_file)))
  message("  Feature distribution:")
  print(table(df$refined_feature))
  message("  Direction distribution:")
  print(table(df$direction))
}
