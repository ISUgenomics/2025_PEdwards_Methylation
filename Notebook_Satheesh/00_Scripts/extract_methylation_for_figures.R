#!/usr/bin/env Rscript
# Extract per-sample methylation levels at key loci for Figures 2 and 3.
# Input : BSobj_2026.RData (BSseq object)
# Output: methylation_for_figures.csv (long format)

library(bsseq)
library(GenomicRanges)
library(dplyr)

# Load BSseq object
message("Loading BSseq object...")
load("BSobj_2026.RData")   # loads object: BSobj

stopifnot(exists("BSobj"), is(BSobj, "BSseq"))
message("BSseq object loaded: ", nrow(BSobj), " CpG sites x ",
        ncol(BSobj), " samples")
message("Sample names: ", paste(sampleNames(BSobj), collapse = " | "))

# Sample metadata
# Order must match the sample order used in makeBSseqData()
sample_meta <- data.frame(
  sample_label = c(
    "P (2017)",  "T (2018)",  "AP (2019)", "M (2017)",  "AZ (2019)",
    "L (2017)",  "AA (2018)", "AU (2019)", "U (2018)",  "AR (2019)",
    "AM (2019)", "G (2017)",  "AO (2019)", "W (2018)",  "AT (2019)",
    "J (2017)",  "I (2017)",  "AV (2019)", "O (2017)",  "A (2017)",
    "AC (2018)", "S (2018)",  "X (2018)",  "Y (2018)",  "AL (2019)"
  ),
  sample_id = c(
    "P",  "T",  "AP", "M",  "AZ",
    "L",  "AA", "AU", "U",  "AR",
    "AM", "G",  "AO", "W",  "AT",
    "J",  "I",  "AV", "O",  "A",
    "AC", "S",  "X",  "Y",  "AL"
  ),
  group = c(
    "Peak",    "Decline", "Low",     "Peak",    "Low",
    "Peak",    "Decline", "Low",     "Decline", "Low",
    "Low",     "Peak",    "Low",     "Decline", "Low",
    "Peak",    "Peak",    "Low",     "Peak",    "Peak",
    "Decline", "Decline", "Decline", "Decline", "Low"
  ),
  stringsAsFactors = FALSE
)

if (!all(sampleNames(BSobj) == sample_meta$sample_label)) {
  warning(
    "Sample name mismatch detected!\n",
    "BSobj names : ", paste(sampleNames(BSobj), collapse=", "), "\n",
    "Expected    : ", paste(sample_meta$sample_label, collapse=", ")
  )
}

# Individual CpG sites (Fig 2A)
# Atp5g3 promoter: DMLs 909 and 910 bp upstream of TSS. Both hypomethylated in
# Peak vs Decline; the 910 bp site also in Low vs Decline.
dml_loci <- GRanges(
  seqnames = c("CM105785.1",        "CM105785.1"),
  ranges   = IRanges(
    start  = c(51064559,             51064560),
    end    = c(51064559,             51064560)
  ),
  locus    = c("Atp5g3_909bp_TSS",  "Atp5g3_910bp_TSS"),
  gene     = c("Atp5g3",            "Atp5g3"),
  figure   = c("Fig2A",             "Fig2A"),
  comparison_sig = c("Peak-Decline", "Peak-Decline + Low-Decline")
)

# DMR regions (Figs 2B, 2C, 3A, 3B)
dmr_loci <- GRanges(
  seqnames = c(
    "CM105767.1",
    "CM105766.1",
    "CM105779.1",
    "CM105765.1",
    "CM105765.1",
    "CM105787.1",
    "CM105787.1"
  ),
  ranges = IRanges(
    start = c(16281578, 95637506, 11606140, 109904477, 110053310, 31299066, 31346726),
    end   = c(16281629, 95637810, 11606212, 109904529, 110053375, 31299137, 31346820)
  ),
  locus  = c(
    "Cyp1a1_promoter_52bp",
    "Cyp2j13_promoter_305bp",
    "Mif_promoter_73bp",
    "Sgk1_exon_53bp",
    "Sgk1_intergenic_66bp",
    "Igf1r_intron_72bp",
    "Igf1r_intron_95bp"
  ),
  gene   = c("Cyp1a1", "Cyp2j13", "Mif", "Sgk1", "Sgk1", "Igf1r", "Igf1r"),
  figure = c("Fig2B",  "Fig2B",   "Fig2C","Fig3A","Fig3A","Fig3B", "Fig3B"),
  comparison_sig = c(
    "Low-Decline",
    "Peak-Decline",
    "Low-Decline",
    "Low-Decline",
    "Peak-Low",
    "Peak-Decline",
    "Low-Decline + Peak-Low"
  )
)

extract_single_cpg <- function(bsobj, query_gr) {

  results <- vector("list", length(query_gr))

  for (i in seq_along(query_gr)) {

    hits <- findOverlaps(query_gr[i], bsobj)

    if (length(hits) == 0) {
      warning("No CpG in BSobj matches position: ", start(query_gr[i]),
              " on ", seqnames(query_gr[i]))
      next
    }

    cpg_idx <- subjectHits(hits)[1]  # use first match if multiple

    Cov <- as.numeric(getCoverage(bsobj)[cpg_idx, ])
    M   <- as.numeric(getCoverage(bsobj, type = "M")[cpg_idx, ])

    meth_pct <- ifelse(Cov > 0, 100 * M / Cov, NA_real_)

    results[[i]] <- data.frame(
      figure         = query_gr$figure[i],
      gene           = query_gr$gene[i],
      locus          = query_gr$locus[i],
      comparison_sig = query_gr$comparison_sig[i],
      chr            = as.character(seqnames(query_gr[i])),
      region_start   = start(query_gr[i]),
      region_end     = end(query_gr[i]),
      n_cpgs         = 1L,
      sample_label   = sampleNames(bsobj),
      methylation_pct= meth_pct,
      coverage       = Cov,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(results)
}

extract_region_mean <- function(bsobj, query_gr) {

  results <- vector("list", length(query_gr))

  for (i in seq_along(query_gr)) {

    bsobj_sub <- subsetByOverlaps(bsobj, query_gr[i])
    n_cpgs    <- nrow(bsobj_sub)

    if (n_cpgs == 0) {
      warning("No CpGs found in region: ", query_gr$locus[i])
      next
    }
    message("  ", query_gr$locus[i], " — ", n_cpgs, " CpGs")

    # Weighted mean per sample: sum(M) / sum(Cov) avoids small-N bias
    Cov_mat <- getCoverage(bsobj_sub)          # CpGs × samples
    M_mat   <- getCoverage(bsobj_sub, type="M")

    Cov_sum <- colSums(Cov_mat, na.rm = TRUE)
    M_sum   <- colSums(M_mat,   na.rm = TRUE)

    meth_pct <- ifelse(Cov_sum > 0, 100 * M_sum / Cov_sum, NA_real_)

    results[[i]] <- data.frame(
      figure         = query_gr$figure[i],
      gene           = query_gr$gene[i],
      locus          = query_gr$locus[i],
      comparison_sig = query_gr$comparison_sig[i],
      chr            = as.character(seqnames(query_gr[i])),
      region_start   = start(query_gr[i]),
      region_end     = end(query_gr[i]),
      n_cpgs         = n_cpgs,
      sample_label   = sampleNames(bsobj),
      methylation_pct= meth_pct,
      coverage       = Cov_sum,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(results)
}

message("\n--- Extracting individual CpG sites (Fig 2A) ---")
dml_data <- extract_single_cpg(BSobj, dml_loci)

message("\n--- Extracting DMR region means (Figs 2B, 2C, 3A, 3B) ---")
dmr_data <- extract_region_mean(BSobj, dmr_loci)

# Combine, annotate groups, and save
all_data <- dplyr::bind_rows(dml_data, dmr_data) %>%
  dplyr::left_join(sample_meta, by = "sample_label") %>%
  dplyr::mutate(
    group = factor(group, levels = c("Peak", "Decline", "Low"))
  ) %>%
  dplyr::select(
    figure, gene, locus, comparison_sig,
    chr, region_start, region_end, n_cpgs,
    sample_label, sample_id, group,
    methylation_pct, coverage
  ) %>%
  dplyr::arrange(figure, locus, group, sample_id)

out_file <- "methylation_for_figures.csv"
write.csv(all_data, out_file, row.names = FALSE)
message("\nSaved: ", out_file, "  (", nrow(all_data), " rows)")

message("\n=== Summary of extracted loci ===")
summary_tbl <- all_data %>%
  dplyr::group_by(figure, gene, locus, comparison_sig) %>%
  dplyr::summarise(
    n_samples      = sum(!is.na(methylation_pct)),
    n_na           = sum(is.na(methylation_pct)),
    mean_meth_pct  = round(mean(methylation_pct, na.rm = TRUE), 1),
    n_cpgs_in_region = unique(n_cpgs),
    .groups = "drop"
  )
print(as.data.frame(summary_tbl), row.names = FALSE)

message("\nGroup means per locus:")
group_means <- all_data %>%
  dplyr::group_by(figure, gene, locus, group) %>%
  dplyr::summarise(
    mean_meth = round(mean(methylation_pct, na.rm = TRUE), 1),
    n         = sum(!is.na(methylation_pct)),
    .groups   = "drop"
  )
print(as.data.frame(group_means), row.names = FALSE)
