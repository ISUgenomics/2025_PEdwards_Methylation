#!/usr/bin/env Rscript
# Figures 2 and 3, with a common y-axis range across all panels.
# Input: methylation_for_figures.csv

library(ggplot2)
library(ggsignif)
library(patchwork)
library(dplyr)

df <- read.csv("methylation_for_figures.csv", stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("Peak", "Decline", "Low"))

cat("Loaded", nrow(df), "rows\n")
cat("Loci:", paste(unique(df$locus), collapse = ", "), "\n\n")

GROUP_COLORS <- c("Peak" = "#2E7D8C", "Decline" = "#D4860A", "Low" = "#5B6E8C")
GROUP_LABELS <- c(
  "Peak"    = "Peak\n(2017)",
  "Decline" = "Decline\n(2018)",
  "Low"     = "Low\n(2019)"
)

# Shared y-axis: same visible range and ticks for all panels
Y_MIN        <- 0
Y_DATA_MAX   <- 100
Y_HEADROOM   <- 125   # space for significance brackets
Y_LIMITS     <- c(Y_MIN, Y_HEADROOM)
Y_BREAKS     <- c(0, 25, 50, 75, 100)

theme_vole <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(
        face = "italic", size = base_size, hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.subtitle      = element_text(
        size = base_size - 1.5, hjust = 0.5,
        colour = "grey40",
        margin = margin(t = 0, b = 5)
      ),
      plot.tag           = element_text(
        size = 12, face = "bold",
        margin = margin(r = 4, b = 4)
      ),
      plot.tag.position  = "topleft",
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(size = base_size - 0.5),
      axis.text.x        = element_text(size = base_size - 1.5, lineheight = 0.85),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      plot.margin        = margin(t = 6, r = 6, b = 6, l = 6)
    )
}

make_panel <- function(locus_name,
                       title,
                       subtitle,
                       sig_pairs  = NULL,
                       panel_tag  = NULL,
                       y_limits   = Y_LIMITS,
                       y_breaks   = Y_BREAKS) {

  sub <- df %>%
    filter(locus == locus_name, !is.na(methylation_pct))

  if (nrow(sub) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No data:\n", locus_name),
                 colour = "red") +
        theme_void() +
        labs(tag = panel_tag)
    )
  }

  p <- ggplot(sub, aes(x = group, y = methylation_pct, fill = group, colour = group)) +
    geom_violin(alpha = 0.28, trim = TRUE, scale = "width", linewidth = 0.25) +
    geom_boxplot(width = 0.13, fill = "white", outlier.shape = NA,
                 linewidth = 0.55, colour = "black") +
    geom_jitter(width = 0.07, size = 1.55, alpha = 0.72,
                shape = 21, fill = "white", stroke = 0.35) +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_colour_manual(values = GROUP_COLORS) +
    scale_x_discrete(labels = GROUP_LABELS) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      y        = "% Methylation",
      tag      = panel_tag
    ) +
    theme_vole()

  if (!is.null(sig_pairs) && length(sig_pairs) > 0) {
    # Fixed bracket heights so all panels look the same
    n_sig <- length(sig_pairs)
    y_positions <- seq(from = 108, by = 8, length.out = n_sig)
    y_positions <- pmin(y_positions, y_limits[2] - 3)

    p <- p +
      geom_signif(
        comparisons      = sig_pairs,
        test             = "wilcox.test",
        test.args        = list(exact = FALSE),
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
        y_position       = y_positions,
        tip_length       = 0.012,
        textsize         = 2.8,
        colour           = "black",
        vjust            = 0.4
      )
  }

  p
}

# Figure 2: promoter DMCs and DMRs
message("Building Figure 2...")

p2A_1 <- make_panel(
  locus_name = "Atp5g3_909bp_TSS",
  title      = "Atp5g3",
  subtitle   = "909 bp from TSS",
  sig_pairs  = list(c("Peak", "Decline")),
  panel_tag  = "A"
)

p2A_2 <- make_panel(
  locus_name = "Atp5g3_910bp_TSS",
  title      = "Atp5g3",
  subtitle   = "910 bp from TSS",
  sig_pairs  = list(c("Peak", "Decline"), c("Low", "Decline"))
)

p2B_1 <- make_panel(
  locus_name = "Cyp1a1_promoter_52bp",
  title      = "Cyp1a1",
  subtitle   = "Promoter DMR (52 bp)",
  sig_pairs  = list(c("Low", "Decline")),
  panel_tag  = "B"
)

p2B_2 <- make_panel(
  locus_name = "Cyp2j13_promoter_305bp",
  title      = "Cyp2j13",
  subtitle   = "Promoter DMR (305 bp)",
  sig_pairs  = list(c("Peak", "Decline"))
)

p2C <- make_panel(
  locus_name = "Mif_promoter_73bp",
  title      = "Mif",
  subtitle   = "Promoter DMR (73 bp)",
  sig_pairs  = list(c("Low", "Decline")),
  panel_tag  = "C"
)

fig2 <- (p2A_1 | p2A_2 | p2B_1 | p2B_2 | p2C) +
  plot_layout(ncol = 5, widths = rep(1, 5)) +
  plot_annotation(
    title = "Figure 2. Differentially methylated sites in gene promoters",
    theme = theme(
      plot.title = element_text(size = 10, face = "bold",
                                hjust = 0.5, margin = margin(b = 8))
    )
  )

ggsave(
  "Figure2_promoter_DMCs_DMRs_v2.tiff",
  fig2,
  width       = 12.5,
  height      = 4.8,
  dpi         = 300,
  compression = "lzw",
  device      = "tiff"
)
ggsave(
  "Figure2_promoter_DMCs_DMRs_v2.png",
  fig2,
  width       = 12.5,
  height      = 4.8,
  dpi         = 300,
  device      = "png"
)
message("  Saved: Figure2_promoter_DMCs_DMRs_v2.tiff")
message("  Saved: Figure2_promoter_DMCs_DMRs_v2.png")

# Figure 3: gene body DMRs
message("Building Figure 3...")

p3A_1 <- make_panel(
  locus_name = "Sgk1_exon_53bp",
  title      = "Sgk1",
  subtitle   = "Exon DMR (53 bp)",
  sig_pairs  = list(c("Low", "Decline")),
  panel_tag  = "A"
)

p3A_2 <- make_panel(
  locus_name = "Sgk1_intergenic_66bp",
  title      = "Sgk1",
  subtitle   = "Intergenic DMR (66 bp)",
  sig_pairs  = list(c("Peak", "Low"))
)

p3B_1 <- make_panel(
  locus_name = "Igf1r_intron_72bp",
  title      = "Igf1r",
  subtitle   = "Intron DMR (72 bp)",
  sig_pairs  = list(c("Peak", "Decline")),
  panel_tag  = "B"
)

p3B_2 <- make_panel(
  locus_name = "Igf1r_intron_95bp",
  title      = "Igf1r",
  subtitle   = "Intron DMR (95 bp)",
  sig_pairs  = list(c("Low", "Decline"), c("Peak", "Low"))
)

fig3 <- (p3A_1 | p3A_2 | p3B_1 | p3B_2) +
  plot_layout(ncol = 4, widths = rep(1, 4)) +
  plot_annotation(
    title = "Figure 3. Differentially methylated regions in gene bodies",
    theme = theme(
      plot.title = element_text(size = 10, face = "bold",
                                hjust = 0.5, margin = margin(b = 8))
    )
  )

ggsave(
  "Figure3_genebody_DMRs_v2.tiff",
  fig3,
  width       = 10.0,
  height      = 4.8,
  dpi         = 300,
  compression = "lzw",
  device      = "tiff"
)
ggsave(
  "Figure3_genebody_DMRs_v2.png",
  fig3,
  width       = 10.0,
  height      = 4.8,
  dpi         = 300,
  device      = "png"
)
message("  Saved: Figure3_genebody_DMRs_v2.tiff")
message("  Saved: Figure3_genebody_DMRs_v2.png")

message("\n=== Group means per locus ===")
summary_tbl <- df %>%
  filter(!is.na(methylation_pct)) %>%
  group_by(figure, locus, group) %>%
  summarise(
    n    = n(),
    mean = round(mean(methylation_pct), 1),
    sd   = round(sd(methylation_pct), 1),
    .groups = "drop"
  )

print(as.data.frame(summary_tbl), row.names = FALSE)

message("\nDone.")
