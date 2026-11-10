#!/usr/bin/env Rscript
# Figures 1-6 (PNG + PDF). Reads only combined/ tables. The prior CHM13
# PRO-seq analysis nominated the fixed promoter; it is not re-fit here.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(scales)})
args <- commandArgs(trailingOnly = TRUE)
comb <- if (length(args) >= 1) args[1] else "combined"
figd <- file.path(comb, "figures"); dir.create(figd, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.delim(file.path(comb, f), stringsAsFactors = FALSE, check.names = FALSE)
res <- rd("all_array_results.tsv"); prof <- rd("per_position_array_weighted_profile.tsv")
gsum <- rd("group_summary.tsv"); qc <- rd("anchor_qc.tsv"); ne <- rd("non_estimable_arrays.tsv")

pal <- c(telomeric_recombination = "#C0392B", telomeric_no_recombination = "#2E86C1", non_telomeric_no_recombination = "#7D6608", ALL = "grey30")
group_labels <- c(
  telomeric_recombination = "Recombination-positive",
  telomeric_no_recombination = "No recombination",
  non_telomeric_no_recombination = "Non-telomeric, no recombination"
)
res$group <- factor(res$recombination_class, levels = names(pal)[1:3])
res$group <- droplevels(res$group)
P0 <- 78; P1 <- 127
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92")))
save2 <- function(p, name, w, h) {
  ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 200, device = ragg::agg_png)
  ggsave(file.path(figd, paste0(name, ".pdf")), p, width = w, height = h,
         device = grDevices::pdf, useDingbats = FALSE)
}
shade <- annotate("rect", xmin = P0 - 0.5, xmax = P1 + 0.5, ymin = -Inf, ymax = Inf, fill = "gold", alpha = 0.25)

# ---- Fig 1: per-position divergence + PRO-seq annotation + callable-array count
nar <- max(prof$n_arrays_with_alignment)
p1a <- ggplot(prof, aes(alignment_column_1based)) + shade +
  geom_ribbon(aes(ymin = array_q25_divergence, ymax = array_q75_divergence), fill = "grey80", alpha = 0.6) +
  geom_line(aes(y = array_mean_rolling11_divergence, colour = "mean (11-col rolling)"), linewidth = 0.7) +
  geom_line(aes(y = array_median_rolling11_divergence, colour = "median (11-col rolling)"), linewidth = 0.7) +
  geom_line(aes(y = array_mean_divergence, colour = "mean (per column)"), linewidth = 0.25, alpha = 0.6) +
  scale_colour_manual(values = c("mean (11-col rolling)" = "#1B4F72", "median (11-col rolling)" = "#E67E22", "mean (per column)" = "grey40"), name = NULL) +
  labs(y = "Pairwise divergence", x = NULL, title = sprintf("Array-weighted L1PA10 divergence (n = %d arrays)", nar),
       subtitle = "Shaded: fixed promoter, alignment columns 78-127. Grey band: IQR across arrays of per-column divergence.") +
  theme(legend.position = "top")
p1b <- ggplot(prof, aes(alignment_column_1based)) + shade +
  geom_line(aes(y = projected_chr21_NEG_PROseq), colour = "grey55", linewidth = 0.3) +
  geom_line(aes(y = rolling_7col_PROseq), colour = "#8E44AD", linewidth = 0.8) +
  labs(y = "CHM13 chr21\nNEG PRO-seq", x = NULL,
       title = "Projected CHM13 NEG-strand PRO-seq nominated the fixed promoter",
       subtitle = "Fixed before the HPRC analysis; not used to place anchors or re-select/reposition windows")
p1c <- ggplot(prof, aes(alignment_column_1based, n_arrays_callable_ge80pct)) + shade + geom_area(fill = "#5D6D7E", alpha = 0.7) +
  labs(x = "L1PA10 alignment column (300-column reference)", y = "Arrays callable\n(>=80% of copies)")
p1 <- p1a / p1b / p1c + plot_layout(heights = c(3, 1.6, 1.2))
save2(p1, "fig1_per_position_divergence", 10, 9)

# ---- Fig 2 / 4: effect distributions by group (array = unit)
effect_plot <- function(df, col, scol, title, ylab, sumname) {
  d <- df[df[[scol]] == "estimable", ]; d$e <- d[[col]]
  gs <- gsum[gsum$analysis == sumname & gsum$group != "ALL" & gsum$n_estimable > 0, ]
  gs$group <- factor(gs$group, levels = levels(res$group))
  lab <- paste0("n=", gs$n_estimable, "\nmean ", sprintf("%.4f", gs$mean_effect), " [", sprintf("%.4f", gs$boot_ci95_lo), ", ", sprintf("%.4f", gs$boot_ci95_hi), "]\n",
                sprintf("%.0f%%", 100 * gs$frac_positive), " > 0; sign-flip p=", signif(gs$signflip_p_one_sided_mean_gt0, 2))
  ggplot(d, aes(group, e, fill = group, colour = group)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
    geom_violin(alpha = 0.25, colour = NA, scale = "width") + geom_jitter(width = 0.18, size = 0.5, alpha = 0.35) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", colour = "black", alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "gold", colour = "black") +
    geom_text(data = data.frame(group = gs$group, lab = lab), aes(group, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.1, size = 2.8) +
    scale_fill_manual(values = pal, guide = "none") + scale_colour_manual(values = pal, guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
    labs(x = NULL, y = ylab, title = title, subtitle = "Effect > 0: promoter less divergent than matched controls",
         caption = "Each point is one array; diamond = array mean; box = median/IQR")
}
save2(effect_plot(res, "effect", "primary_status", "Primary 300-column L1PA10 promoter effect", "Effect (matched-control median - promoter divergence)", "primary_300col_matched"),
      "fig2_primary_effect_by_group", 7.5, 6)
save2(effect_plot(res, "bb_effect", "bb_status", "Conserved-backbone sensitivity: promoter effect", "Backbone effect (distal-backbone median - promoter divergence)", "backbone_distal_matched"),
      "fig4_backbone_effect_by_group", 7.5, 6)

# ---- Fig 3: promoter divergence vs matched-control divergence
d3 <- res[res$primary_status == "estimable", ]
lim <- max(c(d3$promoter_divergence, d3$control_median_divergence), na.rm = TRUE)
p3 <- ggplot(d3, aes(control_median_divergence, promoter_divergence, colour = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) + geom_point(alpha = 0.4, size = 1) +
  scale_colour_manual(values = pal, labels = group_labels, name = NULL) + coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = "Median divergence of composition-matched control windows", y = "Promoter divergence (columns 78-127)",
       title = "Within-array promoter vs matched-control divergence", subtitle = "Below the diagonal: promoter less divergent than its matched controls") +
  theme(legend.position = "top")
p3b <- ggplot(d3, aes(control_median_divergence - promoter_divergence, fill = group)) + geom_histogram(bins = 50, alpha = 0.8, position = "identity") +
  geom_vline(xintercept = 0, linetype = 2) + scale_fill_manual(values = pal, guide = "none") +
  facet_wrap(~group, ncol = 1, scales = "free_y", labeller = as_labeller(group_labels)) +
  labs(x = "Effect", y = "Arrays")
save2(p3 + p3b + plot_layout(widths = c(2, 1.2)), "fig3_promoter_vs_matched_control", 12, 6.5)

# ---- Fig 5: QC
mk <- function(v, nm) data.frame(metric = nm, value = v)
p5a <- ggplot(qc, aes(n_anchors)) + geom_histogram(binwidth = 1, fill = "#2E86C1") + labs(x = "L1PA10 anchors per array", y = "Arrays", title = "Anchors")
p5b <- ggplot(qc, aes(n_units_complete)) + geom_histogram(binwidth = 1, fill = "#27AE60") + labs(x = "Complete units (3.5-20 kb)", y = "Arrays", title = "Complete units")
p5c <- ggplot(qc, aes(fasta_length / 1000)) + geom_histogram(bins = 50, fill = "#7D3C98") + scale_x_log10() + labs(x = "Input array length (kb, log)", y = "Arrays", title = "Array length")
ms <- merge(qc, res[, c("array_id", "recombination_class")], by = "array_id")
p5d <- ggplot(ms, aes(recombination_class, fill = ifelse(mixed_strand, "mixed-strand", "single strand"))) + geom_bar(position = "stack") +
  scale_fill_manual(values = c("mixed-strand" = "#E74C3C", "single strand" = "grey70"), name = NULL) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x)) + labs(x = NULL, y = "Arrays", title = "Strand structure") +
  theme(legend.position = "bottom")
ne2 <- if (nrow(ne)) as.data.frame(table(analysis = ne$analysis, reason = ne$reason_code)) else data.frame(analysis = character(), reason = character(), Freq = integer())
p5e <- if (nrow(ne2)) ggplot(ne2, aes(reorder(reason, Freq), Freq, fill = analysis)) + geom_col(position = "dodge") + coord_flip() +
  labs(x = NULL, y = "Arrays", title = "Non-estimability reasons") + theme(legend.position = "bottom") else ggplot() + theme_void() + labs(title = "No non-estimable arrays")
p5f <- ggplot(res, aes(n_units_irregular)) + geom_histogram(binwidth = 1, fill = "#CA6F1E") + labs(x = "Irregular adjacent units per array", y = "Arrays", title = "Irregular units")
save2((p5a | p5b | p5c) / (p5d | p5e | p5f), "fig5_qc_overview", 14, 8)

# ---- Fig 6: rasterised heatmap (arrays x positions), sorted by group then effect
hm <- read.delim(gzfile(file.path(comb, "array_by_position_rolling11_divergence.tsv.gz")), check.names = FALSE, stringsAsFactors = FALSE)
ord <- res[order(res$group, ifelse(is.na(res$effect), -Inf, -res$effect), res$array_id), c("array_id", "group", "effect", "primary_status")]
ord$row <- seq_len(nrow(ord))
hm <- hm[match(ord$array_id, hm$array_id), ]
mat <- as.matrix(hm[, -1]); n <- nrow(mat)
long <- data.frame(row = rep(ord$row, times = ncol(mat)), col = rep(seq_len(ncol(mat)), each = n), v = as.vector(mat))
long <- long[is.finite(long$v), ]
cap <- quantile(long$v, 0.99, na.rm = TRUE)
h <- ggplot(long, aes(col, row, fill = pmin(v, cap))) + geom_raster() +
  geom_vline(xintercept = c(P0 - 0.5, P1 + 0.5), colour = "black", linewidth = 0.3) +
  scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "11-col\ndivergence\n(capped at p99)", na.value = "white") +
  labs(x = "L1PA10 alignment column (promoter between black lines)", y = NULL,
       title = sprintf("Per-array divergence across the 300-column alignment (%d arrays)", n),
       subtitle = "Rows sorted by group, then by decreasing promoter effect; blank = not callable / non-estimable") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
side <- ggplot(ord, aes(1, row, fill = group)) + geom_raster() + scale_y_reverse(expand = c(0, 0)) + scale_fill_manual(values = pal, name = "Group") +
  theme_void() + theme(legend.position = "left")
ggsave(file.path(figd, "fig6_array_divergence_heatmap.png"), side + h + plot_layout(widths = c(0.03, 1)), width = 11, height = 11, dpi = 200, device = ragg::agg_png)
ggsave(file.path(figd, "fig6_array_divergence_heatmap.pdf"), side + h + plot_layout(widths = c(0.03, 1)),
       width = 11, height = 11, device = grDevices::pdf, useDingbats = FALSE)
cat("figures written to", figd, "\n")
