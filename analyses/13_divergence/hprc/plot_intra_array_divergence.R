#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
input_path <- file.path(root, "tables", "full_hprc", "array_cluster_mash_divergence.tsv")
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

minimum_units <- 5L
h_levels <- c("H1", "H2", "H3")
h_palette <- c(H1 = "#4E79A7", H2 = "#F28E2B", H3 = "#59A14F")

all_values <- fread(input_path)
plot_values <- all_values[cluster_units >= minimum_units]
plot_values[, h_cluster := factor(h_cluster, levels = h_levels)]
if (uniqueN(plot_values$h_cluster) != 3L) stop("Missing H cluster")

summary_table <- plot_values[, .(
  arrays = .N,
  ACREs = sum(cluster_units),
  pairwise_comparisons = sum(pairwise_comparisons),
  mean_divergence_pct = mean(mean_mash_divergence_pct),
  sd_divergence_pct = sd(mean_mash_divergence_pct),
  median_divergence_pct = median(mean_mash_divergence_pct),
  q25_divergence_pct = quantile(mean_mash_divergence_pct, 0.25),
  q75_divergence_pct = quantile(mean_mash_divergence_pct, 0.75),
  min_divergence_pct = min(mean_mash_divergence_pct),
  max_divergence_pct = max(mean_mash_divergence_pct)
), by = h_cluster]
setorder(summary_table, h_cluster)

comparisons <- list(c("H1", "H2"), c("H1", "H3"), c("H2", "H3"))
rank_tests <- rbindlist(lapply(comparisons, function(pair) {
  x <- plot_values[h_cluster == pair[[1]], mean_mash_divergence_pct]
  y <- plot_values[h_cluster == pair[[2]], mean_mash_divergence_pct]
  test <- wilcox.test(x, y, exact = FALSE)
  # R reports the Mann-Whitney U statistic as W for the two-sample test.
  u_value <- as.numeric(test$statistic)
  data.table(
    comparison = paste(pair, collapse = " vs "),
    cluster_1 = pair[[1]], cluster_2 = pair[[2]],
    arrays_1 = length(x), arrays_2 = length(y),
    median_1_pct = median(x), median_2_pct = median(y),
    median_difference_2_minus_1_pct = median(y) - median(x),
    wilcoxon_W = u_value,
    p_value = test$p.value,
    cliffs_delta_cluster_1_minus_2 = 2 * u_value / (length(x) * length(y)) - 1
  )
}))
rank_tests[, BH_q_value := p.adjust(p_value, method = "BH")]

paired_tests <- rbindlist(lapply(comparisons, function(pair) {
  wide <- dcast(
    plot_values[h_cluster %chin% pair],
    array_id + source_name ~ h_cluster,
    value.var = "mean_mash_divergence_pct"
  )
  wide <- wide[complete.cases(wide[, ..pair])]
  if (nrow(wide) < 2L) {
    return(data.table(
      comparison = paste(pair, collapse = " vs "),
      paired_arrays = nrow(wide),
      median_difference_2_minus_1_pct = if (nrow(wide)) wide[[pair[[2]]]] - wide[[pair[[1]]]] else NA_real_,
      paired_wilcoxon_V = NA_real_, p_value = NA_real_
    ))
  }
  differences <- wide[[pair[[2]]]] - wide[[pair[[1]]]]
  test <- wilcox.test(
    wide[[pair[[2]]]], wide[[pair[[1]]]],
    paired = TRUE, exact = FALSE
  )
  data.table(
    comparison = paste(pair, collapse = " vs "),
    paired_arrays = nrow(wide),
    median_difference_2_minus_1_pct = median(differences),
    paired_wilcoxon_V = as.numeric(test$statistic),
    p_value = test$p.value
  )
}))
paired_tests[, BH_q_value := p.adjust(p_value, method = "BH")]

fwrite(summary_table, file.path(table_dir, "HPRC_intra_array_divergence_cluster_summary.tsv"), sep = "\t")
fwrite(rank_tests, file.path(table_dir, "HPRC_intra_array_divergence_rank_tests.tsv"), sep = "\t")
fwrite(paired_tests, file.path(table_dir, "HPRC_intra_array_divergence_paired_tests.tsv"), sep = "\t")
fwrite(plot_values, file.path(table_dir, "HPRC_intra_array_divergence_plotted_arrays.tsv.gz"), sep = "\t", compress = "gzip")

counts <- setNames(summary_table$arrays, as.character(summary_table$h_cluster))
x_labels <- setNames(
  sprintf("%s\n(n=%s arrays)", h_levels, comma(counts[h_levels])),
  h_levels
)
median_labels <- summary_table[, .(
  h_cluster,
  median_divergence_pct,
  label = sprintf("%.2f%%", median_divergence_pct)
)]

theme_fig <- theme_bw(base_size = 16, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 23, face = "bold", colour = "#111111"),
    plot.subtitle = element_text(size = 13, colour = "#4B5563", margin = margin(b = 8)),
    axis.title = element_text(size = 18, face = "bold", colour = "#111111"),
    axis.text = element_text(size = 14, colour = "#252525"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#E4E7EB", linewidth = 0.45),
    plot.margin = margin(16, 18, 14, 16)
  )

plot <- ggplot(
  plot_values,
  aes(h_cluster, mean_mash_divergence_pct, fill = h_cluster, colour = h_cluster)
) +
  geom_violin(width = 0.86, trim = TRUE, alpha = 0.50, linewidth = 0.75) +
  geom_point(
    position = position_jitter(width = 0.13, height = 0, seed = 20261108),
    size = 1.05, alpha = 0.22, stroke = 0
  ) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, fill = "white",
    colour = "#222222", linewidth = 0.60
  ) +
  geom_point(
    data = median_labels,
    aes(h_cluster, median_divergence_pct),
    shape = 23, size = 4.1, fill = "#F2C14E", colour = "#111111",
    stroke = 0.75, inherit.aes = FALSE
  ) +
  geom_text(
    data = median_labels,
    aes(h_cluster, median_divergence_pct, label = label),
    nudge_y = 0.18, family = "Helvetica", fontface = "bold", size = 4.2,
    colour = "#111111", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = h_palette, guide = "none") +
  scale_colour_manual(values = h_palette, guide = "none") +
  scale_x_discrete(labels = x_labels, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 5, 1),
    labels = label_number(accuracy = 0.1, suffix = "%"),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    title = "Intra-array ACRE sequence divergence",
    subtitle = paste0(
      "Array-level mean of within-cluster pairwise Mash distances\n",
      "Arrays with at least five ACREs per cluster; diamonds indicate medians"
    ),
    x = "NTRiage cluster",
    y = "Mean within-array pairwise\nMash divergence"
  ) +
  theme_fig

png_path <- file.path(figure_dir, "Supplementary_Figure_9_HPRC_intra_array_divergence.png")
pdf_path <- file.path(figure_dir, "Supplementary_Figure_9_HPRC_intra_array_divergence.pdf")
ggsave(
  png_path, plot, width = 8.4, height = 6.6, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  pdf_path, plot, width = 8.4, height = 6.6, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

cat("Created:", pdf_path, "\n")
cat("Created:", png_path, "\n")
print(summary_table)
print(rank_tests)
print(paired_tests)
