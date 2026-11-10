#!/usr/bin/env Rscript

options(warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
out_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(out_root, ".."), mustWork = TRUE)

data_dir <- file.path(out_root, "data")
panel_dir <- file.path(out_root, "panels")
figure_dir <- file.path(out_root, "figures")
intermediate_dir <- file.path(out_root, "intermediate")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

# Additional outputs requested after the original Figure 4 assembly. These are
# written to a separate folder so the curated split-only figure remains intact.
all_arrays_root <- file.path(project_root, "figure_5_all_arrays_alt")
all_arrays_data_dir <- file.path(all_arrays_root, "data")
all_arrays_panel_dir <- file.path(all_arrays_root, "panels")
all_arrays_figure_dir <- file.path(all_arrays_root, "figures")
dir.create(all_arrays_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(all_arrays_panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(all_arrays_figure_dir, recursive = TRUE, showWarnings = FALSE)

ntriage_root <- file.path(
  project_root, "NTRIAGE_OCT26", "manuscript_output", "telomeric_ntriage_2794"
)
round2_root <- file.path(
  project_root, "hprc_breakpoint_ntriage_cluster_analysis",
  "custom_manual_shift_breakpoint_map_round2"
)
nov3_root <- file.path(project_root, "nov3_ntriage_analysis")

pca_path <- file.path(ntriage_root, "tables", "unit_pca_metadata.tsv.gz")
pca_variance_path <- file.path(ntriage_root, "tables", "pca_explained_variance.tsv")
current_units_path <- file.path(
  ntriage_root, "native", "telomeric_2794_core.clusters.tsv"
)
audit_path <- file.path(round2_root, "tables", "round2_row_audit_1196.tsv")
heat_path <- file.path(round2_root, "tables", "round2_units_plus_minus_8_1196.tsv")
nov3_units_path <- file.path(
  nov3_root, "runs", "hprc", "native", "hprc.clusters.tsv"
)
assignments_path <- file.path(
  nov3_root, "results", "hprc", "k_selection", "assignments.tsv.gz"
)
hor_path <- file.path(data_dir, "HPRC_HOR_scores_all_thresholds.tsv.gz")

input_paths <- c(
  pca_path, pca_variance_path, current_units_path, audit_path, heat_path,
  nov3_units_path, assignments_path, hor_path
)
for (path in input_paths) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

h_levels <- c("H1", "H2", "H3")
h_palette <- c(H1 = "#4E79A7", H2 = "#F28E2B", H3 = "#59A14F")

# The November HPRC analysis and the October manuscript run use different ACRE
# numbering in a small number of arrays. Map by source_name + unit_in_array,
# then resolve the few segmentation-shifted terminal units by maximal genomic
# overlap. This avoids silently joining unrelated ACRE IDs across runs.
assignments <- fread(cmd = sprintf("gzip -dc %s", shQuote(assignments_path)))
assignments[, h_cluster := fcase(
  selected_cluster == 1L, "H1",
  selected_cluster == 2L, "H2",
  selected_cluster == 3L, "H3",
  default = NA_character_
)]
if (anyNA(assignments$h_cluster)) stop("Unexpected selected_cluster in HPRC k=3 assignments")

nov3_units <- fread(nov3_units_path)
nov3_units <- merge(
  nov3_units,
  assignments[, .(acre_id, h_cluster)],
  by = "acre_id", all.x = TRUE, sort = FALSE
)
if (anyNA(nov3_units$h_cluster)) stop("Missing H1-H3 assignment for November HPRC unit")
if (uniqueN(nov3_units, by = c("source_name", "unit_in_array")) != nrow(nov3_units)) {
  stop("November source_name + unit_in_array keys are not unique")
}

stable_map <- nov3_units[, .(source_name, unit_in_array, h_cluster)]

resolve_missing_by_overlap <- function(frame, start_col, end_col, label) {
  missing <- which(is.na(frame$h_cluster))
  if (!length(missing)) return(frame)
  for (idx in missing) {
    source_value <- frame$source_name[[idx]]
    start_value <- frame[[start_col]][[idx]]
    end_value <- frame[[end_col]][[idx]]
    candidates <- nov3_units[source_name == source_value]
    overlaps <- pmax(
      0,
      pmin(end_value, candidates$end) - pmax(start_value, candidates$start)
    )
    if (!length(overlaps) || max(overlaps) <= 0) {
      stop("Could not overlap-map ", label, " row ", idx, ": ", source_value)
    }
    winner <- which.max(overlaps)
    frame$h_cluster[[idx]] <- candidates$h_cluster[[winner]]
    frame$mapping_method[[idx]] <- "maximal_coordinate_overlap"
    frame$mapping_overlap_bp[[idx]] <- overlaps[[winner]]
    frame$mapped_nov3_acre_id[[idx]] <- candidates$acre_id[[winner]]
  }
  frame
}

audit <- fread(audit_path)
heat <- fread(heat_path)
pca_variance <- fread(pca_variance_path)
pca_all <- fread(cmd = sprintf("gzip -dc %s", shQuote(pca_path)))
current_units <- fread(current_units_path)

stopifnot(nrow(audit) == 1196L, uniqueN(audit$source_name) == 1196L)
stopifnot(uniqueN(heat$source_name) == 1196L)

pca_all <- merge(
  pca_all,
  current_units[, .(
    acre_id, current_start = start, current_end = end,
    current_native_cluster = cluster
  )],
  by = "acre_id", all.x = TRUE, sort = FALSE
)
if (anyNA(pca_all$current_start)) stop("Current PCA units did not map to current coordinates")
pca_all <- merge(
  pca_all, stable_map,
  by = c("source_name", "unit_in_array"), all.x = TRUE, sort = FALSE
)
pca_all[, `:=`(
  mapping_method = fifelse(!is.na(h_cluster), "source_name_plus_unit_index", NA_character_),
  mapping_overlap_bp = NA_real_,
  mapped_nov3_acre_id = NA_character_
)]
pca_all <- resolve_missing_by_overlap(
  pca_all, "current_start", "current_end", "all-array PCA"
)
if (anyNA(pca_all$h_cluster)) stop("Unresolved all-array PCA H cluster")
if (uniqueN(pca_all$source_name) != 2794L || nrow(pca_all) != 77276L) {
  stop(
    "Expected 77,276 PCA units from 2,794 HPRC arrays; found ",
    nrow(pca_all), " units from ", uniqueN(pca_all$source_name), " arrays"
  )
}

pca <- pca_all[source_name %chin% audit$source_name]
if (uniqueN(pca$source_name) != 1196L) stop("PCA subset did not recover all 1,196 arrays")

heat <- merge(
  heat, stable_map,
  by = c("source_name", "unit_in_array"), all.x = TRUE, sort = FALSE
)
heat[, `:=`(
  mapping_method = fifelse(!is.na(h_cluster), "source_name_plus_unit_index", NA_character_),
  mapping_overlap_bp = NA_real_,
  mapped_nov3_acre_id = NA_character_
)]
heat <- resolve_missing_by_overlap(heat, "start", "end", "transition map")
if (anyNA(heat$h_cluster)) stop("Unresolved transition-map H cluster")
setorder(heat, round2_row, units_from_round2_anchor)

pc1_var <- pca_variance[principal_component == 1, explained_variance_ratio][[1]]
pc2_var <- pca_variance[principal_component == 2, explained_variance_ratio][[1]]

pca[, h_cluster := factor(h_cluster, levels = h_levels)]
pca[, cluster_units := .N, by = h_cluster]
setorder(pca, -cluster_units, h_cluster)
pca_counts <- pca[, .(units = .N, arrays = uniqueN(source_name)), by = h_cluster]
pca_centroids <- pca[, .(
  PC1 = median(PC1), PC2 = median(PC2),
  units = .N, arrays = uniqueN(source_name)
), by = h_cluster]
pca_legend_labels <- setNames(
  sprintf("%s  (n=%s)", h_levels, format(pca_counts$units[match(h_levels, pca_counts$h_cluster)], big.mark = ",")),
  h_levels
)

theme_fig3 <- function(base_size = 15) {
  theme_bw(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 23, face = "bold", colour = "#111111", hjust = 0),
      plot.subtitle = element_text(size = 13, colour = "#4B5563", hjust = 0),
      plot.tag = element_text(size = 24, face = "plain", colour = "#111111"),
      plot.tag.position = c(0.01, 0.99),
      axis.title = element_text(size = 18, face = "bold", colour = "#111111"),
      axis.text = element_text(size = 14, colour = "#252525"),
      legend.title = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 15),
      legend.key.height = grid::unit(0.62, "cm"),
      legend.key.width = grid::unit(0.70, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E4E7EB", linewidth = 0.42),
      plot.margin = margin(12, 14, 11, 12)
    )
}

panel_a <- ggplot(pca, aes(PC1, PC2, colour = h_cluster)) +
  geom_point(size = 0.72, alpha = 0.42, stroke = 0) +
  geom_point(
    data = pca_centroids,
    aes(PC1, PC2, fill = h_cluster),
    shape = 23, size = 5.3, stroke = 1.0,
    colour = "#111111", inherit.aes = FALSE
  ) +
  geom_text(
    data = pca_centroids,
    aes(PC1, PC2, label = h_cluster),
    colour = "#111111", fontface = "bold", size = 5.4,
    nudge_y = 0.62, show.legend = FALSE, inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = h_palette, breaks = h_levels, labels = pca_legend_labels,
    drop = FALSE
  ) +
  scale_fill_manual(values = h_palette, guide = "none") +
  labs(
    tag = "B",
    title = "Human ACRO-unit PCA",
    subtitle = sprintf(
      "%s complete units from 1,196 switch-positive arrays",
      format(nrow(pca), big.mark = ",")
    ),
    x = sprintf("PC1 (%.1f%% variance)", 100 * pc1_var),
    y = sprintf("PC2 (%.1f%% variance)", 100 * pc2_var),
    colour = "NTRiage cluster"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5.2, alpha = 1))) +
  theme_fig3(15) +
  theme(legend.position = "right")

# Alternate panel A: the same H1-H3 remapping applied to every complete ACRE
# from all 2,794 HPRC arrays, rather than only the curated split-positive set.
pca_all[, h_cluster := factor(h_cluster, levels = h_levels)]
pca_all[, cluster_units := .N, by = h_cluster]
setorder(pca_all, -cluster_units, h_cluster)
pca_all_counts <- pca_all[, .(
  units = .N,
  arrays = uniqueN(source_name)
), by = h_cluster]
pca_all_centroids <- pca_all[, .(
  PC1 = median(PC1),
  PC2 = median(PC2),
  units = .N,
  arrays = uniqueN(source_name)
), by = h_cluster]
pca_all_legend_labels <- setNames(
  sprintf(
    "%s  (n=%s)",
    h_levels,
    format(
      pca_all_counts$units[match(h_levels, pca_all_counts$h_cluster)],
      big.mark = ","
    )
  ),
  h_levels
)

panel_a_all <- ggplot(pca_all, aes(PC1, PC2, colour = h_cluster)) +
  geom_point(size = 0.48, alpha = 0.25, stroke = 0) +
  geom_point(
    data = pca_all_centroids,
    aes(PC1, PC2, fill = h_cluster),
    shape = 23, size = 5.3, stroke = 1.0,
    colour = "#111111", inherit.aes = FALSE
  ) +
  geom_text(
    data = pca_all_centroids,
    aes(PC1, PC2, label = h_cluster),
    colour = "#111111", fontface = "bold", size = 5.4,
    nudge_y = 0.62, show.legend = FALSE, inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = h_palette, breaks = h_levels, labels = pca_all_legend_labels,
    drop = FALSE
  ) +
  scale_fill_manual(values = h_palette, guide = "none") +
  labs(
    tag = "B",
    title = "Human ACRO-unit PCA",
    subtitle = sprintf(
      "%s complete units from all 2,794 HPRC arrays",
      format(nrow(pca_all), big.mark = ",")
    ),
    x = sprintf("PC1 (%.1f%% variance)", 100 * pc1_var),
    y = sprintf("PC2 (%.1f%% variance)", 100 * pc2_var),
    colour = "NTRiage cluster"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5.2, alpha = 1))) +
  theme_fig3(15) +
  theme(legend.position = "right")

# The full round-2 audit contains 1,196 arrays. Fourteen edge-unavailable calls
# have no localized internal position, leaving 1,182 arrays for the first two
# violin summaries.
array_unit_counts <- pca[, .(complete_units_n = max(unit_in_array)), by = source_name]
violin_arrays <- merge(
  audit, array_unit_counts,
  by = "source_name", all.x = TRUE, sort = FALSE
)
if (anyNA(violin_arrays$complete_units_n)) stop("Missing complete-unit count")
violin_arrays[, transition_position_fraction := round2_units_left / complete_units_n]
violin_arrays <- violin_arrays[
  round2_transition != "edge-unavailable" &
    round2_units_left > 0 & round2_units_left < complete_units_n
]
if (nrow(violin_arrays) != 1182L) {
  stop("Expected 1,182 localized arrays, found ", nrow(violin_arrays))
}

localized_sources <- violin_arrays$source_name
acre_length_plot <- copy(nov3_units[source_name %chin% localized_sources])
acre_length_plot[, `:=`(
  h_cluster = factor(h_cluster, levels = h_levels),
  acre_length_kb = length / 1000
)]
if (
  uniqueN(acre_length_plot$source_name) != 1182L ||
    any(!is.finite(acre_length_plot$acre_length_kb)) ||
    any(acre_length_plot$acre_length_kb <= 0)
) {
  stop("Invalid localized-array ACRE-length plotting subset")
}
acre_length_summary <- acre_length_plot[, .(
  units = .N,
  arrays = uniqueN(source_name),
  mean_length_kb = mean(acre_length_kb),
  median_length_kb = median(acre_length_kb),
  q25_length_kb = quantile(acre_length_kb, 0.25),
  q75_length_kb = quantile(acre_length_kb, 0.75)
), by = h_cluster]

hor <- fread(cmd = sprintf("gzip -dc %s", shQuote(hor_path)))
# The HOR table was generated before the Figure 4-compatible H2/H3 naming
# correction. Swap those two stored labels so every panel uses one convention.
hor[, cluster := fcase(
  cluster == "H1", "H1",
  cluster == "H2", "H3",
  cluster == "H3", "H2",
  default = NA_character_
)]
if (anyNA(hor$cluster)) stop("Unexpected H cluster in HOR score table")
hor <- merge(
  hor,
  nov3_units[, .(acre_id, source_name)],
  by = "acre_id", all.x = TRUE, sort = FALSE
)
if (anyNA(hor$source_name)) stop("HOR rows did not map to November source arrays")
hor_plot <- hor[source_name %chin% audit$source_name]
hor_plot[, h_cluster := factor(cluster, levels = h_levels)]
hor_plot[, hor_score_90 := score_0.1]
if (uniqueN(hor_plot$source_name) != 1196L) stop("HOR subset did not recover 1,196 arrays")

hor_plot_all <- copy(hor)
hor_plot_all[, h_cluster := factor(cluster, levels = h_levels)]
hor_plot_all[, hor_score_90 := score_0.1]
if (nrow(hor_plot_all) != 77351L || uniqueN(hor_plot_all$source_name) != 2794L) {
  stop(
    "Expected 77,351 HOR-scored units from 2,794 arrays; found ",
    nrow(hor_plot_all), " units from ", uniqueN(hor_plot_all$source_name), " arrays"
  )
}
hor_all_means <- hor_plot_all[, .(
  mean_score = mean(hor_score_90)
), by = h_cluster]
hor_all_summary <- hor_plot_all[, .(
  units = .N,
  arrays = uniqueN(source_name),
  mean_score_pct = mean(hor_score_90),
  median_score_pct = median(hor_score_90),
  pct_units_nonzero = 100 * mean(hor_score_90 > 0)
), by = h_cluster]

theme_violin <- function(base_size = 13) {
  theme_bw(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(
        size = 14.5, face = "bold", colour = "#111111", hjust = 0,
        lineheight = 0.94
      ),
      plot.subtitle = element_text(size = 9.8, colour = "#4B5563", hjust = 0),
      plot.tag = element_text(size = 24, face = "plain", colour = "#111111"),
      plot.tag.position = c(0.01, 0.99),
      axis.title = element_text(size = 12.3, face = "bold", colour = "#111111"),
      axis.text = element_text(size = 10.5, colour = "#252525"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "#E4E7EB", linewidth = 0.38),
      plot.margin = margin(10, 8, 9, 10)
    )
}

panel_c_length <- ggplot(
  acre_length_plot,
  aes(x = h_cluster, y = acre_length_kb, fill = h_cluster, color = h_cluster)
) +
  geom_violin(linewidth = 0.65, trim = TRUE, width = 0.88, alpha = 0.72) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, fill = "white",
    color = "#222222", linewidth = 0.48
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 23, size = 2.6,
    fill = "white", color = "#111111", stroke = 0.55
  ) +
  scale_fill_manual(values = h_palette, drop = FALSE) +
  scale_color_manual(values = h_palette, drop = FALSE) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 5),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    tag = "C",
    title = "ACRE length\nby cluster",
    subtitle = sprintf(
      "%s units; 1,182 localized arrays",
      format(nrow(acre_length_plot), big.mark = ",")
    ),
    x = "NTRiage cluster", y = "ACRE length (kb)"
  ) +
  guides(fill = "none", color = "none") +
  theme_violin(13)

panel_b_strength <- ggplot(
  violin_arrays,
  aes(x = factor("Localized transitions"), y = within_minus_cross_ani)
) +
  geom_violin(
    fill = "#86CCF1", color = "#1D395F", linewidth = 0.70,
    width = 0.84, trim = TRUE
  ) +
  geom_boxplot(
    width = 0.14, outlier.shape = NA, fill = "white",
    color = "#1D395F", linewidth = 0.55
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 21, size = 2.6,
    fill = "#F28E2B", color = "white", stroke = 0.55
  ) +
  geom_hline(yintercept = 4, linetype = "22", linewidth = 0.55, color = "#A33F2B") +
  annotate(
    "text", x = 1, y = 4.05, label = "4-point threshold",
    family = "Helvetica", size = 2.8, hjust = 0.5, vjust = -0.25, color = "#8A3528"
  ) +
  labs(
    tag = "B", title = "Transition\nstrength",
    subtitle = "n=1,182 arrays",
    x = NULL,
    y = "Within-block minus cross-boundary ANI\n(percentage points)"
  ) +
  theme_violin(13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

panel_b_position <- ggplot(
  violin_arrays,
  aes(x = factor("Localized transitions"), y = transition_position_fraction)
) +
  geom_violin(
    fill = "#DDEFD9", color = "#2B6B3D", linewidth = 0.70,
    width = 0.84, trim = TRUE
  ) +
  geom_boxplot(
    width = 0.14, outlier.shape = NA, fill = "white",
    color = "#2B6B3D", linewidth = 0.55
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 21, size = 2.6,
    fill = "#F28E2B", color = "white", stroke = 0.55
  ) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Position within\narray",
    subtitle = "n=1,182 arrays",
    x = NULL,
    y = "Fraction of complete ACRO units\nleft of transition"
  ) +
  theme_violin(13) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

panel_b_hor <- ggplot(
  hor_plot,
  aes(x = h_cluster, y = hor_score_90, fill = h_cluster, color = h_cluster)
) +
  geom_violin(linewidth = 0.65, trim = TRUE, width = 0.88, alpha = 0.72) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, fill = "white",
    color = "#222222", linewidth = 0.48
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 23, size = 2.6,
    fill = "white", color = "#111111", stroke = 0.55
  ) +
  scale_fill_manual(values = h_palette, drop = FALSE) +
  scale_color_manual(values = h_palette, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = label_number(suffix = "%"),
    expand = expansion(mult = c(0.005, 0.02))
  ) +
  labs(
    title = "HOR score\nby cluster",
    subtitle = sprintf(
      "%s units; 1,196 arrays; 90%% ID",
      format(nrow(hor_plot), big.mark = ",")
    ),
    x = "NTRiage cluster", y = "HOR score (%)"
  ) +
  guides(fill = "none", color = "none") +
  theme_violin(13)

panel_b <- (panel_b_strength | panel_b_position | panel_b_hor) +
  plot_layout(widths = c(0.92, 0.92, 1.18))

# Corrected all-cohort HOR panel for the alternate main figure. The original
# split-only panel above is retained for reproducibility, but its diamond is a
# median from the selected 1,196 arrays. Here the violin includes all HPRC
# arrays, and the labeled diamond is the pooled per-unit mean.
panel_b_hor_all <- ggplot(
  hor_plot_all,
  aes(x = h_cluster, y = hor_score_90, fill = h_cluster, color = h_cluster)
) +
  geom_violin(linewidth = 0.65, trim = TRUE, width = 0.88, alpha = 0.72) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, fill = "white",
    color = "#222222", linewidth = 0.48
  ) +
  geom_point(
    data = hor_all_means,
    aes(x = h_cluster, y = mean_score),
    shape = 23, size = 3.0, fill = "#F2C14E", color = "#111111",
    stroke = 0.65, inherit.aes = FALSE
  ) +
  geom_text(
    data = hor_all_means,
    aes(x = h_cluster, y = mean_score, label = sprintf("%.1f%%", mean_score)),
    nudge_y = 4.2, family = "Helvetica", fontface = "bold", size = 3.1,
    color = "#111111", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = h_palette, drop = FALSE) +
  scale_color_manual(values = h_palette, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100), breaks = seq(0, 100, 25),
    labels = label_number(suffix = "%"),
    expand = expansion(mult = c(0.005, 0.02))
  ) +
  labs(
    title = "HOR score\nby cluster",
    subtitle = "All 2,794 arrays; 90% ID; diamonds = means",
    x = "NTRiage cluster", y = "HOR score (%)"
  ) +
  guides(fill = "none", color = "none") +
  theme_violin(13)

panel_b_all <- (panel_b_strength | panel_b_position | panel_b_hor_all) +
  plot_layout(widths = c(0.92, 0.92, 1.18))

panel_c_all <- (panel_c_length | panel_b_position | panel_b_hor_all) +
  plot_layout(widths = c(1.12, 0.92, 1.18))

localized_order <- unique(
  violin_arrays[, .(source_name, original_round2_row = round2_row)]
)
setorder(localized_order, original_round2_row)
localized_order[, localized_row := .I]

heat_localized <- merge(
  heat,
  localized_order[, .(source_name, localized_row)],
  by = "source_name", all = FALSE, sort = FALSE
)
heat_localized[, original_round2_row := round2_row]
if (uniqueN(heat_localized$source_name) != 1182L) {
  stop("Transition map did not recover all 1,182 localized arrays")
}
setorder(heat_localized, localized_row, units_from_round2_anchor)
heat_localized[, h_cluster := factor(h_cluster, levels = h_levels)]
heat_counts <- heat_localized[, .(units = .N), by = h_cluster]
heat_legend_labels <- setNames(
  sprintf("%s  (n=%s)", h_levels, format(heat_counts$units[match(h_levels, heat_counts$h_cluster)], big.mark = ",")),
  h_levels
)

panel_a_transition <- ggplot(
  heat_localized,
  aes(units_from_round2_anchor, localized_row, fill = h_cluster)
) +
  geom_tile(width = 0.96, height = 0.96, linewidth = 0) +
  geom_vline(xintercept = 0, linewidth = 1.0, color = "black") +
  scale_fill_manual(
    values = h_palette, breaks = h_levels, labels = heat_legend_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(-8:-1, 1:8), expand = expansion(mult = c(0.008, 0.008))
  ) +
  scale_y_reverse(
    breaks = unique(c(1L, seq(200L, 1000L, 200L), 1182L)),
    expand = expansion(mult = c(0.0015, 0.0015))
  ) +
  labs(
    tag = "A",
    title = NULL,
    subtitle = "n = 1,182",
    x = "ACRO units from transition anchor",
    y = "Localized array row",
    fill = "NTRiage cluster"
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme_fig3(14) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 12.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10.5),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Standalone supplementary analogue of panel C for the arrays outside the
# curated split set. Because these arrays have no transition anchor, their 16
# central units are aligned at the array midpoint. The complement is retained
# in full, including isolated or mixed labels that did not yield a curated
# transition, so the plot does not overstate uniformity.
non_split_sources <- setdiff(unique(pca_all$source_name), audit$source_name)
if (length(non_split_sources) != 1598L) {
  stop("Expected 1,598 arrays outside the curated split set, found ", length(non_split_sources))
}

non_split_units <- copy(pca_all[source_name %chin% non_split_sources])
setorder(non_split_units, source_name, unit_in_array)
non_split_array_summary <- non_split_units[, .(
  array_n = max(unit_in_array),
  H1_units = sum(h_cluster == "H1"),
  H2_units = sum(h_cluster == "H2"),
  H3_units = sum(h_cluster == "H3"),
  unique_clusters = uniqueN(h_cluster),
  adjacent_label_changes = sum(
    as.character(h_cluster) != shift(as.character(h_cluster)),
    na.rm = TRUE
  )
), by = source_name]

dominant_count_matrix <- as.matrix(
  non_split_array_summary[, .(H1_units, H2_units, H3_units)]
)
non_split_array_summary[, dominant_cluster := h_levels[
  max.col(dominant_count_matrix, ties.method = "first")
]]
non_split_array_summary[, dominant_fraction := pmax(H1_units, H2_units, H3_units) / array_n]
non_split_array_summary[, dominant_order := match(dominant_cluster, h_levels)]
setorder(
  non_split_array_summary,
  dominant_order, -dominant_fraction, adjacent_label_changes, -array_n, source_name
)
non_split_array_summary[, non_split_row := .I]

non_split_units <- merge(
  non_split_units,
  non_split_array_summary[, .(
    source_name, array_n, dominant_cluster, dominant_fraction,
    unique_clusters, adjacent_label_changes, non_split_row
  )],
  by = "source_name", all.x = TRUE, sort = FALSE
)
non_split_units[, midpoint_left_unit := floor(array_n / 2)]
non_split_units[, units_from_midpoint := fifelse(
  unit_in_array <= midpoint_left_unit,
  unit_in_array - midpoint_left_unit - 1L,
  unit_in_array - midpoint_left_unit
)]
non_split_heat <- non_split_units[
  units_from_midpoint >= -8L & units_from_midpoint <= 8L &
    units_from_midpoint != 0L
]
non_split_heat[, midpoint_plot_index := fifelse(
  units_from_midpoint < 0L,
  units_from_midpoint + 9L,
  units_from_midpoint + 8L
)]
setorder(non_split_heat, non_split_row, units_from_midpoint)

non_split_counts <- non_split_heat[, .(units = .N), by = h_cluster]
non_split_legend_labels <- setNames(
  sprintf(
    "%s  (n=%s)",
    h_levels,
    format(
      non_split_counts$units[match(h_levels, non_split_counts$h_cluster)],
      big.mark = ","
    )
  ),
  h_levels
)

supp_non_split <- ggplot(
  non_split_heat,
  aes(midpoint_plot_index, non_split_row, fill = h_cluster)
) +
  geom_tile(width = 0.96, height = 0.96, linewidth = 0) +
  scale_fill_manual(
    values = h_palette, breaks = h_levels, labels = non_split_legend_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:16,
    labels = c(-8:-1, 1:8),
    expand = expansion(mult = c(0.008, 0.008))
  ) +
  scale_y_reverse(
    breaks = unique(c(1L, seq(200L, 1400L, 200L), 1598L)),
    expand = expansion(mult = c(0.0015, 0.0015))
  ) +
  labs(
    title = "H1-H3 clusters in arrays without a curated split",
    subtitle = paste0(
      "1,598 arrays outside the curated split set; 16 central units aligned ",
      "at the array midpoint"
    ),
    x = "ACRO units from array midpoint",
    y = "Non-split array row",
    fill = "NTRiage cluster"
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme_fig3(14) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 12.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10.5),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid = element_blank(),
    legend.position = "right"
  )

right_column <- panel_a_all / panel_c_all +
  plot_layout(heights = c(1.06, 0.94))
figure5_redesign <- (panel_a_transition | right_column) +
  plot_layout(widths = c(1.0, 1.42))

# Retain the historical output paths as aliases of the redesigned main figure.
figure4_alt <- figure5_redesign
figure5_all_arrays_alt <- figure5_redesign

save_plot_pair <- function(plot, stem, width, height, dpi = 400) {
  ggsave(
    file.path(panel_dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in", dpi = dpi,
    device = ragg::agg_png, bg = "white"
  )
  ggsave(
    file.path(panel_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in",
    device = grDevices::pdf, useDingbats = FALSE, bg = "white"
  )
}

save_plot_pair(panel_a_transition, "Figure5A_H1_H3_transition_map_1182", 7.4, 10.4)
save_plot_pair(panel_a_all, "Figure5B_H1_H3_PCA_all_2794_arrays", 8.2, 5.4)
save_plot_pair(panel_c_all, "Figure5C_ACRE_length_position_HOR_violins", 9.2, 5.3)

ggsave(
  file.path(all_arrays_panel_dir, "Figure5B_H1_H3_PCA_all_2794_arrays.png"),
  panel_a_all,
  width = 8.2, height = 5.4, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  file.path(all_arrays_panel_dir, "Figure5B_H1_H3_PCA_all_2794_arrays.pdf"),
  panel_a_all,
  width = 8.2, height = 5.4, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)
ggsave(
  file.path(all_arrays_panel_dir, "Figure5C_ACRE_length_position_HOR_violins.png"),
  panel_c_all,
  width = 9.2, height = 5.3, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  file.path(all_arrays_panel_dir, "Figure5C_ACRE_length_position_HOR_violins.pdf"),
  panel_c_all,
  width = 9.2, height = 5.3, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

supp_non_split_png <- file.path(
  all_arrays_figure_dir,
  "Supplementary_Figure_HPRC_non_split_H1_H3_midpoint_map.png"
)
supp_non_split_pdf <- file.path(
  all_arrays_figure_dir,
  "Supplementary_Figure_HPRC_non_split_H1_H3_midpoint_map.pdf"
)
ggsave(
  supp_non_split_png, supp_non_split,
  width = 7.8, height = 11.0, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  supp_non_split_pdf, supp_non_split,
  width = 7.8, height = 11.0, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

final_png <- file.path(figure_dir, "Figure4_alt_H1_H3_remap.png")
final_pdf <- file.path(figure_dir, "Figure4_alt_H1_H3_remap.pdf")
ggsave(
  final_png, figure4_alt,
  width = 16.0, height = 10.5, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  final_pdf, figure4_alt,
  width = 16.0, height = 10.5, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

redesign_png <- file.path(figure_dir, "Figure5_HPRC_H1_H3_redesign.png")
redesign_pdf <- file.path(figure_dir, "Figure5_HPRC_H1_H3_redesign.pdf")
ggsave(
  redesign_png, figure5_redesign,
  width = 16.0, height = 10.5, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  redesign_pdf, figure5_redesign,
  width = 16.0, height = 10.5, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

all_arrays_final_png <- file.path(
  all_arrays_figure_dir,
  "Figure5_alt_H1_H3_all_arrays_PCA.png"
)
all_arrays_final_pdf <- file.path(
  all_arrays_figure_dir,
  "Figure5_alt_H1_H3_all_arrays_PCA.pdf"
)
ggsave(
  all_arrays_final_png, figure5_all_arrays_alt,
  width = 16.0, height = 10.5, units = "in", dpi = 400,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  all_arrays_final_pdf, figure5_all_arrays_alt,
  width = 16.0, height = 10.5, units = "in",
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)

# Exact plotted subsets and mapping audit.
fwrite(
  pca[, .(
    acre_id, array_id, unit_in_array, source_name, h_cluster,
    mapping_method, mapping_overlap_bp, mapped_nov3_acre_id,
    PC1, PC2, PC3, PC4, PC5
  )],
  file.path(data_dir, "Figure4A_H1_H3_PCA_units_1196.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  pca_centroids,
  file.path(data_dir, "Figure4A_H1_H3_PCA_centroids.tsv"), sep = "\t"
)
fwrite(
  pca_all[, .(
    acre_id, array_id, unit_in_array, source_name, h_cluster,
    mapping_method, mapping_overlap_bp, mapped_nov3_acre_id,
    PC1, PC2, PC3, PC4, PC5
  )],
  file.path(all_arrays_data_dir, "Figure5B_H1_H3_PCA_units_all_2794.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  pca_all_centroids,
  file.path(all_arrays_data_dir, "Figure5B_H1_H3_PCA_centroids_all_2794.tsv"),
  sep = "\t"
)
fwrite(
  hor_all_summary,
  file.path(all_arrays_data_dir, "Figure5C_HOR_all_arrays_cluster_summary.tsv"),
  sep = "\t"
)
fwrite(
  violin_arrays[, .(
    source_name, round2_row, round2_transition, within_minus_cross_ani,
    round2_units_left, complete_units_n, transition_position_fraction,
    breakpoint_local_bp, coordinate_source
  )],
  file.path(data_dir, "Figure4B_transition_strength_position_1182.tsv"),
  sep = "\t"
)
fwrite(
  acre_length_plot[, .(
    acre_id, array_id, unit_in_array, source_name, h_cluster,
    start, end, length, acre_length_kb
  )],
  file.path(data_dir, "Figure5C_ACRE_lengths_H1_H3_1182_arrays.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  acre_length_summary,
  file.path(data_dir, "Figure5C_ACRE_length_cluster_summary_1182_arrays.tsv"),
  sep = "\t"
)
fwrite(
  hor_plot[, .(
    acre_id, source_name, array_id, unit_in_array, h_cluster,
    array_n, hor_score_90
  )],
  file.path(data_dir, "Figure4B_HOR_scores_H1_H3_1196_arrays.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  hor_plot_all[, .(
    acre_id, source_name, array_id, unit_in_array, h_cluster,
    array_n, hor_score_90
  )],
  file.path(data_dir, "Figure4B_HOR_scores_H1_H3_all_2794_arrays.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  heat_localized[, .(
    source_name, acre_id, unit_in_array, start, end,
    units_from_round2_anchor, localized_row, original_round2_row, h_cluster,
    mapping_method, mapping_overlap_bp, mapped_nov3_acre_id,
    coordinate_source, breakpoint_local_bp
  )],
  file.path(data_dir, "Figure5A_H1_H3_units_plus_minus_8_1182.tsv"),
  sep = "\t"
)
fwrite(
  non_split_heat[, .(
    source_name, acre_id, array_id, unit_in_array,
    current_start, current_end, units_from_midpoint, midpoint_plot_index,
    non_split_row,
    h_cluster, dominant_cluster, dominant_fraction,
    unique_clusters, adjacent_label_changes,
    mapping_method, mapping_overlap_bp, mapped_nov3_acre_id
  )],
  file.path(
    all_arrays_data_dir,
    "Supplementary_non_split_H1_H3_midpoint_units_1598.tsv.gz"
  ),
  sep = "\t", compress = "gzip"
)
fwrite(
  non_split_array_summary[, .(
    source_name, non_split_row, array_n,
    H1_units, H2_units, H3_units,
    dominant_cluster, dominant_fraction,
    unique_clusters, adjacent_label_changes
  )],
  file.path(
    all_arrays_data_dir,
    "Supplementary_non_split_H1_H3_array_summary_1598.tsv"
  ),
  sep = "\t"
)

mapping_audit <- rbindlist(list(
  pca_all[mapping_method == "maximal_coordinate_overlap", .(
    panel = "B", source_name, unit_in_array, acre_id,
    mapping_method, mapping_overlap_bp, mapped_nov3_acre_id, h_cluster
  )],
  heat_localized[mapping_method == "maximal_coordinate_overlap", .(
    panel = "A", source_name, unit_in_array, acre_id,
    mapping_method, mapping_overlap_bp, mapped_nov3_acre_id, h_cluster
  )]
), use.names = TRUE)
fwrite(mapping_audit, file.path(data_dir, "H1_H3_overlap_mapping_audit.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "localized_arrays_panel_A",
    "transition_map_cells_panel_A",
    "edge_unavailable_arrays_excluded_panel_A",
    "all_HPRC_arrays_panel_B",
    "complete_units_panel_B",
    "ACRE_length_units_panel_C", "ACRE_length_arrays_panel_C",
    "transition_position_arrays_panel_C",
    "HOR_units_panel_C", "HOR_arrays_panel_C",
    "HOR_H1_mean_pct", "HOR_H2_mean_pct", "HOR_H3_mean_pct",
    "PCA_H1_units", "PCA_H2_units", "PCA_H3_units",
    "map_H1_cells", "map_H2_cells", "map_H3_cells",
    "ACRE_H1_median_length_kb", "ACRE_H2_median_length_kb",
    "ACRE_H3_median_length_kb",
    "PCA_overlap_resolved_units", "map_overlap_resolved_cells",
    "HOR_identity_threshold_pct"
  ),
  value = c(
    nrow(violin_arrays),
    nrow(heat_localized),
    1196L - nrow(violin_arrays),
    uniqueN(pca_all$source_name),
    nrow(pca_all),
    nrow(acre_length_plot),
    uniqueN(acre_length_plot$source_name),
    nrow(violin_arrays),
    nrow(hor_plot_all),
    uniqueN(hor_plot_all$source_name),
    hor_all_summary[h_cluster == "H1", mean_score_pct],
    hor_all_summary[h_cluster == "H2", mean_score_pct],
    hor_all_summary[h_cluster == "H3", mean_score_pct],
    pca_all_counts[h_cluster == "H1", units],
    pca_all_counts[h_cluster == "H2", units],
    pca_all_counts[h_cluster == "H3", units],
    heat_counts[h_cluster == "H1", units],
    heat_counts[h_cluster == "H2", units],
    heat_counts[h_cluster == "H3", units],
    acre_length_summary[h_cluster == "H1", median_length_kb],
    acre_length_summary[h_cluster == "H2", median_length_kb],
    acre_length_summary[h_cluster == "H3", median_length_kb],
    pca_all[mapping_method == "maximal_coordinate_overlap", .N],
    heat_localized[mapping_method == "maximal_coordinate_overlap", .N],
    90L
  )
)
fwrite(summary, file.path(data_dir, "Figure4_alt_summary.tsv"), sep = "\t")
fwrite(summary, file.path(data_dir, "Figure5_redesign_summary.tsv"), sep = "\t")

all_arrays_summary <- data.table(
  metric = c(
    "all_HPRC_arrays_panel_B",
    "all_HPRC_complete_units_panel_B",
    "all_PCA_H1_units", "all_PCA_H2_units", "all_PCA_H3_units",
    "all_PCA_overlap_resolved_units",
    "all_HOR_scored_units_panel_C",
    "all_HOR_H1_mean_pct", "all_HOR_H2_mean_pct", "all_HOR_H3_mean_pct",
    "curated_split_arrays_reference_set",
    "arrays_outside_curated_split_set",
    "outside_set_single_cluster_arrays",
    "outside_set_mixed_label_arrays",
    "outside_set_dominant_H1_arrays",
    "outside_set_dominant_H2_arrays",
    "outside_set_dominant_H3_arrays",
    "supplementary_midpoint_cells_H1",
    "supplementary_midpoint_cells_H2",
    "supplementary_midpoint_cells_H3"
  ),
  value = c(
    uniqueN(pca_all$source_name),
    nrow(pca_all),
    pca_all_counts[h_cluster == "H1", units],
    pca_all_counts[h_cluster == "H2", units],
    pca_all_counts[h_cluster == "H3", units],
    pca_all[mapping_method == "maximal_coordinate_overlap", .N],
    nrow(hor_plot_all),
    hor_all_summary[h_cluster == "H1", mean_score_pct],
    hor_all_summary[h_cluster == "H2", mean_score_pct],
    hor_all_summary[h_cluster == "H3", mean_score_pct],
    uniqueN(audit$source_name),
    nrow(non_split_array_summary),
    non_split_array_summary[unique_clusters == 1L, .N],
    non_split_array_summary[unique_clusters > 1L, .N],
    non_split_array_summary[dominant_cluster == "H1", .N],
    non_split_array_summary[dominant_cluster == "H2", .N],
    non_split_array_summary[dominant_cluster == "H3", .N],
    non_split_counts[h_cluster == "H1", units],
    non_split_counts[h_cluster == "H2", units],
    non_split_counts[h_cluster == "H3", units]
  )
)
fwrite(
  all_arrays_summary,
  file.path(all_arrays_data_dir, "Figure5_all_arrays_alt_summary.tsv"),
  sep = "\t"
)

cat("Created:", final_pdf, "\n")
cat("Created:", final_png, "\n")
cat("Created:", redesign_pdf, "\n")
cat("Created:", redesign_png, "\n")
cat("Created:", all_arrays_final_pdf, "\n")
cat("Created:", all_arrays_final_png, "\n")
cat("Created:", supp_non_split_pdf, "\n")
cat("Created:", supp_non_split_png, "\n")
print(summary)
print(all_arrays_summary)
