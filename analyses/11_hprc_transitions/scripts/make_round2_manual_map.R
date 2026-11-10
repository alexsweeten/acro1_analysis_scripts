#!/usr/bin/env Rscript

options(warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
out_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(out_root, ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)
round1_root <- file.path(analysis_root, "custom_manual_shift_breakpoint_map")
tables_dir <- file.path(out_root, "tables")
figures_dir <- file.path(out_root, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

round1 <- fread(file.path(round1_root, "tables", "manual_shift_row_audit_1272.tsv"))
clusters <- fread(file.path(
  project_root, "NTRIAGE_OCT26", "manuscript_output", "telomeric_ntriage_2794",
  "native", "telomeric_2794_core.clusters.tsv"
))
stopifnot(nrow(round1) == 1272L, uniqueN(round1$source_name) == 1272L)

round1[, round1_row := adjusted_plot_row]
round1[, round2_action := "unchanged"]
remove_rows <- c(
  616L, 619L:626L, 1057L:1063L, 1074L:1098L, 1238L:1272L
)
round1[round1_row %in% remove_rows, round2_action := "removed"]
round1[round1_row == 618L, round2_action := "anchor_red_to_blue"]
round1[round1_row %between% c(627L, 636L), round2_action := "anchor_after_last_orange"]
round1[round1_row == 745L, round2_action := "anchor_after_last_orange"]
round1[round1_row %between% c(1099L, 1103L), round2_action := "all_ACRO_units_left"]

kept <- round1[round2_action != "removed"]
setorder(kept, round1_row)
kept[, round2_row := seq_len(.N)]

metadata <- merge(
  clusters,
  kept[, .(
    source_name, round1_row, round2_row, original_plot_row,
    round1_manual_action = manual_action,
    round1_transition = adjusted_transition,
    round1_units_left = adjusted_units_left,
    round2_action, coordinate_source, breakpoint_local_bp,
    welch_score, within_minus_cross_ani
  )],
  by = "source_name", all = FALSE, sort = FALSE
)
setorder(metadata, round2_row, unit_in_array)
metadata[, round2_units_left := round1_units_left]

# Row 618: use its unique C3-to-C1 boundary (red on left, blue on right).
red_blue <- metadata[round1_row == 618L]
red_blue_boundaries <- red_blue[
  shift(cluster, type = "lead") == 1L & cluster == 3L,
  unit_in_array
]
if (length(red_blue_boundaries) != 1L) stop("Row 618 does not have one unique C3-to-C1 boundary")
metadata[round1_row == 618L, round2_units_left := red_blue_boundaries[[1]]]

# Rows 627-636 and 745: put every C2 (orange) unit to the left of the split.
orange_rows <- c(627L:636L, 745L)
last_orange <- metadata[
  round1_row %in% orange_rows & cluster == 2L,
  .(target_units_left = max(unit_in_array)), by = round1_row
]
if (nrow(last_orange) != length(orange_rows)) stop("At least one orange-adjustment row lacks C2")
metadata <- merge(metadata, last_orange, by = "round1_row", all.x = TRUE, sort = FALSE)
metadata[round1_row %in% orange_rows, round2_units_left := target_units_left]
metadata[, target_units_left := NULL]

# Rows 1099-1103: place the transition anchor after the final complete unit.
all_left_rows <- 1099L:1103L
metadata[round1_row %in% all_left_rows, round2_units_left := .N, by = source_name]

if (metadata[, uniqueN(round2_units_left), by = source_name][V1 != 1L, .N] != 0L) {
  stop("A source array received inconsistent round-2 anchors")
}
if (metadata[, any(round2_units_left < 0L | round2_units_left > .N), by = source_name]$V1 |> any()) {
  stop("A round-2 anchor lies beyond an array edge")
}

metadata[, units_from_round2_anchor := fifelse(
  unit_in_array <= round2_units_left,
  unit_in_array - round2_units_left - 1L,
  unit_in_array - round2_units_left
)]

nearest <- metadata[abs(units_from_round2_anchor) == 1L,
  .(source_name, units_from_round2_anchor, cluster)
]
nearest_wide <- dcast(nearest, source_name ~ units_from_round2_anchor, value.var = "cluster")
if (!"-1" %in% names(nearest_wide)) nearest_wide[, `-1` := NA_integer_]
if (!"1" %in% names(nearest_wide)) nearest_wide[, `1` := NA_integer_]
setnames(nearest_wide, c("-1", "1"), c("round2_left_cluster", "round2_right_cluster"))

audit <- unique(metadata[, .(
  source_name, round1_row, round2_row, original_plot_row,
  round1_manual_action, round2_action, round1_transition,
  round1_units_left, round2_units_left, coordinate_source,
  breakpoint_local_bp, welch_score, within_minus_cross_ani
)])
audit <- merge(audit, nearest_wide, by = "source_name", all.x = TRUE)
audit[, round2_transition := fifelse(
  is.na(round2_left_cluster) | is.na(round2_right_cluster),
  "edge-unavailable",
  paste0("C", round2_left_cluster, "->C", round2_right_cluster)
)]
audit[, anchor_change_from_round1_units := round2_units_left - round1_units_left]
audit[, display_name := sub("_([^_]*)$", ":\\1", source_name)]
audit[, labeled_display_name := paste0(
  round2_row, ". [round-1 row ", round1_row, "] ", display_name
)]
setorder(audit, round2_row)
fwrite(audit, file.path(tables_dir, "round2_row_audit_1196.tsv"), sep = "\t")

removed <- round1[round2_action == "removed", .(
  source_name, removed_round1_row = round1_row, original_plot_row,
  round1_transition = adjusted_transition, removal_action = round2_action
)]
setorder(removed, removed_round1_row)
fwrite(removed, file.path(tables_dir, "round2_removed_arrays_76.tsv"), sep = "\t")

heat <- metadata[abs(units_from_round2_anchor) <= 8L, .(
  source_name, acre_id, unit_in_array, start, end,
  units_from_round2_anchor, native_cluster = cluster,
  round1_row, round2_row, original_plot_row, round2_action,
  round1_units_left, round2_units_left, coordinate_source,
  breakpoint_local_bp, welch_score, within_minus_cross_ani
)]
heat <- merge(
  heat,
  audit[, .(source_name, round2_transition, labeled_display_name)],
  by = "source_name", all.x = TRUE
)
setorder(heat, round2_row, units_from_round2_anchor)
fwrite(heat, file.path(tables_dir, "round2_units_plus_minus_8_1196.tsv"), sep = "\t")

pair_counts <- audit[, .(arrays = .N), by = round2_transition][order(-arrays)]
pair_counts[, fraction := arrays / sum(arrays)]
fwrite(pair_counts, file.path(tables_dir, "round2_transition_counts_1196.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "round1_rows", "removed_rows", "round2_rows", "adjusted_rows_round2",
    "unchanged_rows_round2", "round2_C1_to_C2", "round2_C6_to_C2",
    "round2_C2_to_C2", "round2_edge_unavailable"
  ),
  value = c(
    nrow(round1), nrow(removed), nrow(audit),
    sum(audit$round2_action != "unchanged"),
    sum(audit$round2_action == "unchanged"),
    sum(audit$round2_transition == "C1->C2"),
    sum(audit$round2_transition == "C6->C2"),
    sum(audit$round2_transition == "C2->C2"),
    sum(audit$round2_transition == "edge-unavailable")
  )
)
fwrite(summary, file.path(tables_dir, "round2_summary.tsv"), sep = "\t")

native_palette <- c(
  "1" = "#4E79A7", "2" = "#F28E2B", "3" = "#E15759", "4" = "#76B7B2",
  "5" = "#59A14F", "6" = "#EDC948", "7" = "#B07AA1", "8" = "#FF9DA7"
)
observed_clusters <- sort(unique(as.character(heat$native_cluster)))

overview <- ggplot(
  heat, aes(units_from_round2_anchor, round2_row, fill = factor(native_cluster))
) +
  geom_tile(width = 0.96, height = 0.96, linewidth = 0) +
  geom_vline(xintercept = 0, linewidth = 0.9, color = "black") +
  scale_fill_manual(
    values = native_palette, breaks = observed_clusters, limits = observed_clusters,
    drop = TRUE, name = "Native N-TRIAGE cluster"
  ) +
  scale_x_continuous(breaks = c(-8:-1, 1:8)) +
  scale_y_reverse(
    breaks = unique(c(1L, seq(100L, 1100L, 100L), nrow(audit))),
    expand = expansion(mult = c(0.002, 0.002))
  ) +
  labs(
    title = "Second manually curated N-TRIAGE cluster map for HPRC arrays",
    subtitle = "1,196 retained arrays; row-specific removals and anchor adjustments are recorded in the audit table",
    x = "ACRO units from round-2 transition anchor",
    y = "Round-2 array row (see round2_row_audit_1196.tsv)",
    caption = paste0(
      "Row instructions refer to the prior 1,272-row manually adjusted figure. Native cluster calls are unchanged.\n",
      "White cells indicate unavailable units at array edges; the black line marks the round-2 anchor."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 8.5), legend.position = "right"
  )

overview_png <- file.path(figures_dir, "cluster_map_plus_minus_8_manual_round2_1196.png")
overview_pdf <- file.path(figures_dir, "cluster_map_plus_minus_8_manual_round2_1196.pdf")
ggsave(overview_png, overview, width = 12, height = 16, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(overview_pdf, overview, width = 12, height = 16, bg = "white", device = grDevices::pdf, useDingbats = FALSE)

page_size <- 50L
pages <- split(seq_len(nrow(audit)), ceiling(seq_len(nrow(audit)) / page_size))
labeled_pdf <- file.path(figures_dir, "cluster_map_plus_minus_8_manual_round2_1196_labeled_pages.pdf")
grDevices::pdf(labeled_pdf, width = 12, height = 10, onefile = TRUE, useDingbats = FALSE)
for (page_number in seq_along(pages)) {
  rows <- pages[[page_number]]
  page_index <- audit[round2_row %in% rows]
  page_heat <- heat[round2_row %in% rows]
  ordered_labels <- page_index[order(round2_row), labeled_display_name]
  page_heat[, page_label := factor(labeled_display_name, levels = rev(ordered_labels))]
  page_plot <- ggplot(
    page_heat,
    aes(units_from_round2_anchor, page_label, fill = factor(native_cluster))
  ) +
    geom_tile(color = "white", linewidth = 0.12) +
    geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
    scale_fill_manual(
      values = native_palette, breaks = observed_clusters, limits = observed_clusters,
      drop = TRUE, name = "Native N-TRIAGE cluster"
    ) +
    scale_x_continuous(breaks = c(-8:-1, 1:8)) +
    labs(
      title = "Second manually curated N-TRIAGE clusters around HPRC transitions",
      subtitle = paste0(
        "Round-2 rows ", min(rows), "-", max(rows), " of ", nrow(audit),
        " (page ", page_number, " of ", length(pages), ")"
      ),
      x = "ACRO units from round-2 transition anchor", y = NULL
    ) +
    theme_bw(base_size = 9.5) +
    theme(
      panel.grid = element_blank(), axis.text.y = element_text(size = 5.0),
      plot.title = element_text(face = "bold"), legend.position = "right"
    )
  print(page_plot)
}
grDevices::dev.off()

cat("Wrote round-2 manually curated 1,196-row overview and labeled companion\n")
cat("PNG:", overview_png, "\n")
cat("PDF:", overview_pdf, "\n")
cat("Labeled PDF:", labeled_pdf, "\n")
