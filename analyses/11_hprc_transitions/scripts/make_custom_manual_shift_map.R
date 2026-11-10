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
source_root <- file.path(analysis_root, "full_1275_positive_breakpoint_map")
tables_dir <- file.path(out_root, "tables")
figures_dir <- file.path(out_root, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

index <- fread(file.path(source_root, "tables", "figure_row_index_1275.tsv"))
points <- fread(file.path(source_root, "tables", "positive_breakpoints_1275.tsv"))
clusters <- fread(file.path(
  project_root, "NTRIAGE_OCT26", "manuscript_output", "telomeric_ntriage_2794",
  "native", "telomeric_2794_core.clusters.tsv"
))

stopifnot(
  nrow(index) == 1275L,
  uniqueN(index$source_name) == 1275L,
  uniqueN(points$array_id) == 1275L,
  !anyDuplicated(clusters$acre_id)
)

# User-specified manual edits refer to the original overview row numbers.
index[, manual_action := "unchanged"]
index[plot_row_top_to_bottom %between% c(1L, 7L), manual_action := "anchor_shift_right_1_unit"]
index[plot_row_top_to_bottom %between% c(8L, 10L), manual_action := "removed"]
index[plot_row_top_to_bottom %between% c(633L, 1022L), manual_action := "anchor_shift_left_1_unit"]
index[, anchor_shift_units := fifelse(
  manual_action == "anchor_shift_right_1_unit", 1L,
  fifelse(manual_action == "anchor_shift_left_1_unit", -1L, 0L)
)]
index[, original_plot_row := plot_row_top_to_bottom]
kept <- index[manual_action != "removed"]
setorder(kept, original_plot_row)
kept[, adjusted_plot_row := seq_len(.N)]

metadata <- merge(
  clusters,
  points[, .(
    source_name = array_id, coordinate_source, breakpoint_local_bp,
    breakpoint_fraction_core, welch_score, within_minus_cross_ani
  )],
  by = "source_name", all = FALSE, sort = FALSE
)
metadata <- merge(
  metadata,
  kept[, .(
    source_name, original_plot_row, adjusted_plot_row, manual_action,
    anchor_shift_units, original_transition = transition
  )],
  by = "source_name", all = FALSE, sort = FALSE
)
setorder(metadata, adjusted_plot_row, unit_in_array)
metadata[, unit_midpoint := (start + end) / 2]

# Reconstruct the original boundary count used in the source figure, then move
# that boundary by the requested number of complete ACRO units.
metadata[, original_units_left := sum(unit_midpoint < breakpoint_local_bp), by = source_name]
metadata[, adjusted_units_left := original_units_left + anchor_shift_units]
if (metadata[, any(adjusted_units_left < 0L | adjusted_units_left > .N), by = source_name]$V1 |> any()) {
  stop("A manual shift moved an anchor beyond an array edge")
}
metadata[, units_from_adjusted_anchor := fifelse(
  unit_in_array <= adjusted_units_left,
  unit_in_array - adjusted_units_left - 1L,
  unit_in_array - adjusted_units_left
)]

nearest <- metadata[abs(units_from_adjusted_anchor) == 1L,
  .(source_name, units_from_adjusted_anchor, cluster)
]
nearest_wide <- dcast(
  nearest, source_name ~ units_from_adjusted_anchor, value.var = "cluster"
)
if (!"-1" %in% names(nearest_wide)) nearest_wide[, `-1` := NA_integer_]
if (!"1" %in% names(nearest_wide)) nearest_wide[, `1` := NA_integer_]
setnames(nearest_wide, c("-1", "1"), c("adjusted_left_cluster", "adjusted_right_cluster"))

audit <- unique(metadata[, .(
  source_name, original_plot_row, adjusted_plot_row, manual_action,
  anchor_shift_units, original_transition, coordinate_source,
  breakpoint_local_bp, original_units_left, adjusted_units_left,
  welch_score, within_minus_cross_ani
)])
audit <- merge(audit, nearest_wide, by = "source_name", all.x = TRUE)
audit[, adjusted_transition := fifelse(
  is.na(adjusted_left_cluster) | is.na(adjusted_right_cluster),
  "edge-unavailable",
  paste0("C", adjusted_left_cluster, "->C", adjusted_right_cluster)
)]
audit[, display_name := sub("_([^_]*)$", ":\\1", source_name)]
audit[, labeled_display_name := paste0(
  adjusted_plot_row, ". [original row ", original_plot_row, "] ", display_name
)]
setorder(audit, adjusted_plot_row)
fwrite(audit, file.path(tables_dir, "manual_shift_row_audit_1272.tsv"), sep = "\t")

heat <- metadata[abs(units_from_adjusted_anchor) <= 8L, .(
  source_name, acre_id, unit_in_array, start, end,
  units_from_adjusted_anchor, native_cluster = cluster,
  original_plot_row, adjusted_plot_row, manual_action, anchor_shift_units,
  coordinate_source, breakpoint_local_bp, welch_score, within_minus_cross_ani
)]
heat <- merge(
  heat,
  audit[, .(source_name, adjusted_transition, labeled_display_name)],
  by = "source_name", all.x = TRUE
)
setorder(heat, adjusted_plot_row, units_from_adjusted_anchor)
fwrite(heat, file.path(tables_dir, "manual_shift_units_plus_minus_8_1272.tsv"), sep = "\t")

pair_counts <- audit[, .(arrays = .N), by = adjusted_transition][order(-arrays)]
pair_counts[, fraction := arrays / sum(arrays)]
fwrite(pair_counts, file.path(tables_dir, "manual_shift_transition_counts_1272.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "source_rows", "rows_removed", "final_rows", "anchors_shifted_right_one_unit",
    "anchors_shifted_left_one_unit", "unchanged_rows", "adjusted_C1_to_C2",
    "adjusted_C6_to_C2", "adjusted_C2_to_C2", "edge_unavailable"
  ),
  value = c(
    nrow(index), sum(index$manual_action == "removed"), nrow(audit),
    sum(audit$manual_action == "anchor_shift_right_1_unit"),
    sum(audit$manual_action == "anchor_shift_left_1_unit"),
    sum(audit$manual_action == "unchanged"),
    sum(audit$adjusted_transition == "C1->C2"),
    sum(audit$adjusted_transition == "C6->C2"),
    sum(audit$adjusted_transition == "C2->C2"),
    sum(audit$adjusted_transition == "edge-unavailable")
  )
)
fwrite(summary, file.path(tables_dir, "manual_shift_summary.tsv"), sep = "\t")

native_palette <- c(
  "1" = "#4E79A7", "2" = "#F28E2B", "3" = "#E15759", "4" = "#76B7B2",
  "5" = "#59A14F", "6" = "#EDC948", "7" = "#B07AA1", "8" = "#FF9DA7"
)
observed_clusters <- sort(unique(as.character(heat$native_cluster)))

overview <- ggplot(
  heat,
  aes(units_from_adjusted_anchor, adjusted_plot_row, fill = factor(native_cluster))
) +
  geom_tile(width = 0.96, height = 0.96, linewidth = 0) +
  geom_vline(xintercept = 0, linewidth = 0.9, color = "black") +
  scale_fill_manual(
    values = native_palette, breaks = observed_clusters, limits = observed_clusters,
    drop = TRUE, name = "Native N-TRIAGE cluster"
  ) +
  scale_x_continuous(breaks = c(-8:-1, 1:8)) +
  scale_y_reverse(
    breaks = unique(c(1L, seq(100L, 1200L, 100L), nrow(audit))),
    expand = expansion(mult = c(0.002, 0.002))
  ) +
  labs(
    title = "Manually adjusted N-TRIAGE cluster map for HPRC transition-positive arrays",
    subtitle = "Original rows 1-7: anchor +1 unit; rows 8-10: removed; rows 633-1022: anchor -1 unit",
    x = "ACRO units from adjusted transition anchor",
    y = "Adjusted array row (see manual_shift_row_audit_1272.tsv)",
    caption = paste0(
      "The original row order is retained after removing rows 8-10. Manual shifts move the transition anchor, not cluster identities.\n",
      "White cells indicate unavailable units at array edges; the black line marks the adjusted anchor."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 8.5), legend.position = "right"
  )

overview_png <- file.path(figures_dir, "cluster_map_plus_minus_8_manual_shifts_1272.png")
overview_pdf <- file.path(figures_dir, "cluster_map_plus_minus_8_manual_shifts_1272.pdf")
ggsave(overview_png, overview, width = 12, height = 16, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave(overview_pdf, overview, width = 12, height = 16, bg = "white", device = grDevices::pdf, useDingbats = FALSE)

page_size <- 50L
pages <- split(seq_len(nrow(audit)), ceiling(seq_len(nrow(audit)) / page_size))
labeled_pdf <- file.path(figures_dir, "cluster_map_plus_minus_8_manual_shifts_1272_labeled_pages.pdf")
grDevices::pdf(labeled_pdf, width = 12, height = 10, onefile = TRUE, useDingbats = FALSE)
for (page_number in seq_along(pages)) {
  rows <- pages[[page_number]]
  page_index <- audit[adjusted_plot_row %in% rows]
  page_heat <- heat[adjusted_plot_row %in% rows]
  ordered_labels <- page_index[order(adjusted_plot_row), labeled_display_name]
  page_heat[, page_label := factor(labeled_display_name, levels = rev(ordered_labels))]
  page_plot <- ggplot(
    page_heat,
    aes(units_from_adjusted_anchor, page_label, fill = factor(native_cluster))
  ) +
    geom_tile(color = "white", linewidth = 0.12) +
    geom_vline(xintercept = 0, linewidth = 0.8, color = "black") +
    scale_fill_manual(
      values = native_palette, breaks = observed_clusters, limits = observed_clusters,
      drop = TRUE, name = "Native N-TRIAGE cluster"
    ) +
    scale_x_continuous(breaks = c(-8:-1, 1:8)) +
    labs(
      title = "Manually adjusted N-TRIAGE clusters around HPRC transitions",
      subtitle = paste0(
        "Adjusted rows ", min(rows), "-", max(rows), " of ", nrow(audit),
        " (page ", page_number, " of ", length(pages), ")"
      ),
      x = "ACRO units from adjusted transition anchor", y = NULL
    ) +
    theme_bw(base_size = 9.5) +
    theme(
      panel.grid = element_blank(), axis.text.y = element_text(size = 5.0),
      plot.title = element_text(face = "bold"), legend.position = "right"
    )
  print(page_plot)
}
grDevices::dev.off()

cat("Wrote manually adjusted 1,272-row overview and labeled companion\n")
cat("PNG:", overview_png, "\n")
cat("PDF:", overview_pdf, "\n")
cat("Labeled PDF:", labeled_pdf, "\n")
