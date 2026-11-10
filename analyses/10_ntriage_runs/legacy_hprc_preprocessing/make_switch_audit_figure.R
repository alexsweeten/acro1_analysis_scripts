#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
tables <- file.path(root, "tables")
figures <- file.path(root, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

status_palette <- c(
  `No annotated\nrecombination` = "#6A6A6A",
  `Recombination-\nannotated` = "#D55E00"
)
acro1_palette <- c(
  `Total-variation\ndistance` = "#4477AA",
  `|change in S1|\ncanonical` = "#009E73",
  `|change in S2|\ninternal CCC deletion` = "#E69F00",
  `|change in S3|\nterminal deletion` = "#CC79A7"
)

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(color = "#444444", size = rel(0.84), lineheight = 1.05),
      plot.caption = element_text(color = "#555555", size = rel(0.72), hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#222222"),
      legend.position = "none"
    )
}

save_both <- function(plot, stem, width, height) {
  ggsave(file.path(figures, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 400, bg = "white")
  ggsave(file.path(figures, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = "pdf", useDingbats = FALSE, bg = "white")
}

# Panel A: raw pooled rates, with the donor-matched contrast kept visually distinct.
pooled <- fread(file.path(tables, "core_switch_acro1_pooled_counts.tsv"))
matched <- fread(file.path(tables, "core_switch_acro1_donor_matched_metrics.tsv"))
rate_effect <- matched[metric == "switch_rate100"]

rate_plot <- pooled[scope %in% c("no_recombination", "recombination"), .(
  scope,
  collapsed_switch_events,
  evaluable_boundaries,
  switch_events_per_100_evaluable_boundaries
)]
rate_plot[, status := factor(
  fifelse(scope == "no_recombination", "No annotated\nrecombination", "Recombination-\nannotated"),
  levels = names(status_palette)
)]
rate_plot[, count_label := sprintf("%s / %s\nboundaries",
                                   comma(collapsed_switch_events), comma(evaluable_boundaries))]
rate_ymax <- 5.15

p_rate <- ggplot(rate_plot, aes(status, switch_events_per_100_evaluable_boundaries, fill = status)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.35) +
  geom_text(aes(label = count_label), vjust = -0.35, size = 3.15, lineheight = 0.95) +
  annotate("segment", x = 1, xend = 2, y = 4.42, yend = 4.42, linewidth = 0.45, color = "#333333") +
  annotate("segment", x = 1, xend = 1, y = 4.30, yend = 4.42, linewidth = 0.45, color = "#333333") +
  annotate("segment", x = 2, xend = 2, y = 4.30, yend = 4.42, linewidth = 0.45, color = "#333333") +
  annotate(
    "text", x = 1.5, y = 4.63,
    label = sprintf("Donor-matched difference = %.2f\n95%% CI %.2f to %.2f",
                    rate_effect$mean_within_donor_recombination_minus_no_recombination,
                    rate_effect$bootstrap_ci95_lo, rate_effect$bootstrap_ci95_hi),
    size = 3.05, lineheight = 0.95
  ) +
  scale_fill_manual(values = status_palette) +
  scale_y_continuous(limits = c(0, rate_ymax), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Pooled switch rates",
    subtitle = "Raw collapsed switches per 100 evaluable boundaries",
    x = NULL, y = "Switches per 100 boundaries"
  ) +
  theme_manuscript() +
  theme(axis.text.x = element_text(size = 9.3, lineheight = 0.9))

# Panel B: local ACRO1-composition contrasts around switch events.
local <- fread(file.path(tables, "core_switch_local_acro1_tests.tsv"))
local <- local[test %chin% c("donor_matched_TVD", "donor_matched_abs_dS1",
                             "donor_matched_abs_dS2", "donor_matched_abs_dS3")]
metric_labels <- c(
  TVD = "Total-variation\ndistance",
  abs_dS1 = "|change in S1|\ncanonical",
  abs_dS2 = "|change in S2|\ninternal CCC deletion",
  abs_dS3 = "|change in S3|\nterminal deletion"
)
local[, metric_label := factor(metric_labels[metric], levels = rev(unname(metric_labels)))]
local[, `:=`(
  effect_pct = 100 * effect,
  lo_pct = 100 * bootstrap_ci95_lo,
  hi_pct = 100 * bootstrap_ci95_hi,
  effect_label = sprintf("%+.2f", 100 * effect)
)]

p_local <- ggplot(local, aes(effect_pct, metric_label, color = metric_label)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lo_pct, xmax = hi_pct), orientation = "y", width = 0.17,
                linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_text(aes(x = hi_pct + 0.48, label = effect_label), color = "#222222", hjust = 0, size = 3.15) +
  scale_color_manual(values = acro1_palette) +
  scale_x_continuous(limits = c(-3.0, 11.7), breaks = seq(-2.5, 10, by = 2.5)) +
  labs(
    title = "Local ACRO1 changes at switches",
    subtitle = "Recombination-switch minus no-recombination-switch mean; 74 donors",
    x = "Mean within-donor difference (percentage points)", y = NULL
  ) +
  theme_manuscript() +
  theme(axis.text.y = element_text(size = 9.2, lineheight = 0.92))

# Panel C: descriptive dominant-class changes with unambiguous raw denominators.
dominant <- fread(file.path(tables, "core_switch_local_dominant_allele_changes.tsv"))
dominant[, status := factor(
  fifelse(has_recombination == 0, "No annotated\nrecombination", "Recombination-\nannotated"),
  levels = names(status_palette)
)]
dominant[, count_label := sprintf("%s / %s\n(%.1f%%)",
                                  comma(dominant_ACRO1_allele_change_events),
                                  comma(collapsed_switch_events),
                                  100 * dominant_ACRO1_allele_change_fraction)]

p_dominant <- ggplot(dominant, aes(status, dominant_ACRO1_allele_change_fraction, fill = status)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.35) +
  geom_text(aes(y = pmax(dominant_ACRO1_allele_change_fraction, 0.008), label = count_label),
            vjust = -0.35, size = 3.35, lineheight = 0.95) +
  scale_fill_manual(values = status_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.31),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "Dominant ACRO1 class changes",
    subtitle = "Pooled over three units on each side of a collapsed switch",
    x = NULL, y = "Switch events with a class change"
  ) +
  theme_manuscript() +
  theme(axis.text.x = element_text(size = 9.3, lineheight = 0.9))

switch_plot <- p_rate + p_local + p_dominant +
  plot_layout(widths = c(1.02, 1.32, 1.02)) +
  plot_annotation(
    title = "Deep-subfamily switches and local ACRO1 composition in telomeric ACRO arrays",
    subtitle = "Final core-only N-TRIAGE audit: 1,619 collapsed switches among 63,306 evaluable boundaries",
    caption = paste0(
      "Raw pooled rates and event fractions are descriptive. Confidence intervals are 20,000-draw donor-bootstrap intervals from donor-only matching ",
      "(not donor x haplotype-phase matching); the local ACRO1 comparison includes donors with switch events under both labels.\n",
      "Panel C is unmodeled. Recombination status is an array-level annotation, so these associations do not establish causation."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, color = "#444444"),
      plot.caption = element_text(size = 8, color = "#555555", hjust = 0, lineheight = 1.05),
      plot.margin = margin(8, 10, 8, 10)
    )
  )

save_both(switch_plot, "Supplementary_Figure_telomeric_NTRIAGE_switch_audit", 13.4, 6.15)
cat(sprintf("Wrote switch/ACRO1 audit figure to %s\n", figures))
