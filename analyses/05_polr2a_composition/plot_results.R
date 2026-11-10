#!/usr/bin/env Rscript

# Plot the codon-degeneracy-normalized ACRO1--POLR2A CTD JSD null tests.
# This script intentionally uses only base R so the analysis is portable.

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1L) args[[1L]] else "."
summary_path <- file.path(analysis_dir, "results", "jsd_null_summary.tsv")
histogram_path <- file.path(analysis_dir, "results", "jsd_null_histogram.tsv")
figure_dir <- file.path(analysis_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_data <- read.delim(
  summary_path, sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = ""
)
histogram_data <- read.delim(
  histogram_path, sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = ""
)

models <- c(
  "human_genome_iid",
  "exact_mononucleotide_shuffle",
  "dinucleotide_shuffle"
)
model_labels <- c(
  "Human-genome nucleotide background",
  "Exact mononucleotide shuffle",
  "Dinucleotide-preserving shuffle"
)
frames <- c("candidate_frame0", "all_three_forward_frames")
frame_labels <- c(
  "Candidate ACRO1 reading frame",
  "All three forward frames (RAP-style)"
)
ctd_reference <- "full_CTD_1593_1960"
normalization <- "codon_degeneracy_normalized"

selected_summary <- summary_data[
  summary_data$null_model %in% models & summary_data$frame_scheme %in% frames &
    summary_data$ctd_reference == ctd_reference &
    summary_data$normalization == normalization, , drop = FALSE
]
selected_histogram <- histogram_data[
  histogram_data$null_model %in% models & histogram_data$frame_scheme %in% frames &
    histogram_data$ctd_reference == ctd_reference &
    histogram_data$normalization == normalization, , drop = FALSE
]

if (nrow(selected_summary) != 6L) {
  stop("Expected six summary rows, found ", nrow(selected_summary))
}
if (nrow(selected_histogram) != 6L * 200L) {
  stop("Expected 1,200 histogram rows, found ", nrow(selected_histogram))
}

nonzero <- selected_histogram[selected_histogram$count > 0, , drop = FALSE]
x_values <- c(
  nonzero$bin_start, nonzero$bin_end,
  selected_summary$observed_jsd_base2,
  selected_summary$null_mean
)
x_range <- range(x_values, finite = TRUE)
x_pad <- 0.025 * diff(x_range)
x_limits <- c(max(0, x_range[[1L]] - x_pad), min(1, x_range[[2L]] + x_pad))

format_empirical_p <- function(row) {
  hits <- as.numeric(row$null_jsd_less_equal_observed)
  reps <- as.numeric(row$null_replicates)
  p <- as.numeric(row$empirical_p_plus1)
  if (hits == 0) {
    paste0("Pemp = ", format(p, scientific = TRUE, digits = 2),
           "  (0/", format(reps, big.mark = ",", scientific = FALSE), ")")
  } else if (p < 0.001) {
    paste0("Pemp = ", format(p, scientific = TRUE, digits = 2),
           "  (", format(hits, big.mark = ",", scientific = FALSE), "/",
           format(reps, big.mark = ",", scientific = FALSE), ")")
  } else {
    paste0("Pemp = ", formatC(p, format = "f", digits = 4),
           "  (", format(hits, big.mark = ",", scientific = FALSE), "/",
           format(reps, big.mark = ",", scientific = FALSE), ")")
  }
}

draw_figure <- function(device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(
      file.path(figure_dir, "ACRO1_POLR2A_JSD_nulls.png"),
      width = 3600, height = 2250, res = 300, bg = "white"
    )
  } else {
    pdf(
      file.path(figure_dir, "ACRO1_POLR2A_JSD_nulls.pdf"),
      width = 12, height = 7.5, family = "Helvetica", useDingbats = FALSE
    )
  }
  on.exit(dev.off(), add = TRUE)

  layout(matrix(seq_len(6L), nrow = 2L, byrow = TRUE))
  par(
    oma = c(4.3, 5.2, 1.1, 0.7),
    mar = c(2.2, 2.7, 2.9, 0.8),
    mgp = c(1.75, 0.45, 0), tcl = -0.22,
    las = 1, family = "sans", xaxs = "i", yaxs = "i"
  )

  panel_letter <- 0L
  for (frame_index in seq_along(frames)) {
    for (model_index in seq_along(models)) {
      panel_letter <- panel_letter + 1L
      frame <- frames[[frame_index]]
      model <- models[[model_index]]
      hist_row <- selected_histogram[
        selected_histogram$frame_scheme == frame & selected_histogram$null_model == model,
        , drop = FALSE
      ]
      summary_row <- selected_summary[
        selected_summary$frame_scheme == frame & selected_summary$null_model == model,
        , drop = FALSE
      ]
      hist_row <- hist_row[order(hist_row$bin_start), , drop = FALSE]
      bin_width <- hist_row$bin_end - hist_row$bin_start
      density <- hist_row$count / (as.numeric(summary_row$null_replicates) * bin_width)
      y_limit <- max(density) * 1.22

      plot(
        NA, NA, type = "n", xlim = x_limits, ylim = c(0, y_limit),
        axes = FALSE, xlab = "", ylab = ""
      )
      abline(h = pretty(c(0, y_limit), n = 4L), col = "#ECEFF2", lwd = 0.8)
      rect(
        hist_row$bin_start, 0, hist_row$bin_end, density,
        col = "#AFC7E3", border = "#5B84B1", lwd = 0.35
      )
      lines(
        (hist_row$bin_start + hist_row$bin_end) / 2, density,
        col = "#315F8C", lwd = 1.0
      )
      abline(v = as.numeric(summary_row$null_mean), col = "#4D4D4D", lty = 3, lwd = 1.3)
      abline(v = as.numeric(summary_row$observed_jsd_base2), col = "#C83E4D", lwd = 2.2)

      axis(1, at = pretty(x_limits, n = 5L), cex.axis = 0.78)
      if (model_index == 1L) {
        axis(2, at = pretty(c(0, y_limit), n = 4L), cex.axis = 0.75)
        mtext("Null density", side = 2, line = 1.75, cex = 0.78, las = 0)
        mtext(frame_labels[[frame_index]], side = 2, line = 4.2, cex = 0.87, font = 2, las = 0)
      }
      box(col = "#666666", lwd = 0.7)

      if (frame_index == 1L) {
        title(main = model_labels[[model_index]], cex.main = 0.91, font.main = 2, line = 1.05)
      }
      mtext(LETTERS[[panel_letter]], side = 3, adj = 0, line = 1.15, cex = 0.96, font = 2)
      p_label <- format_empirical_p(summary_row)
      text(
        x_limits[[2L]] - 0.018 * diff(x_limits), y_limit * 0.93,
        labels = p_label, adj = c(1, 1), cex = 0.72, col = "#222222"
      )
      text(
        x_limits[[2L]] - 0.018 * diff(x_limits), y_limit * 0.81,
        labels = paste0("Observed JSD = ", formatC(as.numeric(summary_row$observed_jsd_base2), digits = 3, format = "f")),
        adj = c(1, 1), cex = 0.69, col = "#C83E4D"
      )

      if (frame_index == 1L && model_index == 1L) {
        legend(
          "topleft", inset = c(0.025, 0.16),
          legend = c("Null distribution", "Null mean", "Observed"),
          fill = c("#AFC7E3", NA, NA), border = c("#5B84B1", NA, NA),
          lty = c(NA, 3, 1), lwd = c(NA, 1.3, 2.2),
          col = c("#315F8C", "#4D4D4D", "#C83E4D"),
          bty = "n", cex = 0.67, y.intersp = 0.9, x.intersp = 0.7
        )
      }
    }
  }
  mtext(
    "Jensen-Shannon divergence from the full human POLR2A CTD (lower = more similar)",
    side = 1, outer = TRUE, line = 2.15, cex = 0.94
  )
  mtext(
    "Codon-degeneracy-normalized amino-acid composition; 1,000,000 simulations per panel",
    side = 1, outer = TRUE, line = 3.25, cex = 0.76, col = "#444444"
  )
}

draw_figure("png")
draw_figure("pdf")

message("Wrote ", file.path(figure_dir, "ACRO1_POLR2A_JSD_nulls.png"))
message("Wrote ", file.path(figure_dir, "ACRO1_POLR2A_JSD_nulls.pdf"))
