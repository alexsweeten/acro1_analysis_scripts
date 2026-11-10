args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
default_root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."), mustWork = FALSE)
} else {
  normalizePath(".", mustWork = FALSE)
}
root <- if (length(args)) args[[1]] else default_root
d <- read.delim(file.path(root, "tables", "sampled_200_units.tsv"), check.names = FALSE)

metrics <- c("unit_length_l1pa10", "acro1_union_bp", "non_acro1_remainder_bp", "explicit_backbone_bp")
labels <- c(
  "ACRO unit length (bp)\nL1PA10-to-L1PA10",
  "ACRO1 anchor coverage per unit (bp)",
  "Sequence outside exact ACRO1 anchors (bp)",
  "Explicitly annotated backbone (bp)"
)

draw_plot <- function(device, path) {
  if (device == "png") png(path, width = 4800, height = 1300, res = 300)
  if (device == "pdf") pdf(path, width = 16, height = 4.3, useDingbats = FALSE)
  par(mfrow = c(1, 4), mar = c(4.5, 4.6, 3, 1), oma = c(0, 0, 2, 0))
  for (i in seq_along(metrics)) {
    m <- metrics[[i]]
    w <- reshape(d[, c("source_name", "side", m)], idvar = "source_name", timevar = "side", direction = "wide")
    left <- w[[paste0(m, ".left")]]
    right <- w[[paste0(m, ".right")]]
    ylim <- range(c(left, right), na.rm = TRUE)
    plot(c(1, 2), ylim, type = "n", xaxt = "n", xlab = "", ylab = labels[[i]], bty = "l")
    for (j in seq_along(left)) lines(c(1, 2), c(left[j], right[j]), col = adjustcolor("grey40", alpha.f = 0.18), lwd = 0.7)
    boxplot(left, right, at = c(1, 2), add = TRUE, axes = FALSE, outline = FALSE,
            col = c("#D9D9D9", "#F28E2B"), border = "black", boxwex = 0.42)
    axis(1, at = c(1, 2), labels = c("Left\n(non-C2)", "Right\n(C2)"))
    title(main = paste0("Mean R-L = ", round(mean(right - left, na.rm = TRUE), 1), " bp"), cex.main = 0.95)
  }
  mtext("Paired 200-unit HPRC split subsample", outer = TRUE, cex = 1.25, font = 2)
  dev.off()
}

draw_plot("png", file.path(root, "figures", "paired_200_unit_length_acro1_comparison.png"))
draw_plot("pdf", file.path(root, "figures", "paired_200_unit_length_acro1_comparison.pdf"))
