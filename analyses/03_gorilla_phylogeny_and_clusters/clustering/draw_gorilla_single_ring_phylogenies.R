#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))

tree <- read.tree(file.path(out_dir, "gorilla_cluster_hsa_nor_display_rooted.treefile"))
meta <- read.delim(
  file.path(out_dir, "gorilla_cluster_hsa_nor_phylogeny_annotations.tsv"),
  stringsAsFactors = FALSE
)
palette_table <- read.delim(
  file.path(out_dir, "gorilla_cluster_hsa_nor_phylogeny_palette.tsv"),
  stringsAsFactors = FALSE
)

stopifnot(!anyDuplicated(meta$tip_id), setequal(tree$tip.label, meta$tip_id))
meta <- meta[match(tree$tip.label, meta$tip_id), , drop = FALSE]

ring_geometry <- function() {
  state <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x <- state$xx[seq_len(state$Ntip)]
  y <- state$yy[seq_len(state$Ntip)]
  list(angles = atan2(y, x), tip_radius = median(sqrt(x^2 + y^2)))
}

draw_annular_ring <- function(colors) {
  geometry <- ring_geometry()
  angles <- geometry$angles
  inner <- geometry$tip_radius * 1.022
  outer <- geometry$tip_radius * 1.088
  half_width <- pi / length(angles) * 1.025
  for (index in seq_along(angles)) {
    theta <- seq(angles[index] - half_width, angles[index] + half_width, length.out = 5)
    polygon(
      c(inner * cos(theta), rev(outer * cos(theta))),
      c(inner * sin(theta), rev(outer * sin(theta))),
      col = colors[index], border = NA
    )
  }
  symbols(0, 0, circles = inner, inches = FALSE, add = TRUE, fg = "white", lwd = 0.7)
  symbols(0, 0, circles = outer, inches = FALSE, add = TRUE, fg = "white", lwd = 0.7)
}

draw_legend_panel <- function(labels, colors, counts, legend_title) {
  par(mar = c(0.5, 0.2, 5.0, 0.5), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend(
    x = 0.02, y = 0.88,
    legend = paste0(labels, "  (n=", counts, ")"),
    pch = 22, pt.bg = colors, col = "#181818", pt.lwd = 0.8,
    pt.cex = if (length(labels) > 8) 1.55 else 1.75,
    bty = "n", cex = if (length(labels) > 8) 0.88 else 0.96,
    y.intersp = if (length(labels) > 8) 1.10 else 1.25,
    title = legend_title
  )
  text(
    0.02, 0.09,
    paste0(
      Ntip(tree), " gorilla ACREs\n",
      "Black internal branches\n",
      "Branch lengths suppressed for legibility"
    ),
    adj = c(0, 0), cex = 0.75, col = "#5F696E"
  )
}

draw_single_ring_tree <- function(ring_key, metadata_column, title_label, legend_title) {
  key <- palette_table[palette_table$ring == ring_key, , drop = FALSE]
  labels <- key$label
  colors <- setNames(key$color, labels)
  values <- meta[[metadata_column]]
  tip_colors <- unname(colors[values])
  counts <- as.integer(table(factor(values, levels = labels)))
  stopifnot(!anyNA(tip_colors), sum(counts) == Ntip(tree))

  layout(matrix(c(1, 2), nrow = 1), widths = c(5.3, 1.35))
  par(mar = c(0.4, 0.4, 5.6, 0.2), family = "sans", xpd = NA, fg = "#111111")
  terminal_edge <- tree$edge[, 2] <= Ntip(tree)
  plot(
    tree,
    type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
    rotate.tree = 4,
    edge.color = ifelse(terminal_edge, "#333333", "#050505"),
    edge.width = ifelse(terminal_edge, 0.30, 0.82),
    no.margin = FALSE, x.lim = c(-1.24, 1.24), y.lim = c(-1.24, 1.24)
  )
  draw_annular_ring(tip_colors)
  title(
    main = paste0("Gorilla ACRE phylogeny: ", title_label),
    cex.main = 1.45, col.main = "#111111", line = 3.4, font.main = 2
  )
  mtext(
    paste0(Ntip(tree), " ACREs  |  single exterior ring: ", title_label),
    side = 3, line = 1.15, cex = 0.84, col = "#56636A"
  )
  draw_legend_panel(labels, unname(colors[labels]), counts, legend_title)
}

render <- function(stem, ring_key, metadata_column, title_label, legend_title) {
  ragg::agg_png(
    file.path(out_dir, paste0(stem, ".png")),
    width = 6200, height = 4700, res = 300, background = "white"
  )
  draw_single_ring_tree(ring_key, metadata_column, title_label, legend_title)
  dev.off()

  pdf(
    file.path(out_dir, paste0(stem, ".pdf")),
    width = 20.67, height = 15.67, useDingbats = FALSE
  )
  draw_single_ring_tree(ring_key, metadata_column, title_label, legend_title)
  dev.off()
}

render(
  "gorilla_phylogeny_cluster_ring_radial", "inner_cluster", "gorilla_cluster",
  "gorilla cluster", "Gorilla cluster"
)
render(
  "gorilla_phylogeny_hsa_ring_radial", "middle_hsa", "synteny_group",
  "HSA synteny", "HSA synteny"
)
render(
  "gorilla_phylogeny_nor_ring_radial", "outer_nor", "nor_status",
  "NOR status", "NOR status"
)

cat("Rendered three single-ring phylogenies with", Ntip(tree), "tips each\n")
