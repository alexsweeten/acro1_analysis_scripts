#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
root_dir <- dirname(out_dir)

tree <- read.tree(file.path(root_dir, "acro_gorilla_only", "gorilla_acres.treefile"))
meta <- read.delim(
  gzfile(file.path(out_dir, "gorilla_cluster_assignments.tsv.gz")),
  stringsAsFactors = FALSE
)
stopifnot(
  !anyDuplicated(meta$tip_id),
  setequal(tree$tip.label, meta$tip_id),
  all(meta$gorilla_cluster %in% paste0("G", 1:6))
)

# IQ-TREE reports an unrooted topology. Select the internal node with the
# smallest maximum edge-count distance to any tip only to balance the radial
# display; this does not change the unrooted splits.
unit_tree <- tree
unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))
node_distances <- dist.nodes(unit_tree)
internal_nodes <- (Ntip(unit_tree) + 1):(Ntip(unit_tree) + Nnode(unit_tree))
eccentricity <- apply(
  node_distances[internal_nodes, seq_len(Ntip(unit_tree)), drop = FALSE],
  1,
  max
)
display_root <- internal_nodes[which.min(eccentricity)]
tree <- ladderize(root(tree, node = display_root, resolve.root = TRUE), right = TRUE)
meta <- meta[match(tree$tip.label, meta$tip_id), , drop = FALSE]

natural_hsa_order <- function(values) {
  number <- suppressWarnings(as.numeric(sub("^hsa([0-9]+).*$", "\\1", values)))
  suffix <- sub("^hsa[0-9]+", "", values)
  is_y <- grepl("^hsaY$", values, ignore.case = TRUE)
  table <- data.frame(
    value = values,
    group = ifelse(!is.na(number), 1, ifelse(is_y, 2, 3)),
    number = ifelse(is.na(number), Inf, number),
    suffix = suffix,
    stringsAsFactors = FALSE
  )
  table$value[order(table$group, table$number, table$suffix, table$value)]
}

cluster_levels <- paste0("G", 1:6)
cluster_colors <- c(
  G1 = "#4E79A7", G2 = "#F28E2B", G3 = "#59A14F",
  G4 = "#E15759", G5 = "#B07AA1", G6 = "#76B7B2"
)
hsa_levels <- natural_hsa_order(unique(meta$synteny_group))
hsa_palette <- c(
  "#3B82C4", "#E6862A", "#8C8C8C", "#E05A5A", "#B0A52B",
  "#45A778", "#17A6B5", "#C674B2", "#2D8C87", "#6D4C9F", "#D45B93"
)
stopifnot(length(hsa_levels) <= length(hsa_palette))
hsa_colors <- setNames(hsa_palette[seq_along(hsa_levels)], hsa_levels)
nor_levels <- c("NOR+", "NOR-", "non-acro")
nor_colors <- c("NOR+" = "#D81B60", "NOR-" = "#1769AA", "non-acro" = "#A7ADB2")

tip_cluster_colors <- unname(cluster_colors[meta$gorilla_cluster])
tip_hsa_colors <- unname(hsa_colors[meta$synteny_group])
tip_nor_colors <- unname(nor_colors[meta$nor_status])
stopifnot(
  !anyNA(tip_cluster_colors), !anyNA(tip_hsa_colors), !anyNA(tip_nor_colors)
)

cluster_counts <- table(factor(meta$gorilla_cluster, levels = cluster_levels))
hsa_counts <- table(factor(meta$synteny_group, levels = hsa_levels))
nor_counts <- table(factor(meta$nor_status, levels = nor_levels))

ring_geometry <- function() {
  state <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x <- state$xx[seq_len(state$Ntip)]
  y <- state$yy[seq_len(state$Ntip)]
  list(
    angles = atan2(y, x),
    tip_radius = median(sqrt(x^2 + y^2))
  )
}

draw_annular_ring <- function(colors, inner_multiplier, outer_multiplier) {
  geometry <- ring_geometry()
  angles <- geometry$angles
  inner <- geometry$tip_radius * inner_multiplier
  outer <- geometry$tip_radius * outer_multiplier
  half_width <- pi / length(angles) * 1.025
  for (index in seq_along(angles)) {
    theta <- seq(angles[index] - half_width, angles[index] + half_width, length.out = 5)
    polygon(
      c(inner * cos(theta), rev(outer * cos(theta))),
      c(inner * sin(theta), rev(outer * sin(theta))),
      col = colors[index], border = NA
    )
  }
  symbols(0, 0, circles = inner, inches = FALSE, add = TRUE, fg = "white", lwd = 0.6)
  symbols(0, 0, circles = outer, inches = FALSE, add = TRUE, fg = "white", lwd = 0.6)
}

draw_legend_panel <- function() {
  par(mar = c(0.5, 0.2, 4.5, 0.5), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend(
    x = 0.02, y = 0.97,
    legend = paste0(cluster_levels, "  (n=", as.integer(cluster_counts), ")"),
    pch = 22, pt.bg = unname(cluster_colors[cluster_levels]), col = "#181818",
    pt.lwd = 0.7, pt.cex = 1.55, bty = "n", cex = 0.86, y.intersp = 1.15,
    title = "Inner ring: gorilla cluster"
  )
  legend(
    x = 0.02, y = 0.66,
    legend = paste0(hsa_levels, "  (n=", as.integer(hsa_counts), ")"),
    pch = 22, pt.bg = unname(hsa_colors[hsa_levels]), col = "#181818",
    pt.lwd = 0.7, pt.cex = 1.45, bty = "n", cex = 0.80, y.intersp = 1.05,
    title = "Middle ring: HSA synteny"
  )
  legend(
    x = 0.02, y = 0.20,
    legend = paste0(nor_levels, "  (n=", as.integer(nor_counts), ")"),
    pch = 22, pt.bg = unname(nor_colors[nor_levels]), col = "#181818",
    pt.lwd = 0.7, pt.cex = 1.55, bty = "n", cex = 0.85, y.intersp = 1.15,
    title = "Outer ring: NOR status"
  )
  text(
    0.02, 0.045,
    "Black internal branches\nBranch lengths suppressed for legibility",
    adj = c(0, 0), cex = 0.70, col = "#5F696E"
  )
}

draw_tree <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.3, 1.35))
  par(mar = c(0.4, 0.4, 5.6, 0.2), family = "sans", xpd = NA, fg = "#111111")

  terminal_edge <- tree$edge[, 2] <= Ntip(tree)
  edge_colors <- ifelse(terminal_edge, "#333333", "#050505")
  edge_widths <- ifelse(terminal_edge, 0.30, 0.82)
  plot(
    tree,
    type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
    rotate.tree = 4, edge.color = edge_colors, edge.width = edge_widths,
    no.margin = FALSE, x.lim = c(-1.24, 1.24), y.lim = c(-1.24, 1.24)
  )

  draw_annular_ring(tip_cluster_colors, 1.016, 1.055)
  draw_annular_ring(tip_hsa_colors, 1.066, 1.108)
  draw_annular_ring(tip_nor_colors, 1.119, 1.163)

  title(
    main = "Gorilla ACRE phylogeny with cluster, chromosome, and NOR annotations",
    cex.main = 1.38, col.main = "#111111", line = 3.4, font.main = 2
  )
  mtext(
    paste0(
      Ntip(tree),
      " ACREs  |  rings from inside out: G1-G6 cluster, HSA synteny, NOR status"
    ),
    side = 3, line = 1.15, cex = 0.82, col = "#56636A"
  )
  draw_legend_panel()
}

ragg::agg_png(
  file.path(out_dir, "gorilla_cluster_hsa_nor_phylogeny_radial.png"),
  width = 6200, height = 4700, res = 300, background = "white"
)
draw_tree()
dev.off()

pdf(
  file.path(out_dir, "gorilla_cluster_hsa_nor_phylogeny_radial.pdf"),
  width = 20.67, height = 15.67, useDingbats = FALSE
)
draw_tree()
dev.off()

annotation_table <- meta[, c(
  "tip_id", "acre_id", "array_id", "original_chrom", "synteny_group",
  "nor_status", "gorilla_cluster"
)]
write.table(
  annotation_table,
  file.path(out_dir, "gorilla_cluster_hsa_nor_phylogeny_annotations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.tree(
  tree,
  file = file.path(out_dir, "gorilla_cluster_hsa_nor_display_rooted.treefile")
)

palette_table <- rbind(
  data.frame(
    ring = "inner_cluster", label = cluster_levels,
    color = unname(cluster_colors[cluster_levels]), count = as.integer(cluster_counts)
  ),
  data.frame(
    ring = "middle_hsa", label = hsa_levels,
    color = unname(hsa_colors[hsa_levels]), count = as.integer(hsa_counts)
  ),
  data.frame(
    ring = "outer_nor", label = nor_levels,
    color = unname(nor_colors[nor_levels]), count = as.integer(nor_counts)
  )
)
write.table(
  palette_table,
  file.path(out_dir, "gorilla_cluster_hsa_nor_phylogeny_palette.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(
  "Rendered", Ntip(tree), "tips with cluster, HSA, and NOR rings; display root node",
  display_root, "\n"
)
