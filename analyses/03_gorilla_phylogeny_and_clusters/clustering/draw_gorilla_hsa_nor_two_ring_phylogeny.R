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

hsa_key <- palette_table[palette_table$ring == "middle_hsa", , drop = FALSE]
nor_key <- palette_table[palette_table$ring == "outer_nor", , drop = FALSE]
hsa_colors <- setNames(hsa_key$color, hsa_key$label)
nor_colors <- setNames(nor_key$color, nor_key$label)
tip_hsa_colors <- unname(hsa_colors[meta$synteny_group])
tip_nor_colors <- unname(nor_colors[meta$nor_status])
stopifnot(!anyNA(tip_hsa_colors), !anyNA(tip_nor_colors))

hsa_counts <- as.integer(table(factor(meta$synteny_group, levels = hsa_key$label)))
nor_counts <- as.integer(table(factor(meta$nor_status, levels = nor_key$label)))

ring_geometry <- function() {
  state <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x <- state$xx[seq_len(state$Ntip)]
  y <- state$yy[seq_len(state$Ntip)]
  list(angles = atan2(y, x), tip_radius = median(sqrt(x^2 + y^2)))
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
  symbols(0, 0, circles = inner, inches = FALSE, add = TRUE, fg = "white", lwd = 0.65)
  symbols(0, 0, circles = outer, inches = FALSE, add = TRUE, fg = "white", lwd = 0.65)
}

draw_legend_panel <- function() {
  par(mar = c(0.5, 0.2, 5.0, 0.5), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend(
    x = 0.02, y = 0.93,
    legend = paste0(hsa_key$label, "  (n=", hsa_counts, ")"),
    pch = 22, pt.bg = hsa_key$color, col = "#181818", pt.lwd = 0.8,
    pt.cex = 1.5, bty = "n", cex = 0.84, y.intersp = 1.08,
    title = "Inner ring: HSA synteny"
  )
  legend(
    x = 0.02, y = 0.39,
    legend = paste0(nor_key$label, "  (n=", nor_counts, ")"),
    pch = 22, pt.bg = nor_key$color, col = "#181818", pt.lwd = 0.8,
    pt.cex = 1.7, bty = "n", cex = 0.92, y.intersp = 1.22,
    title = "Outer ring: NOR status"
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

draw_tree <- function() {
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
  draw_annular_ring(tip_hsa_colors, 1.018, 1.063)
  draw_annular_ring(tip_nor_colors, 1.075, 1.122)
  title(
    main = "Gorilla ACRE phylogeny: HSA synteny and NOR status",
    cex.main = 1.43, col.main = "#111111", line = 3.4, font.main = 2
  )
  mtext(
    paste0(
      Ntip(tree),
      " ACREs  |  inner ring: HSA synteny  |  outer ring: NOR status"
    ),
    side = 3, line = 1.15, cex = 0.84, col = "#56636A"
  )
  draw_legend_panel()
}

draw_linear_legend_panel <- function() {
  par(mar = c(0.5, 0.2, 5.0, 0.5), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend(
    x = 0.02, y = 0.93,
    legend = paste0(hsa_key$label, "  (n=", hsa_counts, ")"),
    pch = 22, pt.bg = hsa_key$color, col = "#181818", pt.lwd = 0.8,
    pt.cex = 1.5, bty = "n", cex = 0.84, y.intersp = 1.08,
    title = "First tip strip: HSA synteny"
  )
  legend(
    x = 0.02, y = 0.39,
    legend = paste0(nor_key$label, "  (n=", nor_counts, ")"),
    pch = 22, pt.bg = nor_key$color, col = "#181818", pt.lwd = 0.8,
    pt.cex = 1.7, bty = "n", cex = 0.92, y.intersp = 1.22,
    title = "Second tip strip: NOR status"
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

draw_linear_tree <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.3, 1.35))
  par(mar = c(0.7, 0.6, 5.8, 0.2), family = "sans", xpd = NA, fg = "#111111")
  terminal_edge <- tree$edge[, 2] <= Ntip(tree)
  maximum_depth <- max(node.depth(tree, method = 1)) - 1
  plot(
    tree,
    type = "phylogram", direction = "rightwards",
    use.edge.length = FALSE, show.tip.label = FALSE,
    edge.color = ifelse(terminal_edge, "#333333", "#050505"),
    edge.width = ifelse(terminal_edge, 0.28, 0.78),
    no.margin = FALSE,
    x.lim = c(0, maximum_depth + 90),
    y.lim = c(0.5, Ntip(tree) + 23)
  )
  state <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_x <- max(state$xx[seq_len(state$Ntip)])
  tip_y <- state$yy[seq_len(state$Ntip)]
  gap <- 8
  strip_width <- 23
  first_left <- tip_x + gap
  first_right <- first_left + strip_width
  second_left <- first_right + 5
  second_right <- second_left + strip_width
  rect(
    first_left, tip_y - 0.5, first_right, tip_y + 0.5,
    col = tip_hsa_colors, border = NA
  )
  rect(
    second_left, tip_y - 0.5, second_right, tip_y + 0.5,
    col = tip_nor_colors, border = NA
  )
  rect(first_left, min(tip_y) - 0.5, first_right, max(tip_y) + 0.5, border = "white", lwd = 0.7)
  rect(second_left, min(tip_y) - 0.5, second_right, max(tip_y) + 0.5, border = "white", lwd = 0.7)
  text((first_left + first_right) / 2, Ntip(tree) + 9, "HSA", cex = 0.78, font = 2)
  text((second_left + second_right) / 2, Ntip(tree) + 9, "NOR", cex = 0.78, font = 2)
  title(
    main = "Gorilla ACRE phylogeny: HSA synteny and NOR status",
    cex.main = 1.43, col.main = "#111111", line = 3.5, font.main = 2
  )
  mtext(
    paste0(
      Ntip(tree),
      " ACREs  |  first tip strip: HSA synteny  |  second tip strip: NOR status"
    ),
    side = 3, line = 1.35, cex = 0.84, col = "#56636A"
  )
  draw_linear_legend_panel()
}

stem <- "gorilla_phylogeny_hsa_inner_nor_outer_radial"
ragg::agg_png(
  file.path(out_dir, paste0(stem, ".png")),
  width = 6200, height = 4700, res = 300, background = "white"
)
draw_tree()
dev.off()

pdf(
  file.path(out_dir, paste0(stem, ".pdf")),
  width = 20.67, height = 15.67, useDingbats = FALSE
)
draw_tree()
dev.off()

linear_stem <- "gorilla_phylogeny_hsa_first_nor_second_linear"
ragg::agg_png(
  file.path(out_dir, paste0(linear_stem, ".png")),
  width = 5200, height = 10000, res = 300, background = "white"
)
draw_linear_tree()
dev.off()

pdf(
  file.path(out_dir, paste0(linear_stem, ".pdf")),
  width = 17.33, height = 33.33, useDingbats = FALSE
)
draw_linear_tree()
dev.off()

cat("Rendered HSA-inner/NOR-outer phylogeny with", Ntip(tree), "tips\n")
