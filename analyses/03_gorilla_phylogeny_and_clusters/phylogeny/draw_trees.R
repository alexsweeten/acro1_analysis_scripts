#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))

tree <- read.tree(file.path(out_dir, "gorilla_acres.treefile"))
# IQ-TREE writes an arbitrary origin for this unrooted topology. Choose the
# internal node with the smallest maximum edge-count distance to all tips as a
# neutral display root. This changes no unrooted split or support value.
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
meta <- read.delim(file.path(out_dir, "gorilla_tip_metadata.tsv"), stringsAsFactors = FALSE)
stopifnot(setequal(tree$tip.label, meta$tip_id))
meta <- meta[match(tree$tip.label, meta$tip_id), , drop = FALSE]

cluster_table <- read.delim(file.path(out_dir, "cluster_palette.tsv"), stringsAsFactors = FALSE)
cluster_palette <- setNames(cluster_table$color, as.character(cluster_table$de_novo_cluster))
chromosome_table <- read.delim(file.path(out_dir, "chromosome_palette.tsv"), stringsAsFactors = FALSE)
chromosome_palette <- setNames(chromosome_table$color, chromosome_table$synteny_group)

homogeneous_edge_colors <- function(tree, categories, palette) {
  children <- split(tree$edge[, 2], tree$edge[, 1])
  node_state <- rep(NA_character_, Ntip(tree) + Nnode(tree))
  node_state[seq_len(Ntip(tree))] <- categories
  resolve <- function(node) {
    if (node <= Ntip(tree)) return(node_state[node])
    child_states <- vapply(children[[as.character(node)]], resolve, character(1))
    if (all(!is.na(child_states)) && length(unique(child_states)) == 1) {
      node_state[node] <<- child_states[1]
    }
    node_state[node]
  }
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
  resolve(root)
  colors <- rep("#C7CDD1", nrow(tree$edge))
  resolved <- node_state[tree$edge[, 2]]
  colors[!is.na(resolved)] <- unname(palette[resolved[!is.na(resolved)]])
  colors
}

draw_tree <- function(color_by = c("cluster", "chromosome"), layout_type = c("radial", "rectangular")) {
  color_by <- match.arg(color_by)
  layout_type <- match.arg(layout_type)
  if (color_by == "cluster") {
    categories <- as.character(meta$de_novo_cluster)
    palette <- cluster_palette
    counts <- table(factor(categories, levels = names(palette)))
    legend_text <- paste0("G", names(palette), " (n=", as.integer(counts), ")")
    legend_title <- "Gorilla-only cluster"
    title_text <- paste0("All 1,232 gorilla ACREs colored by de novo cluster (k = ", length(palette), ")")
  } else {
    categories <- meta$synteny_group
    palette <- chromosome_palette
    counts <- table(factor(categories, levels = names(palette)))
    legend_text <- paste0(names(palette), " (n=", as.integer(counts), ")")
    legend_title <- "Chromosome (HSA synteny)"
    title_text <- "All 1,232 gorilla ACREs colored by chromosome"
  }
  tip_colors <- unname(palette[categories])
  edge_colors <- homogeneous_edge_colors(tree, categories, palette)
  terminal <- tree$edge[, 2] <= Ntip(tree)

  layout(matrix(c(1, 2), nrow = 1), widths = c(5.4, 1.45))
  par(mar = c(1.2, 1.2, 5.3, 0.5), family = "sans", fg = "#263238")
  if (layout_type == "radial") {
    plot(
      tree, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
      edge.color = edge_colors, edge.width = ifelse(terminal, 1.15, 0.48),
      no.margin = FALSE, rotate.tree = 7
    )
    tiplabels(pch = 21, bg = tip_colors, col = "white", lwd = 0.15, cex = 0.25)
    subtitle <- "Unrooted ML topology, center-rooted for display; homogeneous descendant branches inherit color"
  } else {
    plot(
      tree, type = "phylogram", use.edge.length = FALSE, show.tip.label = FALSE,
      edge.color = edge_colors, edge.width = ifelse(terminal, 1.0, 0.42),
      no.margin = FALSE
    )
    tiplabels(pch = 21, bg = tip_colors, col = "white", lwd = 0.12, cex = 0.20)
    subtitle <- "Unrooted ML topology, center-rooted only for a balanced rectangular display"
  }
  title(main = title_text, cex.main = 1.18, col.main = "#17202A", line = 3.5)
  mtext(subtitle, side = 3, line = 1.0, cex = 0.72, col = "#546E7A")

  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  legend(
    "center", legend = legend_text, pt.bg = unname(palette), col = "white",
    pch = 21, pt.cex = 1.3, bty = "n", cex = 0.86,
    title = legend_title, xjust = 0.5, yjust = 0.5
  )
}

render <- function(stem, color_by, layout_type, width, height, pdf_width, pdf_height) {
  ragg::agg_png(
    file.path(out_dir, paste0(stem, ".png")), width = width, height = height,
    res = 220, background = "white"
  )
  draw_tree(color_by, layout_type)
  dev.off()
  pdf(
    file.path(out_dir, paste0(stem, ".pdf")), width = pdf_width,
    height = pdf_height, useDingbats = FALSE
  )
  draw_tree(color_by, layout_type)
  dev.off()
}

render("tree_gorilla_clusters_radial", "cluster", "radial", 5000, 4200, 25, 21)
render("tree_gorilla_clusters_rectangular", "cluster", "rectangular", 4300, 6200, 21.5, 31)
render("tree_gorilla_chromosomes_radial", "chromosome", "radial", 5000, 4200, 25, 21)
render("tree_gorilla_chromosomes_rectangular", "chromosome", "rectangular", 4300, 6200, 21.5, 31)

cat("Rendered four all-gorilla tree views as PNG and vector PDF\n")
