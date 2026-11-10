#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Could not determine plotting-script path")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
project <- dirname(dirname(script_path))

tree <- read.tree(file.path(project, "results", "panel_sequences.contree"))
manifest <- read.delim(
  file.path(project, "inputs", "panel_manifest.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
rownames(manifest) <- manifest$tip_label

marmoset_tip <- manifest$tip_label[manifest$species_key == "marmoset-partial"]
if (length(marmoset_tip) != 1) stop("Expected one marmoset partial tip")
tree <- root(tree, outgroup = marmoset_tip, resolve.root = TRUE)
write.tree(tree, file.path(project, "results", "panel_sequences.rooted.contree"))

# The partial marmoset tract has a much longer terminal branch than the
# complete ingroup units. Preserve the inferred tree on disk, but shorten that
# single terminal branch in the figures so the ingroup topology remains
# legible at publication scale.
plot_tree <- tree
marmoset_index <- match(marmoset_tip, plot_tree$tip.label)
marmoset_edge <- which(plot_tree$edge[, 2] == marmoset_index)
marmoset_original_length <- plot_tree$edge.length[marmoset_edge]
marmoset_display_length <- min(marmoset_original_length, 0.06)
plot_tree$edge.length[marmoset_edge] <- marmoset_display_length
display_note <- "Display note: marmoset terminal branch shortened for readability; topology and all other branch lengths are unchanged."
write.table(
  data.frame(
    tip = marmoset_tip,
    original_terminal_branch = marmoset_original_length,
    displayed_terminal_branch = marmoset_display_length
  ),
  file.path(project, "results", "figure_branch_scaling.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

meta <- manifest[tree$tip.label, , drop = FALSE]
display_group <- meta$species_key
display_group[display_group %in% c("siamang-hap1", "siamang-hap2")] <- "siamang"

palette <- c(
  "human" = "#4C78A8", "gorilla" = "#F58518",
  "bonobo" = "#54A24B", "chimpanzee" = "#E45756",
  "sumatran-orangutan" = "#B279A2", "bornean-orangutan" = "#FF9DA6",
  "siamang" = "#1B9E77", "marmoset-partial" = "#666666"
)
branch_color <- "#000000"
species_display <- c(
  "human" = "Human", "gorilla" = "Gorilla", "bonobo" = "Bonobo",
  "chimpanzee" = "Chimpanzee", "sumatran-orangutan" = "Sumatran orangutan",
  "bornean-orangutan" = "Bornean orangutan", "siamang" = "Siamang",
  "marmoset-partial" = "Marmoset (partial tract)"
)

tip_colors <- unname(palette[display_group])
tip_pch <- ifelse(display_group == "marmoset-partial", 17, 16)
legend_order <- names(species_display)
legend_counts <- table(display_group)
legend_text <- vapply(
  legend_order,
  function(group) paste0(species_display[[group]], " (n=", unname(legend_counts[[group]]), ")"),
  character(1)
)
max_depth <- max(node.depth.edgelength(plot_tree))

tip_text <- paste(
  meta$species_display,
  gsub("_", " ", meta$chromosome, fixed = TRUE),
  format(meta$start, big.mark = ",", scientific = FALSE, trim = TRUE),
  sep = " | "
)

draw_species <- function(path, type = c("phylogram", "fan"), device = c("png", "pdf")) {
  type <- match.arg(type); device <- match.arg(device)
  if (device == "png") {
    if (type == "fan") png(path, width = 10000, height = 10000, res = 220)
    else png(path, width = 8000, height = 14000, res = 220)
    edge_width <- 8.0; point_cex <- ifelse(type == "fan", 19.35, 21.15); legend_cex <- 1.2
  } else {
    if (type == "fan") pdf(path, width = 34, height = 34, useDingbats = FALSE)
    else pdf(path, width = 24, height = 40, useDingbats = FALSE)
    edge_width <- 4.0; point_cex <- ifelse(type == "fan", 8.1, 9.0); legend_cex <- 0.85
  }
  par(mar = if (type == "fan") c(5, 5, 7, 5) else c(3, 3, 7, 3), xpd = NA)
  plot(
    plot_tree, type = type, show.tip.label = FALSE, edge.width = edge_width,
    edge.color = branch_color,
    no.margin = FALSE,
    x.lim = if (type == "phylogram") c(0, max_depth * 1.42) else NULL
  )
  tiplabels(pch = tip_pch, col = tip_colors, bg = tip_colors, cex = point_cex)
  title("Cross-primate ACRO-unit phylogeny: 50 units per species", cex.main = 1.5, font.main = 2)
  mtext(display_note, side = 1, line = 1.2, cex = ifelse(device == "png", 0.82, 0.72), col = "#333333")
  legend(
    if (type == "fan") "topleft" else "topright",
    legend = legend_text, col = unname(palette[legend_order]),
    pch = c(rep(16, 7), 17), pt.cex = 3.1, cex = legend_cex,
    bty = "n", y.intersp = 1.1
  )
  dev.off()
}

draw_labeled <- function(path, type = c("phylogram", "fan"), device = c("png", "pdf")) {
  type <- match.arg(type); device <- match.arg(device)
  display_tree <- plot_tree; display_tree$tip.label <- tip_text
  if (device == "png") {
    if (type == "fan") png(path, width = 12000, height = 12000, res = 200)
    else png(path, width = 10000, height = 18000, res = 200)
    edge_width <- 7.0; point_cex <- ifelse(type == "fan", 9.45, 10.35)
    label_cex <- ifelse(type == "fan", 0.36, 0.48); legend_cex <- 0.95
  } else {
    if (type == "fan") pdf(path, width = 46, height = 46, useDingbats = FALSE)
    else pdf(path, width = 36, height = 54, useDingbats = FALSE)
    edge_width <- 3.6; point_cex <- ifelse(type == "fan", 4.5, 5.22)
    label_cex <- ifelse(type == "fan", 0.2, 0.26); legend_cex <- 0.65
  }
  par(mar = if (type == "fan") c(8, 8, 8, 8) else c(3, 3, 8, 25), xpd = NA)
  plot(
    display_tree, type = type, show.tip.label = TRUE,
    tip.color = tip_colors, cex = label_cex,
    label.offset = ifelse(type == "fan", 0.0017, 0.0008),
    edge.width = edge_width, edge.color = branch_color, no.margin = FALSE,
    x.lim = if (type == "phylogram") c(0, max_depth * 2.05) else NULL
  )
  tiplabels(pch = tip_pch, col = tip_colors, bg = tip_colors, cex = point_cex)
  title("Cross-primate ACRO-unit phylogeny with genomic locations", cex.main = 1.45, font.main = 2)
  mtext(display_note, side = 1, line = 1.2, cex = ifelse(device == "png", 0.72, 0.58), col = "#333333")
  legend(
    if (type == "fan") "topleft" else "topright",
    legend = legend_text, col = unname(palette[legend_order]),
    pch = c(rep(16, 7), 17), pt.cex = 2.8, cex = legend_cex,
    bty = "n", y.intersp = 1.1
  )
  dev.off()
}

for (device in c("png", "pdf")) {
  draw_species(file.path(project, "figures", paste0("phylogeny_50_species_rectangular.", device)), "phylogram", device)
  draw_species(file.path(project, "figures", paste0("phylogeny_50_species_radial.", device)), "fan", device)
  draw_labeled(file.path(project, "figures", paste0("phylogeny_50_species_chromosome_rectangular.", device)), "phylogram", device)
  draw_labeled(file.path(project, "figures", paste0("phylogeny_50_species_chromosome_radial.", device)), "fan", device)
}

audit_groups <- c(
  "human", "gorilla", "bonobo", "chimpanzee",
  "sumatran-orangutan", "bornean-orangutan", "siamang"
)
audit <- do.call(rbind, lapply(audit_groups, function(group) {
  labels <- tree$tip.label[display_group == group]
  data.frame(group = group, n_tips = length(labels), monophyletic = is.monophyletic(tree, labels))
}))
write.table(audit, file.path(project, "results", "species_monophyly.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
