#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "phylogeny_origin_acro1_full"
results <- file.path(root, "results")
figures <- file.path(root, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

arch <- read.delim(file.path(results, "orthologous_locus_architecture.tsv"), check.names = FALSE)
ev <- read.delim(file.path(results, "integrated_evidence_matrix.tsv"), check.names = FALSE)

display_names <- c(
  Human = "Human", Chimpanzee = "Chimpanzee", Bonobo = "Bonobo", Gorilla = "Gorilla",
  Sumatran_orangutan = "Sumatran orangutan", Bornean_orangutan = "Bornean orangutan",
  Siamang_hap1 = "Siamang haplotype 1", Siamang_hap2 = "Siamang haplotype 2",
  Hoolock_gibbon = "Hoolock gibbon", Nomascus_gibbon = "Nomascus gibbon",
  Rhesus_macaque = "Rhesus macaque", Marmoset = "Common marmoset",
  Golden_lion_tamarin = "Golden lion tamarin", Squirrel_monkey = "Squirrel monkey"
)
plot_order <- rev(c(
  "Marmoset", "Golden_lion_tamarin", "Squirrel_monkey", "Rhesus_macaque",
  "Nomascus_gibbon", "Hoolock_gibbon", "Siamang_hap2", "Siamang_hap1",
  "Bornean_orangutan", "Sumatran_orangutan", "Gorilla", "Bonobo", "Chimpanzee", "Human"
))
arch$display <- factor(display_names[arch$assembly_key], levels = display_names[plot_order])

metrics <- rbind(
  data.frame(display = arch$display, metric = "Inferred ortholog interval", value = arch$ortholog_locus_span_bp / 1000),
  data.frame(display = arch$display, metric = "Recognizable ACRO1-like sequence", value = arch$ACRO1_bp_in_inferred_locus / 1000)
)
metrics$metric <- factor(metrics$metric, levels = c("Inferred ortholog interval", "Recognizable ACRO1-like sequence"))
p1 <- ggplot(metrics, aes(value, display, fill = metric)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("#b8bec7", "#6a3d9a")) +
  labs(x = "Length (kb)", y = NULL, title = "HSA19-syntenic ACRO1 locus across primates") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
ggsave(file.path(figures, "HSA19_ortholog_sizes.png"), p1, width = 10.5, height = 6.3, dpi = 300)
ggsave(file.path(figures, "HSA19_ortholog_sizes.pdf"), p1, width = 10.5, height = 6.3)

state_cols <- c(
  empty_or_unrecognizable = "#f2f2f2", isolated_ACRO1_fragment = "#d9c2e9",
  short_ACRO1_like_array = "#b99ad6", short_partial_ACRE = "#9b72bd",
  moderate_ACRO1_like_array = "#8055a5", ACRO1_array_with_ACRE_components = "#663b8c",
  complete_ACRE_like = "#4b266f", expanded_ACRO1_like_array = "#30104e"
)
ev$display <- factor(display_names[ev$assembly_key], levels = display_names[plot_order])
ev$locus_state <- factor(ev$locus_state, levels = names(state_cols))
p2 <- ggplot(ev, aes(x = "HSA19 ortholog", y = display, fill = locus_state)) +
  geom_tile(color = "white", linewidth = 0.7) +
  scale_fill_manual(values = state_cols, drop = TRUE, name = "Sequence state") +
  labs(x = NULL, y = NULL, title = "ACRO1/ACRE state at the orthologous locus") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"), plot.title = element_text(face = "bold"), legend.position = "right")
ggsave(file.path(figures, "HSA19_state_matrix.png"), p2, width = 9.5, height = 6.3, dpi = 300)
ggsave(file.path(figures, "HSA19_state_matrix.pdf"), p2, width = 9.5, height = 6.3)

plot_tree <- function(tree_path, title, out_name, width = 8.5, height = 7.0) {
  if (!file.exists(tree_path)) return(invisible(NULL))
  tr <- read.tree(tree_path)
  pdf(file.path(figures, paste0(out_name, ".pdf")), width = width, height = height)
  par(mar = c(1, 1, 3, 1))
  plot(tr, type = "phylogram", cex = 0.72, no.margin = FALSE, label.offset = max(node.depth.edgelength(tr)) * 0.012)
  nodelabels(tr$node.label, frame = "none", cex = 0.55, adj = c(1.2, -0.2))
  title(title, font.main = 2)
  dev.off()
  png(file.path(figures, paste0(out_name, ".png")), width = width, height = height, units = "in", res = 300)
  par(mar = c(1, 1, 3, 1))
  plot(tr, type = "phylogram", cex = 0.72, no.margin = FALSE, label.offset = max(node.depth.edgelength(tr)) * 0.012)
  nodelabels(tr$node.label, frame = "none", cex = 0.55, adj = c(1.2, -0.2))
  title(title, font.main = 2)
  dev.off()
}

tree_dir <- file.path(results, "phylogenies")
plot_tree(file.path(tree_dir, "unique_flank_tree.treefile"), "Unique-flank maximum-likelihood tree", "unique_flank_tree")
plot_tree(file.path(tree_dir, "ACRE_backbone_tree.treefile"), "Alignable ACRE-backbone maximum-likelihood tree", "ACRE_backbone_tree", 7.5, 5.5)
plot_tree(file.path(tree_dir, "ACRO1_locus_consensus_tree.treefile"), "ACRO1 locus-consensus tree (exploratory)", "ACRO1_locus_consensus_tree")

