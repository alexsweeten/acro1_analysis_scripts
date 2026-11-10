#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
out_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
tree <- read.tree(file.path(out_dir, "gorilla_acres.treefile"))
meta <- read.delim(file.path(out_dir, "gorilla_tip_metadata.tsv"), stringsAsFactors = FALSE)
stopifnot(setequal(tree$tip.label, meta$tip_id))
meta <- meta[match(tree$tip.label, meta$tip_id), , drop = FALSE]

children <- split(tree$edge[, 2], tree$edge[, 1])
descendants <- vector("list", Ntip(tree) + Nnode(tree))
for (tip in seq_len(Ntip(tree))) descendants[[tip]] <- tip
resolve_descendants <- function(node) {
  if (node <= Ntip(tree)) return(descendants[[node]])
  descendants[[node]] <<- unlist(
    lapply(children[[as.character(node)]], resolve_descendants), use.names = FALSE
  )
  descendants[[node]]
}
root <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]
invisible(resolve_descendants(root))
all_tips <- seq_len(Ntip(tree))

parse_support <- function(label) {
  if (is.na(label) || label == "") return(NA_real_)
  suppressWarnings(as.numeric(strsplit(label, "/", fixed = TRUE)[[1]][1]))
}

smallest_unrooted_side <- function(group_tips) {
  best_size <- length(all_tips)
  best_edge <- NA_integer_
  for (edge_index in seq_len(nrow(tree$edge))) {
    child <- tree$edge[edge_index, 2]
    side <- descendants[[child]]
    other_size <- length(all_tips) - length(side)
    if (all(group_tips %in% side) && length(side) < best_size) {
      best_size <- length(side)
      best_edge <- edge_index
    }
    if (all(!group_tips %in% side) && other_size < best_size) {
      best_size <- other_size
      best_edge <- -edge_index
    }
  }
  support <- NA_real_
  if (!is.na(best_edge) && best_edge > 0) {
    child <- tree$edge[best_edge, 2]
    if (child > Ntip(tree)) support <- parse_support(tree$node.label[child - Ntip(tree)])
  }
  list(size = best_size, support = support)
}

patristic <- cophenetic.phylo(tree)
nearest_index <- apply(patristic + diag(Inf, Ntip(tree)), 1, which.min)

group_rows <- list()
row_index <- 0
for (trait in c("de_novo_cluster", "synteny_group")) {
  values <- as.character(meta[[trait]])
  for (group in sort(unique(values))) {
    indices <- which(values == group)
    envelope <- smallest_unrooted_side(indices)
    within <- if (length(indices) > 1) {
      mean(patristic[indices, indices][upper.tri(patristic[indices, indices])])
    } else {
      NA_real_
    }
    outside <- setdiff(all_tips, indices)
    row_index <- row_index + 1
    group_rows[[row_index]] <- data.frame(
      trait = trait,
      group = group,
      tips = length(indices),
      smallest_unrooted_envelope_tips = envelope$size,
      envelope_purity = length(indices) / envelope$size,
      monophyletic_unrooted = envelope$size == length(indices),
      envelope_sh_alrt = envelope$support,
      mean_within_patristic_distance = within,
      mean_to_other_tips_patristic_distance = mean(patristic[indices, outside]),
      nearest_neighbor_same_group_fraction = mean(values[nearest_index[indices]] == group),
      stringsAsFactors = FALSE
    )
  }
}
group_stats <- do.call(rbind, group_rows)
write.table(
  group_stats, file.path(out_dir, "tree_group_clade_stats.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

fitch_score <- function(states) {
  node_sets <- vector("list", Ntip(tree) + Nnode(tree))
  for (tip in seq_len(Ntip(tree))) node_sets[[tip]] <- states[tip]
  score <- 0L
  recurse <- function(node) {
    if (node <= Ntip(tree)) return(node_sets[[node]])
    child_sets <- lapply(children[[as.character(node)]], recurse)
    current <- child_sets[[1]]
    if (length(child_sets) > 1) {
      for (index in 2:length(child_sets)) {
        common <- intersect(current, child_sets[[index]])
        if (length(common)) {
          current <- common
        } else {
          current <- union(current, child_sets[[index]])
          score <<- score + 1L
        }
      }
    }
    current
  }
  recurse(root)
  score
}

set.seed(1)
signal_rows <- list()
for (trait in c("de_novo_cluster", "synteny_group")) {
  states <- as.character(meta[[trait]])
  observed <- fitch_score(states)
  permutations <- replicate(999, fitch_score(sample(states, replace = FALSE)))
  signal_rows[[trait]] <- data.frame(
    trait = trait,
    states = length(unique(states)),
    observed_fitch_changes = observed,
    permuted_mean_changes = mean(permutations),
    permuted_sd_changes = sd(permutations),
    lower_tail_permutation_p = (1 + sum(permutations <= observed)) / 1000,
    nearest_neighbor_same_state_fraction = mean(states[nearest_index] == states),
    stringsAsFactors = FALSE
  )
}
write.table(
  do.call(rbind, signal_rows), file.path(out_dir, "tree_phylogenetic_signal.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

node_support <- vapply(tree$node.label, parse_support, numeric(1))
support_summary <- data.frame(
  tips = Ntip(tree),
  internal_nodes = Nnode(tree),
  nodes_with_sh_alrt = sum(!is.na(node_support)),
  median_sh_alrt = median(node_support, na.rm = TRUE),
  fraction_sh_alrt_ge_80 = mean(node_support >= 80, na.rm = TRUE),
  total_branch_length = sum(tree$edge.length),
  stringsAsFactors = FALSE
)
write.table(
  support_summary, file.path(out_dir, "tree_support_summary.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

cat("Wrote tree group, phylogenetic-signal, and branch-support summaries\n")
