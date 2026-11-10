#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: filter_reciprocal_alignments.R CHROMOSOME_ROOT")
root <- normalizePath(args[[1]], mustWork = TRUE)
config <- read.delim(file.path(root, "species_config.tsv"), stringsAsFactors = FALSE)
query_ids <- config$species_id[config$species_id != "human"]

min_mapq <- 30
min_block_bp <- 5000
min_identity <- 0.95
min_endpoint_match_coverage <- 0.80
min_reciprocal_overlap <- 0.80

empty_paf <- function() {
  data.frame(query_name = character(), query_length = numeric(), query_start = numeric(),
             query_end = numeric(), strand = character(), target_name = character(),
             target_length = numeric(), target_start = numeric(), target_end = numeric(),
             matching_bases = numeric(), block_length = numeric(), mapq = numeric(),
             tags = character(), stringsAsFactors = FALSE)
}

read_paf <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(empty_paf())
  x <- read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  fill = TRUE, stringsAsFactors = FALSE)
  core <- x[, 1:12]
  names(core) <- names(empty_paf())[1:12]
  if (ncol(x) > 12) {
    core$tags <- apply(x[, 13:ncol(x), drop = FALSE], 1, function(z) {
      paste(z[!is.na(z) & nzchar(z)], collapse = ";")
    })
  } else {
    core$tags <- ""
  }
  core
}

filter_candidates <- function(x) {
  if (!nrow(x)) return(x)
  x$identity <- x$matching_bases / x$block_length
  x$target_span <- x$target_end - x$target_start
  x$query_span <- x$query_end - x$query_start
  x$target_match_coverage <- x$matching_bases / x$target_span
  x$query_match_coverage <- x$matching_bases / x$query_span
  x$is_primary <- grepl("(^|;)tp:A:P($|;)", x$tags)
  x[x$is_primary & x$mapq >= min_mapq & x$block_length >= min_block_bp &
      x$identity >= min_identity &
      x$target_match_coverage >= min_endpoint_match_coverage &
      x$query_match_coverage >= min_endpoint_match_coverage, ]
}

overlap_fraction <- function(a_start, a_end, b_start, b_end) {
  pmax(0, pmin(a_end, b_end) - pmax(a_start, b_start)) /
    pmax(a_end - a_start, b_end - b_start)
}

dir.create(file.path(root, "alignments"), recursive = TRUE, showWarnings = FALSE)
all_blocks <- list()
all_summary <- list()

for (query_id in query_ids) {
  forward_path <- file.path(root, "raw_alignments", paste0("human_target_vs_", query_id, ".paf"))
  reverse_path <- file.path(root, "raw_alignments", paste0(query_id, "_target_vs_human.paf"))
  forward_all <- read_paf(forward_path)
  reverse_all <- read_paf(reverse_path)
  forward <- filter_candidates(forward_all)
  reverse <- filter_candidates(reverse_all)

  if (nrow(forward)) {
    forward$human_start <- forward$target_start
    forward$human_end <- forward$target_end
    forward$primate_start <- forward$query_start
    forward$primate_end <- forward$query_end
    forward$forward_id <- seq_len(nrow(forward))
  }
  if (nrow(reverse)) {
    reverse$human_start <- reverse$query_start
    reverse$human_end <- reverse$query_end
    reverse$primate_start <- reverse$target_start
    reverse$primate_end <- reverse$target_end
    reverse$reverse_id <- seq_len(nrow(reverse))
  }

  pairs <- data.frame()
  if (nrow(forward) && nrow(reverse)) {
    candidates <- list()
    k <- 0L
    for (i in seq_len(nrow(forward))) {
      for (j in seq_len(nrow(reverse))) {
        if (forward$strand[[i]] != reverse$strand[[j]]) next
        human_overlap <- overlap_fraction(forward$human_start[[i]], forward$human_end[[i]],
                                          reverse$human_start[[j]], reverse$human_end[[j]])
        primate_overlap <- overlap_fraction(forward$primate_start[[i]], forward$primate_end[[i]],
                                            reverse$primate_start[[j]], reverse$primate_end[[j]])
        reciprocal_score <- min(human_overlap, primate_overlap)
        if (reciprocal_score <= 0) next
        k <- k + 1L
        candidates[[k]] <- data.frame(forward_id = i, reverse_id = j,
                                      human_reciprocal_overlap = human_overlap,
                                      primate_reciprocal_overlap = primate_overlap,
                                      reciprocal_score = reciprocal_score)
      }
    }
    if (length(candidates)) {
      pairs <- do.call(rbind, candidates)
      best_reverse_for_forward <- tapply(seq_len(nrow(pairs)), pairs$forward_id,
                                         function(ix) ix[[which.max(pairs$reciprocal_score[ix])]])
      best_forward_for_reverse <- tapply(seq_len(nrow(pairs)), pairs$reverse_id,
                                         function(ix) ix[[which.max(pairs$reciprocal_score[ix])]])
      mutual_rows <- intersect(as.integer(best_reverse_for_forward), as.integer(best_forward_for_reverse))
      pairs <- pairs[mutual_rows, ]
      pairs <- pairs[pairs$reciprocal_score >= min_reciprocal_overlap, ]
    }
  }

  kept <- if (nrow(pairs)) forward[pairs$forward_id, ] else forward[0, ]
  filtered_paf_path <- file.path(root, "alignments", paste0("human_vs_", query_id, ".paf"))
  if (nrow(kept)) {
    write.table(kept[, names(empty_paf())[1:12]], filtered_paf_path, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    block_table <- data.frame(
      comparison = paste0("human_vs_", query_id),
      human_start = kept$human_start,
      human_end = kept$human_end,
      primate_start = kept$primate_start,
      primate_end = kept$primate_end,
      strand = kept$strand,
      matching_bases = kept$matching_bases,
      block_length = kept$block_length,
      mapq = kept$mapq,
      identity = kept$identity,
      human_match_coverage = kept$target_match_coverage,
      primate_match_coverage = kept$query_match_coverage,
      human_reciprocal_overlap = pairs$human_reciprocal_overlap,
      primate_reciprocal_overlap = pairs$primate_reciprocal_overlap,
      reciprocal_score = pairs$reciprocal_score
    )
    all_blocks[[query_id]] <- block_table
  } else {
    file.create(filtered_paf_path, showWarnings = FALSE)
  }

  all_summary[[query_id]] <- data.frame(
    comparison = paste0("human_vs_", query_id),
    forward_raw_alignments = nrow(forward_all),
    reverse_raw_alignments = nrow(reverse_all),
    forward_quality_candidates = nrow(forward),
    reverse_quality_candidates = nrow(reverse),
    reciprocal_one_to_one_blocks = nrow(kept),
    reciprocal_forward_blocks = sum(kept$strand == "+"),
    reciprocal_inverted_blocks = sum(kept$strand == "-"),
    summed_block_length_bp = sum(kept$block_length),
    min_mapq = min_mapq,
    min_block_bp = min_block_bp,
    min_identity = min_identity,
    min_endpoint_match_coverage = min_endpoint_match_coverage,
    min_reciprocal_overlap = min_reciprocal_overlap
  )
}

if (length(all_blocks)) {
  write.table(do.call(rbind, all_blocks), file.path(root, "reciprocal_unique_blocks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(data.frame(), file.path(root, "reciprocal_unique_blocks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
write.table(do.call(rbind, all_summary), file.path(root, "reciprocal_alignment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
