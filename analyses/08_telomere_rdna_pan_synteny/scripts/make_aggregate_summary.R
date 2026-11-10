#!/usr/bin/env Rscript

all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", all_args[grep("^--file=", all_args)][1])
root <- normalizePath(dirname(script_path), mustWork = TRUE)
chromosomes <- c("hsa13", "hsa14", "hsa15", "hsa21", "hsa22")

read_with_chr <- function(chromosome, filename) {
  path <- file.path(root, chromosome, filename)
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  cbind(chromosome = chromosome, tab)
}

reciprocal_summary <- do.call(
  rbind,
  lapply(chromosomes, read_with_chr, filename = "reciprocal_alignment_summary.tsv")
)
write.table(
  reciprocal_summary,
  file.path(root, "all_chromosomes_reciprocal_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

display_summary <- do.call(
  rbind,
  lapply(chromosomes, read_with_chr, filename = "alignment_summary.tsv")
)
write.table(
  display_summary,
  file.path(root, "all_chromosomes_display_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

blocks <- do.call(
  rbind,
  lapply(chromosomes, read_with_chr, filename = "reciprocal_unique_blocks.tsv")
)

displayed <- do.call(
  rbind,
  lapply(chromosomes, read_with_chr, filename = "displayed_blocks.tsv")
)

validation <- data.frame(
  check = c(
    "reciprocal_block_mapq_gte_30",
    "reciprocal_block_length_gte_5000",
    "reciprocal_identity_gte_0.95",
    "reciprocal_endpoint_match_coverage_gte_0.80",
    "reciprocal_overlap_gte_0.80",
    "displayed_target_censat_overlap_is_zero",
    "displayed_primate_censat_overlap_is_zero"
  ),
  passed = c(
    all(blocks$mapq >= 30),
    all(blocks$block_length >= 5000),
    all(blocks$identity >= 0.95),
    all(blocks$human_match_coverage >= 0.80 & blocks$primate_match_coverage >= 0.80),
    all(blocks$human_reciprocal_overlap >= 0.80 & blocks$primate_reciprocal_overlap >= 0.80),
    all(displayed$target_censat_overlap_fraction == 0),
    all(displayed$query_censat_overlap_fraction == 0)
  ),
  observed_minimum = c(
    min(blocks$mapq),
    min(blocks$block_length),
    min(blocks$identity),
    min(c(blocks$human_match_coverage, blocks$primate_match_coverage)),
    min(c(blocks$human_reciprocal_overlap, blocks$primate_reciprocal_overlap)),
    max(displayed$target_censat_overlap_fraction),
    max(displayed$query_censat_overlap_fraction)
  ),
  stringsAsFactors = FALSE
)
write.table(
  validation,
  file.path(root, "validation_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (!all(validation$passed)) {
  stop("One or more validation checks failed; inspect validation_summary.tsv")
}

message("Wrote aggregate summaries and all validation checks passed.")
