#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(Rcpp)
})

root <- normalizePath(".")
out_dir <- file.path(root, "marmoset_acro1_fragment_analysis")
seq_dir <- file.path(out_dir, "sequences")
unit_fa <- file.path(root, "promoter_puriofying_analysis_chm13only", "acro_units_structural.fa")

query_files <- sort(list.files(seq_dir, pattern = "\\.fa$", full.names = TRUE))
queries <- DNAStringSet(lapply(query_files, function(path) readDNAStringSet(path)[[1]]))
names(queries) <- sub("\\.fa$", "", basename(query_files))

units <- readDNAStringSet(unit_fa)
unit1 <- units[[1]]

# Representative intervals in structural unit 1. Coordinates below are
# 1-based, inclusive, converted from the actual aligned-fragment coordinates
# in chr21_ACRO_repeatmasker.bed relative to chr21:66240-72574.
component_ranges <- data.frame(
  component = c(
    "L1PA10", "COMP_1624", "MER21B", "COMP_38_main", "ACRO1_short",
    "COMP_38_partial", "ACRO1_long", "COMP_37", "L1MB8", "COMP_1625"
  ),
  class = c(
    "L1PA10", "COMP", "MER21B", "COMP", "ACRO1",
    "COMP", "ACRO1", "COMP", "L1MB8", "COMP"
  ),
  start = c(1, 325, 473, 934, 2252, 2522, 2840, 4438, 5087, 5416),
  end = c(286, 471, 933, 2251, 2515, 2839, 4437, 5085, 5415, 6333),
  stringsAsFactors = FALSE
)

components <- DNAStringSet(lapply(seq_len(nrow(component_ranges)), function(i) {
  subseq(unit1, start = component_ranges$start[i], end = component_ranges$end[i])
}))
names(components) <- component_ranges$component

cppFunction(code = '
Rcpp::List smith_waterman_cpp(std::string query, std::string subject,
                              int match_score = 2, int mismatch_score = -3,
                              int gap_score = -4) {
  const int n = query.size();
  const int m = subject.size();
  std::vector<int> h((n + 1) * (m + 1), 0);
  std::vector<unsigned char> trace((n + 1) * (m + 1), 0);
  int best = 0, best_i = 0, best_j = 0;
  for (int i = 1; i <= n; ++i) {
    for (int j = 1; j <= m; ++j) {
      int diag = h[(i - 1) * (m + 1) + (j - 1)] +
        (query[i - 1] == subject[j - 1] ? match_score : mismatch_score);
      int up = h[(i - 1) * (m + 1) + j] + gap_score;
      int left = h[i * (m + 1) + (j - 1)] + gap_score;
      int val = 0;
      unsigned char direction = 0;
      if (diag > val) { val = diag; direction = 1; }
      if (up > val) { val = up; direction = 2; }
      if (left > val) { val = left; direction = 3; }
      h[i * (m + 1) + j] = val;
      trace[i * (m + 1) + j] = direction;
      if (val > best) { best = val; best_i = i; best_j = j; }
    }
  }
  int i = best_i, j = best_j;
  int columns = 0, aligned_bases = 0, matches = 0;
  while (i > 0 && j > 0 && h[i * (m + 1) + j] > 0) {
    unsigned char direction = trace[i * (m + 1) + j];
    if (direction == 1) {
      ++columns; ++aligned_bases;
      if (query[i - 1] == subject[j - 1]) ++matches;
      --i; --j;
    } else if (direction == 2) {
      ++columns; --i;
    } else if (direction == 3) {
      ++columns; --j;
    } else {
      break;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("score") = best,
    Rcpp::Named("query_start") = i + 1,
    Rcpp::Named("query_end") = best_i,
    Rcpp::Named("subject_start") = j + 1,
    Rcpp::Named("subject_end") = best_j,
    Rcpp::Named("alignment_columns") = columns,
    Rcpp::Named("aligned_bases") = aligned_bases,
    Rcpp::Named("matches") = matches
  );
}
')

alignment_stats <- function(query, subject, strand) {
  aln <- smith_waterman_cpp(as.character(query), as.character(subject))
  aligned_bases <- aln$aligned_bases
  data.frame(
    strand = strand,
    score = aln$score,
    query_start = aln$query_start,
    query_end = aln$query_end,
    subject_start = aln$subject_start,
    subject_end = aln$subject_end,
    alignment_columns = aln$alignment_columns,
    aligned_bases = aligned_bases,
    matches = aln$matches,
    identity = if (aligned_bases > 0) aln$matches / aligned_bases else NA_real_,
    stringsAsFactors = FALSE
  )
}

rows <- list()
k <- 1L
for (q_name in names(queries)) {
  for (i in seq_along(components)) {
    subject <- components[[i]]
    fwd <- alignment_stats(queries[[q_name]], subject, "+")
    rev <- alignment_stats(reverseComplement(queries[[q_name]]), subject, "-")
    best <- if (fwd$score >= rev$score) fwd else rev
    rows[[k]] <- cbind(
      data.frame(
        query = q_name,
        query_length = length(queries[[q_name]]),
        component = names(components)[i],
        class = component_ranges$class[i],
        component_length = length(subject),
        stringsAsFactors = FALSE
      ),
      best
    )
    k <- k + 1L
  }
}

results <- do.call(rbind, rows)
results$query_coverage <- results$aligned_bases / results$query_length
results$component_coverage <- results$aligned_bases / results$component_length
results <- results[order(results$query, -results$score), ]
write.table(
  results,
  file = file.path(out_dir, "component_local_alignments.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

writeXStringSet(components, file.path(out_dir, "human_ACRO_unit1_components.fa"))

top <- do.call(rbind, lapply(split(results, results$query), function(x) head(x, 5)))
write.table(
  top,
  file = file.path(out_dir, "component_local_alignments_top5.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Retain substantial, non-redundant component matches. A shorter same-class
# match nested by >=50% within a higher-scoring match is redundant.
robust <- results[
  results$aligned_bases >= 100 & results$identity >= 0.70,
]
selected_rows <- list()
kk <- 1L
for (q_name in unique(robust$query)) {
  q_rows <- robust[robust$query == q_name, ]
  q_rows <- q_rows[order(-q_rows$score), ]
  kept <- q_rows[0, ]
  for (i in seq_len(nrow(q_rows))) {
    candidate <- q_rows[i, ]
    same_class <- kept[kept$class == candidate$class, ]
    redundant <- FALSE
    if (nrow(same_class) > 0) {
      overlaps <- pmax(
        0,
        pmin(same_class$query_end, candidate$query_end) -
          pmax(same_class$query_start, candidate$query_start) + 1
      )
      candidate_span <- candidate$query_end - candidate$query_start + 1
      redundant <- any(overlaps / candidate_span >= 0.5)
    }
    if (!redundant) kept <- rbind(kept, candidate)
  }
  selected_rows[[kk]] <- kept
  kk <- kk + 1L
}
segments <- do.call(rbind, selected_rows)

segments$query_start_original <- ifelse(
  segments$strand == "+", segments$query_start,
  segments$query_length - segments$query_end + 1
)
segments$query_end_original <- ifelse(
  segments$strand == "+", segments$query_end,
  segments$query_length - segments$query_start + 1
)
name_match <- regexec("^([^_]+)_([0-9]+)-([0-9]+)$", segments$query)
name_parts <- regmatches(segments$query, name_match)
segments$accession <- vapply(name_parts, `[[`, character(1), 2)
segments$locus_start_bed0 <- as.integer(vapply(name_parts, `[[`, character(1), 3))
segments$locus_end_bed0 <- as.integer(vapply(name_parts, `[[`, character(1), 4))
segments$genome_start_bed0 <- segments$locus_start_bed0 + segments$query_start_original - 1
segments$genome_end_bed0 <- segments$locus_start_bed0 + segments$query_end_original
segments$genome_start_1based <- segments$genome_start_bed0 + 1
segments$genome_end_1based <- segments$genome_end_bed0
segments <- segments[order(segments$accession, segments$genome_start_bed0), ]
write.table(
  segments,
  file = file.path(out_dir, "robust_component_segments.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

locus_rows <- lapply(split(results, results$query), function(x) {
  best_acro <- x[x$class == "ACRO1", ][1, ]
  best_non <- x[x$class != "ACRO1", ][1, ]
  substantial_non <- best_non$aligned_bases >= 100 && best_non$identity >= 0.70
  data.frame(
    query = x$query[1],
    query_length = x$query_length[1],
    best_ACRO1_component = best_acro$component,
    ACRO1_aligned_bases = best_acro$aligned_bases,
    ACRO1_identity = best_acro$identity,
    ACRO1_query_coverage = best_acro$query_coverage,
    best_non_ACRO1_component = best_non$component,
    non_ACRO1_aligned_bases = best_non$aligned_bases,
    non_ACRO1_identity = best_non$identity,
    non_ACRO1_query_coverage = best_non$query_coverage,
    interpretation = if (substantial_non) {
      paste("ACRO1 plus", best_non$component)
    } else {
      "ACRO1 only; no substantial non-ACRO1 component match"
    },
    stringsAsFactors = FALSE
  )
})
locus_summary <- do.call(rbind, locus_rows)
write.table(
  locus_summary,
  file = file.path(out_dir, "locus_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(top[, c(
  "query", "query_length", "component", "class", "strand", "score",
  "aligned_bases", "identity", "query_coverage", "component_coverage"
)], row.names = FALSE)
