#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else normalizePath(".")
data_root <- Sys.getenv("ACRO1_DATA_ROOT", file.path(root, "data"))
out_dir <- if (length(args) >= 2) args[[2]] else file.path(root, "acro_telomere_rdna_position_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_bed <- function(path) {
  x <- read.delim(
    path, header = FALSE, sep = "\t", quote = "", comment.char = "",
    fill = TRUE, stringsAsFactors = FALSE
  )
  x <- x[!is.na(x$V2) & !is.na(x$V3), ]
  names(x)[1:6] <- c("chrom", "start", "end", "annotation", "score", "strand")
  x$start <- as.numeric(x$start)
  x$end <- as.numeric(x$end)
  x
}

merge_touching <- function(x) {
  if (nrow(x) == 0) return(x)
  x <- x[order(x$chrom, x$start, x$end), ]
  rows <- list()
  k <- 1L
  for (chrom in unique(x$chrom)) {
    z <- x[x$chrom == chrom, ]
    cur_start <- z$start[1]
    cur_end <- z$end[1]
    source_records <- 1L
    if (nrow(z) > 1) {
      for (i in 2:nrow(z)) {
        if (z$start[i] <= cur_end) {
          cur_end <- max(cur_end, z$end[i])
          source_records <- source_records + 1L
        } else {
          rows[[k]] <- data.frame(
            chrom = chrom, start = cur_start, end = cur_end,
            source_records = source_records
          )
          k <- k + 1L
          cur_start <- z$start[i]
          cur_end <- z$end[i]
          source_records <- 1L
        }
      }
    }
    rows[[k]] <- data.frame(
      chrom = chrom, start = cur_start, end = cur_end,
      source_records = source_records
    )
    k <- k + 1L
  }
  do.call(rbind, rows)
}

human_inputs <- data.frame(
  assembly = c("T2T-CHM13v2.0", "HG002v1.1"),
  path = c(
    file.path(root, "chm13", "chm13v2.0.cenSatv2.1.bed"),
    file.path(data_root, "primate_genomes", "hg002v1.1.cenSatv2.0.bed")
  ),
  pattern = c(
    "^chr(13|14|15|21|22)$",
    "^chr(13|14|15|21|22)_(MATERNAL|PATERNAL)$"
  ),
  stringsAsFactors = FALSE
)

array_rows <- list()
for (i in seq_len(nrow(human_inputs))) {
  d <- read_bed(human_inputs$path[i])
  d <- d[
    grepl(human_inputs$pattern[i], d$chrom) &
      grepl("ACRO1", d$annotation, fixed = TRUE),
  ]
  merged <- merge_touching(d)
  merged$assembly <- human_inputs$assembly[i]
  merged$length_bp <- merged$end - merged$start
  merged$telomere_distance_bp <- merged$start
  array_rows[[i]] <- merged
}
human_arrays <- do.call(rbind, array_rows)
human_arrays <- human_arrays[, c(
  "assembly", "chrom", "start", "end", "length_bp",
  "telomere_distance_bp", "source_records"
)]
write.table(
  human_arrays, file.path(out_dir, "human_acrocentric_ACRO_arrays.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

thresholds <- c(25000, 50000, 75000, 100000, 150000, 200000)
windows <- c(1000000, 2000000, 3000000, 4000000, 5000000)
grid <- expand.grid(
  assembly = unique(human_arrays$assembly),
  large_threshold_bp = thresholds,
  telomere_window_bp = windows,
  stringsAsFactors = FALSE
)
grid$n_large_arrays <- 0L
grid$n_within_window <- 0L
grid$array_fraction_within <- NA_real_
grid$large_array_bp <- 0
grid$bp_within_window <- 0
grid$bp_fraction_within <- NA_real_
for (i in seq_len(nrow(grid))) {
  z <- human_arrays[
    human_arrays$assembly == grid$assembly[i] &
      human_arrays$length_bp >= grid$large_threshold_bp[i],
  ]
  inside <- z$telomere_distance_bp < grid$telomere_window_bp[i]
  grid$n_large_arrays[i] <- nrow(z)
  grid$n_within_window[i] <- sum(inside)
  grid$array_fraction_within[i] <- if (nrow(z)) mean(inside) else NA_real_
  grid$large_array_bp[i] <- sum(z$length_bp)
  grid$bp_within_window[i] <- sum(z$length_bp[inside])
  grid$bp_fraction_within[i] <- if (sum(z$length_bp)) {
    sum(z$length_bp[inside]) / sum(z$length_bp)
  } else NA_real_
}
write.table(
  grid, file.path(out_dir, "human_telomere_threshold_sensitivity.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ACRO and HSat1A burden on either side of the principal rDNA/GAP/rDNA block
# in chimpanzee and bonobo primary acrocentric homologs. Only chromosomes with
# an annotated rDNA block are estimable. Coordinates increase from p telomere
# toward the centromere, so below/above the block are distal/proximal.
pan_inputs <- data.frame(
  species = c("Chimpanzee", "Bonobo"),
  path = c(
    file.path(data_root, "primate_genomes", "PanTro.cenSatv2.0.bed"),
    file.path(data_root, "primate_genomes", "PanPan.cenSatv2.0.bed")
  ),
  stringsAsFactors = FALSE
)

side_rows <- list()
k <- 1L
for (i in seq_len(nrow(pan_inputs))) {
  d <- read_bed(pan_inputs$path[i])
  d <- d[grepl("_hsa(13|14|15|21|22)$", d$chrom), ]
  for (chrom in unique(d$chrom)) {
    z <- d[d$chrom == chrom, ]
    rdna <- z[z$annotation == "rDNA", ]
    if (nrow(rdna) == 0) next
    block_start <- min(rdna$start)
    block_end <- max(rdna$end)
    for (feature in c("ACRO", "HSat1A", "ACRO_plus_HSat1A")) {
      is_feature <- if (feature == "ACRO") {
        grepl("ACRO1", z$annotation, fixed = TRUE)
      } else if (feature == "HSat1A") {
        z$annotation == "HSat1A"
      } else {
        grepl("ACRO1", z$annotation, fixed = TRUE) | z$annotation == "HSat1A"
      }
      q <- z[is_feature, ]
      distal_bp <- sum(pmax(0, pmin(q$end, block_start) - q$start))
      proximal_bp <- sum(pmax(0, q$end - pmax(q$start, block_end)))
      overlap_bp <- sum(pmax(0, pmin(q$end, block_end) - pmax(q$start, block_start)))
      side_rows[[k]] <- data.frame(
        species = pan_inputs$species[i], chrom = chrom, feature = feature,
        rdna_block_start = block_start, rdna_block_end = block_end,
        distal_bp = distal_bp, proximal_bp = proximal_bp,
        overlap_bp = overlap_bp, total_classified_bp = distal_bp + proximal_bp,
        proximal_fraction = if (distal_bp + proximal_bp > 0) {
          proximal_bp / (distal_bp + proximal_bp)
        } else NA_real_,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
}
side_table <- do.call(rbind, side_rows)
write.table(
  side_table, file.path(out_dir, "pan_rdna_side_by_chromosome.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

side_summary <- aggregate(
  cbind(distal_bp, proximal_bp, overlap_bp) ~ species + feature,
  side_table, sum
)
side_summary$total_classified_bp <- side_summary$distal_bp + side_summary$proximal_bp
side_summary$proximal_fraction <- side_summary$proximal_bp / side_summary$total_classified_bp
side_summary$n_rdna_chromosomes <- as.integer(table(side_table$species)[
  side_summary$species
] / length(unique(side_table$feature)))
write.table(
  side_summary, file.path(out_dir, "pan_rdna_side_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Summarize the orientation of the accepted CenSat-masked reciprocal
# unique-sequence blocks in the rDNA-bearing homolog set used for the figures.
synteny_root <- file.path(
  root, "synteny_unique_masked_hsa13_14_15_21_22_all_primates"
)
orientation_rows <- list()
k <- 1L
for (homolog in c("hsa13", "hsa14", "hsa21", "hsa22")) {
  path <- file.path(synteny_root, homolog, "reciprocal_unique_blocks.tsv")
  blocks <- read.delim(path, sep = "\t", stringsAsFactors = FALSE)
  blocks <- blocks[blocks$comparison %in% c(
    "human_vs_chimpanzee", "human_vs_bonobo"
  ), ]
  if (nrow(blocks) == 0) next
  blocks$species <- ifelse(
    blocks$comparison == "human_vs_chimpanzee", "Chimpanzee", "Bonobo"
  )
  for (species in unique(blocks$species)) {
    z <- blocks[blocks$species == species, ]
    orientation_rows[[k]] <- data.frame(
      homolog = homolog,
      species = species,
      forward_blocks = sum(z$strand == "+"),
      inverted_blocks = sum(z$strand == "-"),
      forward_matching_bp = sum(z$matching_bases[z$strand == "+"]),
      inverted_matching_bp = sum(z$matching_bases[z$strand == "-"]),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
orientation <- do.call(rbind, orientation_rows)
orientation_summary <- aggregate(
  cbind(
    forward_blocks, inverted_blocks,
    forward_matching_bp, inverted_matching_bp
  ) ~ species,
  orientation, sum
)
orientation_summary$inverted_matching_fraction <- with(
  orientation_summary,
  inverted_matching_bp / (forward_matching_bp + inverted_matching_bp)
)
orientation_summary$homologs_with_blocks <- as.integer(table(orientation$species)[
  orientation_summary$species
])
write.table(
  orientation,
  file.path(out_dir, "pan_synteny_orientation_by_homolog.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  orientation_summary,
  file.path(out_dir, "pan_synteny_orientation_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("CHM13 >=100 kb sensitivity:\n")
print(grid[
  grid$assembly == "T2T-CHM13v2.0" & grid$large_threshold_bp == 100000,
], row.names = FALSE)
cat("\nPan side summary:\n")
print(side_summary, row.names = FALSE)
cat("\nPan reciprocal-synteny orientation summary:\n")
print(orientation_summary, row.names = FALSE)
