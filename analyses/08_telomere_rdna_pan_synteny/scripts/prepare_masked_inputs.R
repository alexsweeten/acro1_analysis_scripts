#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: prepare_masked_inputs.R SOURCE_ROOT DEST_ROOT")
source_root <- normalizePath(args[[1]], mustWork = TRUE)
dest_root <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(dest_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dest_root, "masks"), recursive = TRUE, showWarnings = FALSE)

config <- read.delim(file.path(source_root, "species_config.tsv"), header = TRUE, sep = "\t",
                     quote = "", na.strings = "NA", stringsAsFactors = FALSE)
source_project_root <- dirname(source_root)
config$fasta_actual <- vapply(config$fasta, function(x) {
  normalizePath(file.path(source_project_root, x), mustWork = TRUE)
}, character(1))
config$bed_actual <- vapply(config$bed, function(x) {
  normalizePath(file.path(source_project_root, x), mustWork = TRUE)
}, character(1))

output_config <- config[, c("species_id", "label", "assembly", "sequence", "fasta", "bed",
                            "original_endpoint_bp", "cluster_start_bp", "cluster_end_bp")]
output_config$fasta <- file.path("..", "..", "ACREAGE", "primate_genomes", basename(config$fasta_actual))
output_config$bed <- file.path("..", "..", "ACREAGE", "primate_genomes", basename(config$bed_actual))
write.table(output_config, file.path(dest_root, "species_config.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

merge_intervals <- function(start, end) {
  if (!length(start)) return(data.frame(start = numeric(), end = numeric()))
  o <- order(start, end)
  start <- start[o]
  end <- end[o]
  out_start <- start[[1]]
  out_end <- end[[1]]
  for (i in seq.int(2, length(start))) {
    last <- length(out_start)
    if (start[[i]] <= out_end[[last]]) {
      out_end[[last]] <- max(out_end[[last]], end[[i]])
    } else {
      out_start <- c(out_start, start[[i]])
      out_end <- c(out_end, end[[i]])
    }
  }
  data.frame(start = out_start, end = out_end)
}

mask_summary <- list()
for (i in seq_len(nrow(config))) {
  row <- config[i, ]
  bed <- read.delim(row$bed_actual, header = FALSE, sep = "\t", quote = "",
                    comment.char = "", skip = 1, fill = TRUE, stringsAsFactors = FALSE)
  bed <- bed[bed[[1]] == row$sequence & bed[[2]] < row$original_endpoint_bp & bed[[4]] != "ct", 1:4]
  names(bed) <- c("sequence", "start", "end", "annotation")
  bed$start <- pmax(0, bed$start)
  bed$end <- pmin(row$original_endpoint_bp, bed$end)
  bed <- bed[bed$end > bed$start, ]
  merged <- merge_intervals(bed$start, bed$end)
  mask_bed <- data.frame(
    sequence = rep(paste0(row$sequence, ":1-", row$original_endpoint_bp), nrow(merged)),
    start = merged$start,
    end = merged$end
  )
  write.table(mask_bed, file.path(dest_root, "masks", paste0(row$species_id, ".censat_mask.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  masked_bp <- sum(merged$end - merged$start)
  mask_summary[[i]] <- data.frame(
    species_id = row$species_id,
    sequence = row$sequence,
    endpoint_bp = row$original_endpoint_bp,
    mask_intervals = nrow(merged),
    masked_bp = masked_bp,
    retained_unmasked_bp = row$original_endpoint_bp - masked_bp,
    retained_unmasked_fraction = (row$original_endpoint_bp - masked_bp) / row$original_endpoint_bp
  )
}
write.table(do.call(rbind, mask_summary), file.path(dest_root, "masking_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
