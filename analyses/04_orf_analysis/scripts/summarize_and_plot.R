args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[[1]] else "primate_acro1_orf_test"

arrays <- read.delim(file.path(outdir, "array_orf_results.tsv"), na.strings = "NA", check.names = FALSE)
tested <- arrays[arrays$status == "ok", ]
tested$q_atg_within_assembly <- ave(
  tested$p_max_atg_open_aa,
  tested$key,
  FUN = function(x) p.adjust(x, method = "BH")
)
tested$q_stopfree_within_assembly <- ave(
  tested$p_max_stopfree_aa,
  tested$key,
  FUN = function(x) p.adjust(x, method = "BH")
)
tested$q_atg_all_arrays <- p.adjust(tested$p_max_atg_open_aa, method = "BH")
tested$q_stopfree_all_arrays <- p.adjust(tested$p_max_stopfree_aa, method = "BH")
write.table(
  tested,
  file.path(outdir, "array_orf_results_with_fdr.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

tracts <- read.delim(file.path(outdir, "tract_orf_results.tsv"), check.names = FALSE)
tracts$atg_fraction_of_tract <- 3 * tracts$actual_max_atg_open_aa / tracts$tract_length
tracts$stopfree_fraction_of_tract <- 3 * tracts$actual_max_stopfree_aa / tracts$tract_length

keys <- unique(tracts$key)
tract_summary <- do.call(rbind, lapply(keys, function(key) {
  z <- tracts[tracts$key == key, ]
  data.frame(
    key = key,
    species = z$species[[1]],
    assembly = z$assembly[[1]],
    de_novo_tracts = nrow(z),
    median_tract_length_bp = median(z$tract_length),
    tracts_with_zero_stop_frame = sum(z$actual_zero_stop_frames > 0),
    fraction_tracts_with_zero_stop_frame = mean(z$actual_zero_stop_frames > 0),
    median_atg_open_aa = median(z$actual_max_atg_open_aa),
    median_atg_fraction_of_tract = median(z$atg_fraction_of_tract),
    median_stopfree_aa = median(z$actual_max_stopfree_aa),
    median_stopfree_fraction_of_tract = median(z$stopfree_fraction_of_tract)
  )
}))
write.table(
  tract_summary,
  file.path(outdir, "actual_tract_species_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

dir.create(file.path(outdir, "figures"), showWarnings = FALSE)
species_summary <- read.delim(file.path(outdir, "species_summary.tsv"), na.strings = "NA", check.names = FALSE)
species_summary <- species_summary[species_summary$status == "ok", ]
order_keys <- species_summary$key
tract_summary <- tract_summary[match(order_keys, tract_summary$key), ]
labels <- c("CHM13", "HG002", "Gorilla", "Chimp", "Bonobo", "Sumatran\norangutan", "Bornean\norangutan")
colors <- c("#3A6EA5", "#73A2C6", "#8B5E3C", "#6A994E", "#386641", "#E09F3E", "#C56A1A")

draw_plot <- function() {
  old <- par(mfrow = c(1, 3), mar = c(7.5, 4.5, 3, 1), las = 1, family = "sans")
  on.exit(par(old))
  barplot(
    tract_summary$fraction_tracts_with_zero_stop_frame,
    names.arg = labels, col = colors, border = NA, ylim = c(0, 1), las = 2,
    ylab = "Fraction of de novo tracts", main = "At least one stop-free frame"
  )
  abline(h = 0.9, lty = 3, col = "grey45")
  barplot(
    tract_summary$median_atg_open_aa,
    names.arg = labels, col = colors, border = NA, las = 2,
    ylab = "Amino acids", main = "Median ATG-initiated open segment"
  )
  midpoint <- barplot(
    rbind(species_summary$mean_observed_atg_open_aa_10x, species_summary$mean_null_atg_open_aa_10x),
    beside = TRUE, names.arg = labels, col = c("#294C60", "#C8C8C8"), border = NA,
    las = 2, ylab = "Amino acids", main = "Array consensus versus dinucleotide null"
  )
  legend("topright", legend = c("Observed", "Null mean"), fill = c("#294C60", "#C8C8C8"), bty = "n", cex = 0.85)
}

png(file.path(outdir, "figures", "de_novo_acro1_orf_summary.png"), width = 2400, height = 950, res = 180)
draw_plot()
dev.off()
pdf(file.path(outdir, "figures", "de_novo_acro1_orf_summary.pdf"), width = 13.3, height = 5.3)
draw_plot()
dev.off()

cat(file.path(outdir, "actual_tract_species_summary.tsv"), "\n")
