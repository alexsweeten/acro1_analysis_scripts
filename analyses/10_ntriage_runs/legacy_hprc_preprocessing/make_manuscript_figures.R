#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
tables <- file.path(root, "tables")
figures <- file.path(root, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

cluster_levels <- paste0("C", 1:8)
cluster_palette <- c(
  C1 = "#4477AA", C2 = "#EE6677", C3 = "#228833", C4 = "#CCBB44",
  C5 = "#66CCEE", C6 = "#AA3377", C7 = "#BBBBBB", C8 = "#000000"
)
status_palette <- c(`No recombination` = "#6A6A6A", Recombination = "#D55E00")

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(color = "#444444", size = rel(0.86)),
      plot.caption = element_text(color = "#555555", size = rel(0.72), hjust = 0),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.background = element_blank(), strip.text = element_text(face = "bold")
    )
}

save_both <- function(plot, stem, width, height) {
  ggsave(file.path(figures, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 400, bg = "white")
  ggsave(file.path(figures, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = "pdf", useDingbats = FALSE, bg = "white")
}

# Main figure: sequence-space clusters, paired recombination effects, and ACRO1 classes.
pca <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(tables, "unit_pca_metadata.tsv.gz"))))
pca_var <- fread(file.path(tables, "pca_explained_variance.tsv"))
set.seed(20261026)
pca[, cluster_label := factor(paste0("C", cluster), levels = cluster_levels)]
pca_plot_data <- pca[, .SD[sample.int(.N, min(.N, 1400L))], by = cluster_label]
pca_centroid <- pca[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = cluster_label]
pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * pca_var[principal_component == 1, explained_variance_ratio])
pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * pca_var[principal_component == 2, explained_variance_ratio])

p_pca <- ggplot(pca_plot_data, aes(PC1, PC2, color = cluster_label)) +
  geom_point(size = 0.38, alpha = 0.24, stroke = 0) +
  geom_text_repel(data = pca_centroid, aes(label = cluster_label), color = "black", fontface = "bold",
                  size = 3.2, seed = 20261026, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = cluster_palette, drop = FALSE) +
  labs(title = "N-TRIAGE sequence space", subtitle = "Array-balanced visualization; colors are native k = 8 clusters",
       x = pc1_lab, y = pc2_lab, color = "Cluster") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme_manuscript()

effects <- fread(file.path(tables, "recombination_cluster_fraction_tests.tsv"))
effects[, cluster_label := factor(paste0("C", cluster), levels = rev(cluster_levels))]
effects[, `:=`(
  effect_pct = 100 * mean_within_person_fraction_difference,
  lo_pct = 100 * bootstrap_ci95_lo,
  hi_pct = 100 * bootstrap_ci95_hi,
  sig = fifelse(q_bh_across_8_clusters < 0.05, "BH q < 0.05", "not significant")
)]
p_effect <- ggplot(effects, aes(effect_pct, cluster_label, color = paste0("C", cluster))) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo_pct, xmax = hi_pct), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(values = cluster_palette, guide = "none") +
  labs(title = "Matched recombination association",
       subtitle = "Recombination minus no recombination; 230 individuals",
       x = "Mean within-person difference (percentage points)", y = NULL) +
  theme_manuscript()

acro <- fread(file.path(tables, "recombination_acro1_subfamily_tests.tsv"))
acro_long <- rbindlist(list(
  acro[, .(subfamily, status = "No recombination", fraction = mean_fraction_no_recombination)],
  acro[, .(subfamily, status = "Recombination", fraction = mean_fraction_recombination)]
))
acro_long[, status := factor(status, levels = c("No recombination", "Recombination"))]
p_acro <- ggplot(acro_long, aes(subfamily, fraction, fill = status)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = status_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "ACRO1 structural composition",
       subtitle = "Bars: array means; all three matched differences have BH q = 1e-5",
       x = "N-TRIAGE ACRO1 class", y = "Mean fraction of classified monomers", fill = NULL) +
  theme_manuscript() + theme(legend.position = "top")

main_plot <- p_pca + (p_effect / p_acro) +
  plot_layout(widths = c(1.18, 1)) +
  plot_annotation(
    title = "Core-only N-TRIAGE analysis of 2,794 telomeric HPRC ACRO arrays",
    subtitle = "77,276 complete ACRO units; extraction flanks excluded",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "#444444"))
  )
save_both(main_plot, "Figure_telomeric_NTRIAGE_recombination", 13.2, 7.8)

# Population figure: raw individual means plus covariate-adjusted omnibus tests.
pop <- fread(file.path(tables, "superpopulation_individual_means.tsv"))
pop_order <- c("AFR", "AMR", "EAS", "EUR", "SAS")
pop[, pop_label := factor(sprintf("%s (n=%d)", superpopulation, individuals),
                          levels = sprintf("%s (n=%d)", pop_order,
                                           individuals[match(pop_order, superpopulation)]))]
pop_long <- melt(pop, id.vars = c("superpopulation", "pop_label", "individuals"),
                 measure.vars = paste0("cluster_", 1:8, "_fraction"),
                 variable.name = "cluster", value.name = "fraction")
pop_long[, cluster_label := factor(sub("cluster_([0-9]+)_fraction", "C\\1", cluster),
                                   levels = cluster_levels)]
pop_long[, z_within_cluster := as.numeric(scale(fraction)), by = cluster_label]

p_heat <- ggplot(pop_long, aes(cluster_label, pop_label, fill = z_within_cluster)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = percent(fraction, accuracy = 0.1)), size = 2.6) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
                       name = "Within-cluster\nZ score") +
  labs(title = "Cluster composition by superpopulation",
       subtitle = "Cells show unadjusted individual-level means",
       x = "N-TRIAGE cluster", y = NULL) +
  theme_manuscript() + theme(panel.grid = element_blank())

pop_tests <- fread(file.path(tables, "superpopulation_cluster_tests.tsv"))
pop_global <- fread(file.path(tables, "superpopulation_global_composition_test.tsv"))
pop_tests[, cluster_label := factor(paste0("C", cluster), levels = cluster_levels)]
pop_tests[, significant := q_bh_across_8_clusters < 0.05]
p_poptest <- ggplot(pop_tests, aes(cluster_label, partial_R2, fill = cluster_label)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = ifelse(significant, "*", "")), vjust = -0.25, size = 5) +
  scale_fill_manual(values = cluster_palette, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Covariate-adjusted population effect",
       subtitle = "Freedman-Lane tests control recombination fraction and mean log length; * BH q < 0.05",
       x = "N-TRIAGE cluster", y = "Partial R2") +
  theme_manuscript()

population_plot <- p_heat / p_poptest +
  plot_layout(heights = c(1.2, 0.85)) +
  plot_annotation(
    title = "Exploratory population structure in telomeric ACRO arrays",
    subtitle = sprintf("Global adjusted composition test: partial R2 = %.1f%%, permutation P = %.2g (%d individuals)",
                       100 * pop_global$partial_R2[[1]], pop_global$p_value[[1]],
                       pop_global$individuals[[1]]),
    caption = "Population metadata were unavailable for 78 of 231 individuals; AMR (n = 6) and EUR (n = 8) estimates are imprecise.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10, color = "#444444"),
                  plot.caption = element_text(size = 8, color = "#555555", hjust = 0))
  )
save_both(population_plot, "Figure_telomeric_NTRIAGE_population", 10.8, 8.3)

# Supplement: array-level metrics, k sensitivity, and PCA by recombination label.
metrics <- fread(file.path(tables, "recombination_array_metric_tests.tsv"))
arrays <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(tables, "array_level_summary.tsv.gz"))))
metric_labels <- c(
  switches_per_100_units = "Switches per 100 units", switch_count = "Switch count",
  complete_acres = "Complete ACRO units", input_length = "Core length",
  estimated_period = "Estimated period", periodicity_score = "Periodicity score",
  cluster_entropy = "Cluster entropy"
)
metric_sd <- vapply(names(metric_labels), function(x) sd(arrays[[x]], na.rm = TRUE), numeric(1))
metrics[, sd_all := metric_sd[metric]]
metrics[, `:=`(effect_sd = mean_within_person_difference / sd_all,
               lo_sd = bootstrap_ci95_lo / sd_all, hi_sd = bootstrap_ci95_hi / sd_all,
               label = factor(metric_labels[metric], levels = rev(unname(metric_labels))))]
p_metrics <- ggplot(metrics, aes(effect_sd, label)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo_sd, xmax = hi_sd), height = 0.18, linewidth = 0.7, color = "#555555") +
  geom_point(size = 2.5, color = "#D55E00") +
  labs(title = "Matched array-property differences",
       subtitle = "Recombination minus no recombination; effects scaled by the cohort SD",
       x = "Mean within-person difference (SD units)", y = NULL) +
  theme_manuscript()

kdat <- fread(file.path(tables, "k_sensitivity_array_balanced_sample.tsv"))
p_k <- ggplot(kdat, aes(k, silhouette)) +
  geom_line(linewidth = 0.7, color = "#444444") + geom_point(size = 2, color = "#4477AA") +
  geom_point(data = kdat[k == 8], size = 3.2, shape = 21, stroke = 1, fill = "#EE6677") +
  geom_hline(yintercept = unique(kdat$native_k8_silhouette_on_same_sample), linetype = 2,
             color = "#AA3377", linewidth = 0.6) +
  scale_x_continuous(breaks = 2:14) +
  labs(title = "Cluster-number sensitivity", subtitle = "k = 3 has the highest refit silhouette; native k = 8 shown by dashed line",
       x = "Refit k", y = "Silhouette") +
  theme_manuscript()

set.seed(20261026)
pca_recomb <- pca[, .SD[sample.int(.N, min(.N, 5000L))], by = has_recombination]
pca_recomb[, status := factor(ifelse(has_recombination == 1, "Recombination", "No recombination"),
                              levels = c("No recombination", "Recombination"))]
p_recomb_pca <- ggplot(pca_recomb, aes(PC1, PC2, color = status)) +
  geom_point(size = 0.35, alpha = 0.18, stroke = 0) +
  scale_color_manual(values = status_palette) +
  labs(title = "Sequence-space separation by array label", x = pc1_lab, y = pc2_lab, color = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme_manuscript() + theme(legend.position = "top")

supp_plot <- (p_metrics | p_k) / p_recomb_pca +
  plot_layout(heights = c(1, 1.08)) +
  plot_annotation(
    title = "Sensitivity and secondary patterns in the core-only N-TRIAGE analysis",
    caption = "Recombination status is an array-level association label; separation should not be interpreted as evidence that recombination caused the sequence differences.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.caption = element_text(size = 8, color = "#555555", hjust = 0))
  )
save_both(supp_plot, "Supplementary_Figure_telomeric_NTRIAGE_sensitivity", 12.5, 9.2)

# Supplement: interpretable repeat and ACRO1 signatures of the eight clusters.
bio <- fread(file.path(tables, "core_cluster_biological_summary.tsv"))
bio[, cluster_label := factor(cluster_label, levels = cluster_levels)]
components <- fread(file.path(tables, "core_cluster_top_components.tsv"))
components[, component_group := fifelse(top_component == "ACRO1", "ACRO1",
                                 fifelse(top_component == "COMP-subunit_ACRO_rnd-5_family-38",
                                         "family-38", "other"))]
component_plot <- components[, .(fraction = sum(within_cluster_fraction)),
                             by = .(cluster_label, component_group)]
component_plot[, cluster_label := factor(cluster_label, levels = cluster_levels)]
component_plot[, component_group := factor(component_group, levels = c("ACRO1", "family-38", "other"))]
p_components <- ggplot(component_plot, aes(cluster_label, fraction, fill = component_group)) +
  geom_col(width = 0.76, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(ACRO1 = "#D55E00", `family-38` = "#0072B2", other = "#999999")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Highest-coverage component", subtitle = "Per-unit top label, not physical coverage fraction",
       x = "N-TRIAGE cluster", y = "Units", fill = NULL) +
  theme_manuscript() + theme(legend.position = "top")

acro_bio <- melt(bio,
                 id.vars = "cluster_label",
                 measure.vars = c("ACRO1_S1_canonical_fraction",
                                  "ACRO1_S2_internal_CCC_deletion_fraction",
                                  "ACRO1_S3_terminal_deletion_fraction"),
                 variable.name = "allele", value.name = "fraction")
acro_bio[, allele := factor(allele,
  levels = c("ACRO1_S1_canonical_fraction", "ACRO1_S2_internal_CCC_deletion_fraction",
             "ACRO1_S3_terminal_deletion_fraction"),
  labels = c("S1 canonical", "S2 internal CCC deletion", "S3 terminal deletion"))]
p_alleles <- ggplot(acro_bio, aes(cluster_label, fraction, fill = allele)) +
  geom_col(width = 0.76, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("S1 canonical" = "#009E73", "S2 internal CCC deletion" = "#E69F00",
                               "S3 terminal deletion" = "#CC79A7")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "ACRO1 allele composition", x = "N-TRIAGE cluster", y = "Classified monomers", fill = NULL) +
  theme_manuscript() + theme(legend.position = "top")

diagnostic <- fread(file.path(tables, "core_cluster_all_diagnostic_feature_contrasts.tsv"))
selected_features <- c("ACRO1", "MER21A", "MER21B", "family-1624", "family-1625", "family-37",
                       "family-38", "L1PA10_3end", "L1PA12_5end", "L1PA14_3end", "L1PA2_3end")
diagnostic <- diagnostic[feature %in% selected_features]
diagnostic[, feature := factor(feature, levels = rev(selected_features))]
diagnostic[, cluster_label := factor(cluster_label, levels = cluster_levels)]
p_diag <- ggplot(diagnostic, aes(cluster_label, feature, fill = standardized_difference)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
                       limits = c(-4.1, 6.2), oob = squish, name = "Standardized\ndifference") +
  labs(title = "Annotated diagnostic k-mer signatures",
       subtitle = "Cluster mean minus all other units; descriptive because these features helped define clusters",
       x = "N-TRIAGE cluster", y = NULL) +
  theme_manuscript() + theme(panel.grid = element_blank())

position <- melt(bio,
                 id.vars = "cluster_label",
                 measure.vars = c("first_position_enrichment_vs_overall", "last_position_enrichment_vs_overall"),
                 variable.name = "position", value.name = "enrichment")
position[, position := factor(position,
  levels = c("first_position_enrichment_vs_overall", "last_position_enrichment_vs_overall"),
  labels = c("First unit", "Last unit"))]
position[, log2_enrichment := log2(enrichment)]
p_position <- ggplot(position, aes(cluster_label, log2_enrichment, fill = position)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "#666666") +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  scale_fill_manual(values = c("First unit" = "#56B4E9", "Last unit" = "#E69F00")) +
  labs(title = "Array-end enrichment", subtitle = "First/last refer to input FASTA order; genomic orientation was not validated",
       x = "N-TRIAGE cluster", y = "log2 enrichment vs all units", fill = NULL) +
  theme_manuscript() + theme(legend.position = "top")

biology_plot <- (p_components | p_alleles) / (p_diag | p_position) +
  plot_annotation(
    title = "Biological signatures of core-only N-TRIAGE clusters",
    caption = "Diagnostic-frequency values are diagnostic k-mer fractions, not genomic coverage. Cluster signatures are descriptive and not independent validation.",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.caption = element_text(size = 8, color = "#555555", hjust = 0))
  )
save_both(biology_plot, "Supplementary_Figure_telomeric_NTRIAGE_cluster_biology", 13.0, 9.0)

cat(sprintf("Wrote manuscript figures to %s\n", figures))
