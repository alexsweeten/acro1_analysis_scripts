#!/usr/bin/env Rscript

# Quantitative ACRO-HSAT1A co-enrichment analysis across HG002, three African
# apes, and two orangutan outgroups. Species-level phylogenetic tests, a
# chromosome-homolog interaction test, and local adjacency tests are separated.

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork = FALSE)

usage <- function(status = 0L) {
  cat(paste0(
    "Usage:\n",
    "  Rscript run_six_species_analysis.R [options]\n\n",
    "Options:\n",
    "  --input-dir DIR   Directory containing the six CenSat BED and FASTA .fai files\n",
    "                    (default: $ACRO_HSAT1A_INPUT_DIR or SCRIPT_DIR/data)\n",
    "  --output-dir DIR  Output directory\n",
    "                    (default: $ACRO_HSAT1A_OUTPUT_DIR or SCRIPT_DIR/results)\n",
    "  --permutations N  Number of run-label permutations (default: 100000)\n",
    "  --seed N          Random seed (default: 20261021)\n",
    "  -h, --help        Show this help message\n\n",
    "For compatibility, four positional arguments are also accepted in this order:\n",
    "INPUT_DIR OUTPUT_DIR PERMUTATIONS SEED.\n"
  ))
  quit(save = "no", status = status)
}

env_or <- function(name, fallback) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else fallback
}

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) usage(0L)

input_dir <- env_or("ACRO_HSAT1A_INPUT_DIR", file.path(script_dir, "data"))
output_dir <- env_or("ACRO_HSAT1A_OUTPUT_DIR", file.path(script_dir, "results"))
n_perm <- as.integer(env_or("ACRO_HSAT1A_PERMUTATIONS", "100000"))
seed <- as.integer(env_or("ACRO_HSAT1A_SEED", "20261021"))

if (length(args) && !any(startsWith(args, "--"))) {
  if (length(args) > 4L) usage(2L)
  if (length(args) >= 1L) input_dir <- args[[1]]
  if (length(args) >= 2L) output_dir <- args[[2]]
  if (length(args) >= 3L) n_perm <- as.integer(args[[3]])
  if (length(args) >= 4L) seed <- as.integer(args[[4]])
} else if (length(args)) {
  i <- 1L
  while (i <= length(args)) {
    option <- args[[i]]
    if (grepl("=", option, fixed = TRUE)) {
      pieces <- strsplit(option, "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
    } else {
      key <- option
      i <- i + 1L
      if (i > length(args)) {
        cat("Missing value for ", key, "\n", sep = "", file = stderr())
        usage(2L)
      }
      value <- args[[i]]
    }
    if (key == "--input-dir") {
      input_dir <- value
    } else if (key == "--output-dir") {
      output_dir <- value
    } else if (key == "--permutations") {
      n_perm <- as.integer(value)
    } else if (key == "--seed") {
      seed <- as.integer(value)
    } else {
      cat("Unknown option: ", key, "\n", sep = "", file = stderr())
      usage(2L)
    }
    i <- i + 1L
  }
}

if (is.na(n_perm) || n_perm < 1L) stop("--permutations must be a positive integer")
if (is.na(seed)) stop("--seed must be an integer")
input_dir <- path.expand(input_dir)
output_dir <- path.expand(output_dir)

for (pkg in c("ggplot2", "ape", "nlme")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package missing: ", pkg)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

metadata <- data.frame(
  tip = c(
    "Human", "Gorilla", "Chimpanzee", "Bonobo",
    "Sumatran_orangutan", "Bornean_orangutan"
  ),
  species = c(
    "Human (HG002)", "Gorilla", "Chimpanzee", "Bonobo",
    "Sumatran orangutan", "Bornean orangutan"
  ),
  assembly = c("HG002v1.1", "mGorGor1", "mPanTro3", "mPanPan1", "mPonAbe1", "mPonPyg2"),
  bed_name = c(
    "hg002v1.1.cenSatv2.0.bed", "GorGor.cenSatv2.0.bed",
    "PanTro.cenSatv2.0.bed", "PanPan.cenSatv2.0.bed",
    "PonAbe.cenSatv2.0.bed", "PonPyg.cenSatv2.0.bed"
  ),
  fai_name = c(
    "hg002v1.1.fasta.gz.fai", "mGorGor1.dip.cur.20231122.fasta.gz.fai",
    "mPanTro3.dip.cur.20231122.fasta.gz.fai",
    "mPanPan1.dip.cur.20231122.fasta.gz.fai",
    "mPonAbe1.dip.cur.20231205.fasta.gz.fai",
    "mPonPyg2.dip.cur.20231122.fasta.gz.fai"
  ),
  primary_pattern = c(
    "^chr(13|14|15|21|22)_(MATERNAL|PATERNAL)$",
    rep("_hsa(13|14|15|21|22)$", 5)
  ),
  broad_clade = c(
    "Human", "African ape", "African ape", "African ape",
    "Orangutan", "Orangutan"
  ),
  stringsAsFactors = FALSE
)
metadata$bed_path <- file.path(input_dir, metadata$bed_name)
metadata$fai_path <- file.path(input_dir, metadata$fai_name)
missing_inputs <- c(
  metadata$bed_path[!file.exists(metadata$bed_path)],
  metadata$fai_path[!file.exists(metadata$fai_path)]
)
if (length(missing_inputs) > 0) stop("Missing input(s): ", paste(missing_inputs, collapse = ", "))

read_bed <- function(path) {
  x <- read.delim(
    path, header = FALSE, sep = "\t", quote = "", comment.char = "",
    fill = TRUE, stringsAsFactors = FALSE
  )
  x <- x[!is.na(x$V2) & !is.na(x$V3), , drop = FALSE]
  names(x)[1:4] <- c("chrom", "start", "end", "annotation")
  x$start <- as.numeric(x$start)
  x$end <- as.numeric(x$end)
  x$bp <- x$end - x$start
  x$category <- ifelse(
    grepl("ACRO1", x$annotation, fixed = TRUE), "ACRO",
    ifelse(x$annotation == "HSat1A", "HSAT1A", "OTHER")
  )
  x
}

read_fai <- function(path) {
  x <- read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(chrom = x$V1, sequence_length = as.numeric(x$V2), stringsAsFactors = FALSE)
}

bed_by_tip <- setNames(vector("list", nrow(metadata)), metadata$tip)
fai_by_tip <- setNames(vector("list", nrow(metadata)), metadata$tip)
species_rows <- list()
chromosome_rows <- list()

for (i in seq_len(nrow(metadata))) {
  m <- metadata[i, ]
  d <- read_bed(m$bed_path)
  f <- read_fai(m$fai_path)
  bed_by_tip[[m$tip]] <- d
  fai_by_tip[[m$tip]] <- f
  primary_contigs <- f$chrom[grepl(m$primary_pattern, f$chrom)]
  d_primary <- d[d$chrom %in% primary_contigs, , drop = FALSE]

  species_rows[[i]] <- data.frame(
    tip = m$tip,
    species = m$species,
    assembly = m$assembly,
    broad_clade = m$broad_clade,
    assembly_bp = sum(f$sequence_length),
    genomewide_acro_bp = sum(d$bp[d$category == "ACRO"]),
    genomewide_hsat1a_bp = sum(d$bp[d$category == "HSAT1A"]),
    primary_acrocentric_acro_bp = sum(d_primary$bp[d_primary$category == "ACRO"]),
    primary_acrocentric_hsat1a_bp = sum(d_primary$bp[d_primary$category == "HSAT1A"]),
    primary_chromosome_haplotypes = length(primary_contigs),
    stringsAsFactors = FALSE
  )

  for (chrom in sort(primary_contigs)) {
    z <- d[d$chrom == chrom, , drop = FALSE]
    homolog <- if (m$tip == "Human") {
      paste0("hsa", sub("chr([0-9]+)_.*", "\\1", chrom))
    } else {
      sub(".*_hsa", "hsa", chrom)
    }
    haplotype <- if (grepl("MATERNAL|_mat_|_hap1_", chrom)) "haplotype_1" else "haplotype_2"
    chromosome_rows[[length(chromosome_rows) + 1L]] <- data.frame(
      tip = m$tip,
      species = m$species,
      broad_clade = m$broad_clade,
      chromosome = chrom,
      homolog = homolog,
      haplotype = haplotype,
      sequence_length = f$sequence_length[match(chrom, f$chrom)],
      acro_bp = sum(z$bp[z$category == "ACRO"]),
      hsat1a_bp = sum(z$bp[z$category == "HSAT1A"]),
      stringsAsFactors = FALSE
    )
  }
}

species <- do.call(rbind, species_rows)
chromosomes <- do.call(rbind, chromosome_rows)
species$genomewide_acro_per_gb <- species$genomewide_acro_bp / species$assembly_bp * 1e9
species$genomewide_hsat1a_per_gb <- species$genomewide_hsat1a_bp / species$assembly_bp * 1e9
species$log_primary_acro <- log10(species$primary_acrocentric_acro_bp + 1)
species$log_primary_hsat1a <- log10(species$primary_acrocentric_hsat1a_bp + 1)
species$log_genome_acro <- log10(species$genomewide_acro_bp + 1)
species$log_genome_hsat1a <- log10(species$genomewide_hsat1a_bp + 1)

write.table(
  species, file.path(output_dir, "six_species_species_burdens.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  chromosomes, file.path(output_dir, "six_species_chromosome_haplotype_burdens.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Enumerate all permutations of n indices (n <= 6 here).
all_permutations <- function(n) {
  x <- as.matrix(expand.grid(rep(list(seq_len(n)), n)))
  x[apply(x, 1, function(z) length(unique(z)) == n), , drop = FALSE]
}

exact_spearman <- function(x, y) {
  permutations <- all_permutations(length(x))
  observed <- cor(rank(x), rank(y))
  null <- apply(permutations, 1, function(idx) cor(rank(x), rank(y[idx])))
  c(
    rho = observed,
    p_one_sided = mean(null >= observed - 1e-15),
    p_two_sided = mean(abs(null) >= abs(observed) - 1e-15),
    permutations = nrow(permutations)
  )
}

primary_rank <- exact_spearman(
  species$primary_acrocentric_acro_bp,
  species$primary_acrocentric_hsat1a_bp
)
genome_rank <- exact_spearman(species$genomewide_acro_bp, species$genomewide_hsat1a_bp)
positive <- species$primary_acrocentric_hsat1a_bp > 0
positive_rank <- exact_spearman(
  species$primary_acrocentric_acro_bp[positive],
  species$primary_acrocentric_hsat1a_bp[positive]
)

# Topology-based PGLS. Branch lengths are generated by Grafen's method because
# no external divergence-time tree was supplied. A unit-edge tree is included
# as a branch-length sensitivity analysis.
library(ape)
library(nlme)
topology <- read.tree(text =
  "((Sumatran_orangutan,Bornean_orangutan),(Gorilla,(Human,(Chimpanzee,Bonobo))));"
)

fit_pgls <- function(tree, label) {
  d <- species[match(tree$tip.label, species$tip), , drop = FALSE]
  fit <- gls(
    log_primary_hsat1a ~ log_primary_acro,
    data = d,
    correlation = corBrownian(1, tree, form = ~ tip),
    method = "ML"
  )
  co <- summary(fit)$tTable["log_primary_acro", ]
  data.frame(
    test = paste0("PGLS Brownian, ", label),
    scope = "all six species",
    statistic_name = "log-log slope",
    statistic = unname(co[["Value"]]),
    standard_error = unname(co[["Std.Error"]]),
    p_one_sided = ifelse(co[["Value"]] >= 0, co[["p-value"]] / 2, 1 - co[["p-value"]] / 2),
    p_two_sided = unname(co[["p-value"]]),
    n = nrow(d),
    stringsAsFactors = FALSE
  )
}

tree_grafen <- compute.brlen(topology, method = "Grafen")
tree_unit <- topology
tree_unit$edge.length <- rep(1, nrow(tree_unit$edge))
pgls_results <- rbind(
  fit_pgls(tree_grafen, "Grafen branch lengths"),
  fit_pgls(tree_unit, "unit branch lengths")
)

# Positive-HSAT1A clade sensitivity (Human + African apes).
positive_topology <- compute.brlen(
  read.tree(text = "(Gorilla,(Human,(Chimpanzee,Bonobo)));"), method = "Grafen"
)
d_positive <- species[match(positive_topology$tip.label, species$tip), , drop = FALSE]
positive_pgls <- gls(
  log_primary_hsat1a ~ log_primary_acro,
  data = d_positive,
  correlation = corBrownian(1, positive_topology, form = ~ tip),
  method = "ML"
)
positive_co <- summary(positive_pgls)$tTable["log_primary_acro", ]

species_test_results <- rbind(
  data.frame(
    test = "Exact tip permutation Spearman",
    scope = "all six species; primary acrocentrics",
    statistic_name = "rho", statistic = primary_rank[["rho"]],
    standard_error = NA_real_, p_one_sided = primary_rank[["p_one_sided"]],
    p_two_sided = primary_rank[["p_two_sided"]], n = 6
  ),
  data.frame(
    test = "Exact tip permutation Spearman",
    scope = "all six species; genome-wide",
    statistic_name = "rho", statistic = genome_rank[["rho"]],
    standard_error = NA_real_, p_one_sided = genome_rank[["p_one_sided"]],
    p_two_sided = genome_rank[["p_two_sided"]], n = 6
  ),
  data.frame(
    test = "Exact tip permutation Spearman",
    scope = "HSAT1A-positive species only; primary acrocentrics",
    statistic_name = "rho", statistic = positive_rank[["rho"]],
    standard_error = NA_real_, p_one_sided = positive_rank[["p_one_sided"]],
    p_two_sided = positive_rank[["p_two_sided"]], n = sum(positive)
  ),
  pgls_results,
  data.frame(
    test = "PGLS Brownian, Grafen branch lengths",
    scope = "HSAT1A-positive species only",
    statistic_name = "log-log slope", statistic = unname(positive_co[["Value"]]),
    standard_error = unname(positive_co[["Std.Error"]]),
    p_one_sided = unname(positive_co[["p-value"]]) / 2,
    p_two_sided = unname(positive_co[["p-value"]]), n = sum(positive)
  )
)
species_test_results$interpretation <- c(
  "Non-phylogenetic species permutation",
  "Non-phylogenetic genome-wide sensitivity",
  "Non-phylogenetic positive-clade sensitivity",
  "Primary phylogenetic sensitivity",
  "Branch-length sensitivity",
  "Positive-clade phylogenetic sensitivity"
)
write.table(
  species_test_results,
  file.path(output_dir, "species_level_association_tests.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Aggregate the two haplotypes for each species x homolog. This tests whether
# species-specific deviations among the five chromosome homologs co-vary after
# removing additive species and homolog effects.
homolog <- aggregate(cbind(acro_bp, hsat1a_bp) ~ tip + species + homolog, chromosomes, sum)
homolog$log_acro <- log10(homolog$acro_bp + 1)
homolog$log_hsat1a <- log10(homolog$hsat1a_bp + 1)
Z <- model.matrix(~ tip + homolog, data = homolog)
project_residual <- function(v, design = Z) as.numeric(v - design %*% qr.coef(qr(design), v))
x_resid <- project_residual(homolog$log_acro)
y_resid <- project_residual(homolog$log_hsat1a)
homolog_beta <- sum(x_resid * y_resid) / sum(x_resid^2)
reduced_fitted <- homolog$log_hsat1a - y_resid
species_indices <- split(seq_len(nrow(homolog)), homolog$tip)
homolog_null <- numeric(n_perm)
for (b in seq_len(n_perm)) {
  permuted_residual <- y_resid
  for (ii in species_indices) permuted_residual[ii] <- sample(permuted_residual[ii])
  y_perm <- project_residual(reduced_fitted + permuted_residual)
  homolog_null[[b]] <- sum(x_resid * y_perm) / sum(x_resid^2)
}
homolog_p_one <- (sum(homolog_null >= homolog_beta) + 1) / (n_perm + 1)
homolog_p_two <- (sum(abs(homolog_null) >= abs(homolog_beta)) + 1) / (n_perm + 1)
write.table(
  homolog,
  file.path(output_dir, "species_by_homolog_burdens.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Exact haplotype-pair sign flip in HSAT1A-positive species. This asks whether
# the haplotype with more ACRO also tends to have more HSAT1A within a matched
# species x chromosome homolog.
positive_chromosomes <- chromosomes[chromosomes$tip %in% species$tip[positive], ]
pair_rows <- lapply(
  split(positive_chromosomes, interaction(positive_chromosomes$tip, positive_chromosomes$homolog)),
  function(z) {
    z <- z[order(z$haplotype), ]
    data.frame(
      tip = z$tip[[1]], species = z$species[[1]], homolog = z$homolog[[1]],
      delta_log_acro = diff(log10(z$acro_bp + 1)),
      delta_log_hsat1a = diff(log10(z$hsat1a_bp + 1)),
      stringsAsFactors = FALSE
    )
  }
)
pairs <- do.call(rbind, pair_rows)
pair_observed <- mean(pairs$delta_log_acro * pairs$delta_log_hsat1a)
n_pairs <- nrow(pairs)
pair_null <- numeric(2^n_pairs)
for (mask in 0:(2^n_pairs - 1L)) {
  signs <- ifelse(as.integer(intToBits(mask))[seq_len(n_pairs)] == 1L, 1, -1)
  pair_null[[mask + 1L]] <- mean(pairs$delta_log_acro * pairs$delta_log_hsat1a * signs)
}
pair_p_one <- mean(pair_null >= pair_observed - 1e-15)
pair_p_two <- mean(abs(pair_null) >= abs(pair_observed) - 1e-15)
write.table(
  pairs,
  file.path(output_dir, "matched_haplotype_differences_positive_species.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

within_genome_tests <- data.frame(
  test = c(
    "Species- and homolog-adjusted Freedman-Lane permutation",
    "Exact matched-haplotype sign flip in HSAT1A-positive species"
  ),
  statistic_name = c("partial log-log slope", "mean product of paired log differences"),
  statistic = c(homolog_beta, pair_observed),
  p_one_sided = c(homolog_p_one, pair_p_one),
  p_two_sided = c(homolog_p_two, pair_p_two),
  n = c(nrow(homolog), n_pairs),
  interpretation = c(
    "Two haplotypes summed within each species x homolog cell",
    paste0(sum(sign(pairs$delta_log_acro) == sign(pairs$delta_log_hsat1a)),
      "/", n_pairs, " matched pairs change in the same direction")
  ),
  stringsAsFactors = FALSE
)
write.table(
  within_genome_tests,
  file.path(output_dir, "within_genome_association_tests.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Spatial test: collapse directly touching runs of the same broad category,
# then permute run labels within chromosome. Orangutans have no HSAT1A and are
# recorded as not estimable. This test is distinct from abundance association.
collapse_runs <- function(z) {
  z <- z[order(z$start, z$end), , drop = FALSE]
  current <- z[1, , drop = FALSE]
  out <- list()
  if (nrow(z) >= 2L) {
    for (i in 2:nrow(z)) {
      if (z$category[[i]] == current$category[[1]] && z$start[[i]] == current$end[[1]]) {
        current$end[[1]] <- z$end[[i]]
      } else {
        out[[length(out) + 1L]] <- current
        current <- z[i, , drop = FALSE]
      }
    }
  }
  out[[length(out) + 1L]] <- current
  do.call(rbind, out)
}

thresholds <- c(0, 10000, 100000)
boundary_count <- function(run_list, shuffled = FALSE) {
  total <- numeric(length(thresholds))
  for (z in run_list) {
    category <- z$category
    if (shuffled) category <- sample(category)
    if (length(category) < 2L) next
    gaps <- z$start[-1] - z$end[-nrow(z)]
    cross <- (category[-length(category)] == "ACRO" & category[-1] == "HSAT1A") |
      (category[-length(category)] == "HSAT1A" & category[-1] == "ACRO")
    total <- total + vapply(thresholds, function(t) sum(cross & gaps <= t), numeric(1))
  }
  total
}

adjacency_rows <- list()
for (i in seq_len(nrow(metadata))) {
  m <- metadata[i, ]
  d <- bed_by_tip[[m$tip]]
  primary_contigs <- fai_by_tip[[m$tip]]$chrom[
    grepl(m$primary_pattern, fai_by_tip[[m$tip]]$chrom)
  ]
  d <- d[d$chrom %in% primary_contigs, , drop = FALSE]
  if (!any(d$category == "HSAT1A")) {
    for (threshold in thresholds) {
      adjacency_rows[[length(adjacency_rows) + 1L]] <- data.frame(
        tip = m$tip, species = m$species, max_gap_bp = threshold,
        observed = 0, null_mean = NA_real_, enrichment_fold = NA_real_,
        p_enrichment = NA_real_, p_depletion = NA_real_, estimable = FALSE
      )
    }
    next
  }
  runs <- lapply(split(d, d$chrom), collapse_runs)
  observed <- boundary_count(runs, shuffled = FALSE)
  null <- matrix(0, nrow = n_perm, ncol = length(thresholds))
  for (b in seq_len(n_perm)) null[b, ] <- boundary_count(runs, shuffled = TRUE)
  for (k in seq_along(thresholds)) {
    adjacency_rows[[length(adjacency_rows) + 1L]] <- data.frame(
      tip = m$tip, species = m$species, max_gap_bp = thresholds[[k]],
      observed = observed[[k]], null_mean = mean(null[, k]),
      enrichment_fold = observed[[k]] / mean(null[, k]),
      p_enrichment = (sum(null[, k] >= observed[[k]]) + 1) / (n_perm + 1),
      p_depletion = (sum(null[, k] <= observed[[k]]) + 1) / (n_perm + 1),
      estimable = TRUE
    )
  }
}
adjacency <- do.call(rbind, adjacency_rows)
write.table(
  adjacency,
  file.path(output_dir, "local_adjacency_tests.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Figures.
library(ggplot2)
tip_colors <- c(
  Human = "#1B9E77", Gorilla = "#6A3D9A", Chimpanzee = "#E31A1C",
  Bonobo = "#1F78B4", Sumatran_orangutan = "#FF7F00", Bornean_orangutan = "#B15928"
)
species$primary_acro_mb <- species$primary_acrocentric_acro_bp / 1e6
species$primary_hsat_mb <- species$primary_acrocentric_hsat1a_bp / 1e6
species$label_y <- species$primary_hsat_mb + 1.5
species$label_y[species$tip == "Sumatran_orangutan"] <- 2.4
species$label_y[species$tip == "Bornean_orangutan"] <- 0.8

p_primary <- ggplot(
  species,
  aes(x = primary_acro_mb, y = primary_hsat_mb, color = tip, label = species)
) +
  geom_point(size = 4) +
  geom_text(aes(y = label_y), size = 3.3, show.legend = FALSE) +
  scale_color_manual(values = tip_colors) +
  expand_limits(x = c(0, 15), y = c(0, 58)) +
  labs(
    x = "ACRO annotation on primary acrocentric haplotypes (Mb)",
    y = "HSAT1A annotation on primary acrocentric haplotypes (Mb)",
    color = "Species",
    title = "Quantitative ACRO-HSAT1A co-enrichment across six assemblies",
    subtitle = sprintf(
      "Exact tip permutation rho = %.3f, p = %.4g; phylogenetic PGLS p = %.3f",
      primary_rank[["rho"]], primary_rank[["p_one_sided"]],
      pgls_results$p_two_sided[pgls_results$test == "PGLS Brownian, Grafen branch lengths"]
    )
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "figures", "six_species_primary_acro_hsat1a.png"),
  p_primary, width = 8.8, height = 5.8, dpi = 300)
ggsave(file.path(output_dir, "figures", "six_species_primary_acro_hsat1a.pdf"),
  p_primary, width = 8.8, height = 5.8)

plot_adj <- adjacency[adjacency$max_gap_bp == 0 & adjacency$estimable, ]
p_adj <- ggplot(plot_adj, aes(x = reorder(species, enrichment_fold), y = enrichment_fold, fill = tip)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey45") +
  geom_col(width = 0.68) +
  scale_fill_manual(values = tip_colors) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Observed / permuted-mean direct ACRO-HSAT1A boundaries",
    title = "Local adjacency is distinct from abundance co-enrichment",
    subtitle = "Orangutans are not estimable because HSAT1A is absent"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "figures", "direct_adjacency_by_species.png"),
  p_adj, width = 8.2, height = 4.8, dpi = 300)
ggsave(file.path(output_dir, "figures", "direct_adjacency_by_species.pdf"),
  p_adj, width = 8.2, height = 4.8)

parameters <- data.frame(
  parameter = c(
    "input_directory", "ACRO_definition", "HSAT1A_definition",
    "species_level_tree", "PGLS_primary_branch_lengths",
    "permutations", "seed"
  ),
  value = c(
    normalizePath(input_dir), "annotation contains ACRO1", "annotation equals HSat1A",
    write.tree(topology), "Grafen", n_perm, seed
  ),
  stringsAsFactors = FALSE
)
write.table(parameters, file.path(output_dir, "analysis_parameters.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat("Six-species analysis complete\n")
cat("Primary-acrocentric exact rho:", primary_rank[["rho"]],
  "one-sided p:", primary_rank[["p_one_sided"]], "\n")
cat("Genome-wide exact rho:", genome_rank[["rho"]],
  "one-sided p:", genome_rank[["p_one_sided"]], "\n")
cat("Grafen PGLS two-sided p:",
  pgls_results$p_two_sided[pgls_results$test == "PGLS Brownian, Grafen branch lengths"], "\n")
cat("Species x homolog permutation p:", homolog_p_one, "\n")
cat("Matched haplotype exact p:", pair_p_one, "\n")
cat("Direct adjacency:\n")
print(adjacency[adjacency$max_gap_bp == 0, c(
  "species", "observed", "null_mean", "enrichment_fold", "p_enrichment", "p_depletion"
)], row.names = FALSE)
