#!/usr/bin/env Rscript
# Cluster-aware group contrast: effect ~ group + (1 | individual).
# Inference: (a) lmerTest Satterthwaite, (b) OLS with individual-clustered sandwich SEs,
#            (c) individual-cluster bootstrap of the lmer fixed effect (whole individuals resampled).
# Args: input.tsv output.tsv seed n_boot
suppressPackageStartupMessages({library(lme4); library(lmerTest); library(sandwich)})
args <- commandArgs(trailingOnly = TRUE)
d <- read.delim(args[1], stringsAsFactors = FALSE)
seed <- as.integer(args[3]); B <- as.integer(args[4])
ctrl <- lmerControl(check.conv.singular = "ignore", check.conv.grad = "ignore", check.conv.hess = "ignore")
fit_beta <- function(dat) {
  m <- tryCatch(suppressMessages(suppressWarnings(lmer(effect ~ grp + (1 | ind), data = dat, REML = TRUE, control = ctrl))),
                error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  fixef(m)[["grpA"]]
}
out <- list()
for (an in unique(d$analysis)) {
  s <- d[d$analysis == an, ]
  a_name <- setdiff(unique(s$group), unique(s$reference))[1]; b_name <- unique(s$reference)[1]
  dat <- data.frame(effect = s$effect, grp = factor(ifelse(s$group == a_name, "A", "B"), levels = c("B", "A")),
                    ind = factor(s$individual_id))
  both <- sum(rowSums(table(dat$ind, dat$grp) > 0) == 2)
  m <- suppressMessages(suppressWarnings(lmer(effect ~ grp + (1 | ind), data = dat, REML = TRUE, control = ctrl)))
  co <- summary(m)$coefficients["grpA", ]
  vc <- as.data.frame(VarCorr(m))
  icc <- vc$vcov[vc$grp == "ind"] / sum(vc$vcov)
  base <- list(analysis = an, group = paste(a_name, "-", b_name), n_arrays = nrow(dat), n_clusters = nlevels(dat$ind),
               n_clusters_with_both_groups = both, cluster_unit = "individual")
  out[[length(out) + 1]] <- c(base, list(test = "group_contrast_lmer_random_intercept_Satterthwaite", estimate = unname(co["Estimate"]),
       se = unname(co["Std. Error"]), p_two_sided = unname(co["Pr(>|t|)"]),
       detail = sprintf("REML; individual ICC=%.3f; singular=%s", icc, isSingular(m))))
  ols <- lm(effect ~ grp, data = dat)
  vcl <- vcovCL(ols, cluster = dat$ind, type = "HC1")
  se <- sqrt(vcl["grpA", "grpA"]); est <- coef(ols)[["grpA"]]; G <- nlevels(dat$ind)
  out[[length(out) + 1]] <- c(base, list(test = "group_contrast_OLS_cluster_robust_by_individual", estimate = est, se = se,
       p_two_sided = 2 * pt(-abs(est / se), df = G - 1), detail = sprintf("CR1 sandwich, t with G-1=%d df", G - 1)))
  set.seed(seed + sum(utf8ToInt(an)))
  ids <- levels(dat$ind); by_ind <- split(seq_len(nrow(dat)), dat$ind)
  bs <- numeric(B)
  for (b in seq_len(B)) {
    pick <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(by_ind[pick], use.names = FALSE)
    bd <- dat[idx, ]
    bd$ind <- factor(rep(seq_along(pick), lengths(by_ind[pick])))
    bs[b] <- fit_beta(bd)
  }
  bs <- bs[is.finite(bs)]
  pb <- (1 + min(sum(bs <= 0), sum(bs >= 0))) * 2 / (length(bs) + 1)
  out[[length(out) + 1]] <- c(base, list(test = "group_contrast_lmer_individual_cluster_bootstrap", estimate = unname(co["Estimate"]),
       boot_ci95_lo = unname(quantile(bs, 0.025)), boot_ci95_hi = unname(quantile(bs, 0.975)), p_two_sided = min(1, pb),
       signflip_draws = length(bs), seed_base = seed,
       detail = "percentile CI; two-sided p = 2*min(tail)+plus-one on bootstrap distribution of the fixed effect"))
}
cols <- unique(unlist(lapply(out, names)))
res <- do.call(rbind, lapply(out, function(o) { o <- o[cols]; names(o) <- cols; as.data.frame(lapply(o, function(v) if (is.null(v)) NA else v), stringsAsFactors = FALSE) }))
write.table(res, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
