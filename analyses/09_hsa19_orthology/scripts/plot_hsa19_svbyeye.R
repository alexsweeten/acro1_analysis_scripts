#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath("phylogeny_origin_acro1_full")
target <- if (length(args) >= 2) args[[2]] else "both"
if (!target %in% c("both", "all14", "representative9")) {
  stop("Target must be one of: both, all14, representative9")
}
local_lib <- file.path(root, "R_library")
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(SVbyEye)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
})

work <- file.path(root, "svbyeye")
figures <- file.path(root, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)

all_order <- c(
  "Human", "Chimpanzee", "Bonobo", "Gorilla", "Sumatran_orangutan",
  "Bornean_orangutan", "Siamang_hap1", "Siamang_hap2", "Hoolock_gibbon",
  "Nomascus_gibbon", "Rhesus_macaque", "Marmoset",
  "Golden_lion_tamarin", "Squirrel_monkey"
)

representative_order <- c(
  "Human", "Chimpanzee", "Gorilla", "Sumatran_orangutan", "Siamang_hap1",
  "Nomascus_gibbon", "Rhesus_macaque", "Marmoset", "Squirrel_monkey"
)

display_names <- c(
  Human = "Human",
  Chimpanzee = "Chimpanzee",
  Bonobo = "Bonobo",
  Gorilla = "Gorilla",
  Sumatran_orangutan = "Sumatran orangutan",
  Bornean_orangutan = "Bornean orangutan",
  Siamang_hap1 = "Siamang, haplotype 1",
  Siamang_hap2 = "Siamang, haplotype 2",
  Hoolock_gibbon = "Hoolock gibbon",
  Nomascus_gibbon = "Nomascus gibbon",
  Rhesus_macaque = "Rhesus macaque",
  Marmoset = "Common marmoset",
  Golden_lion_tamarin = "Golden lion tamarin",
  Squirrel_monkey = "Squirrel monkey"
)

feature_colors <- c(
  "Inferred orthologous interval" = "#F5F5F5",
  "ACRO1 profile-HMM tract" = "#7349A4",
  "Non-ACRO ACRE component" = "#E69F00",
  "Associated mobile-element component" = "#7C8798"
)

annotations <- read.delim(
  file.path(work, "HSA19_svbyeye_annotations.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

make_plot <- function(
  paf_name,
  order,
  basename,
  title_suffix,
  height,
  show_title = TRUE,
  legend_multiplier = 1
) {
  paf_file <- file.path(work, paf_name)
  paf <- SVbyEye::readPaf(
    paf.file = paf_file,
    include.paf.tags = TRUE,
    restrict.paf.tags = "cg"
  )
  paf$identity <- paf$n.match / paf$aln.len
  paf <- paf[paf$aln.len >= 3000 & paf$identity >= 0.65, ]
  if (nrow(paf) == 0) stop("No PAF records remained after filtering: ", paf_file)

  plt <- SVbyEye::plotAVA(
    paf.table = paf,
    seqnames.order = order,
    binsize = 5000,
    color.by = "identity",
    perc.identity.breaks = c(70, 80, 85, 90, 95, 98, 99),
    outline.alignments = FALSE
  )

  annot <- annotations[annotations$ID %in% order, ]
  annot$feature <- factor(annot$feature, levels = names(feature_colors))
  annot_gr <- GenomicRanges::GRanges(
    seqnames = annot$seqnames,
    ranges = IRanges::IRanges(start = annot$start, end = annot$end),
    feature = annot$feature,
    ID = annot$ID
  )
  plt <- SVbyEye::addAnnotation(
    ggplot.obj = plt,
    annot.gr = annot_gr,
    shape = "rectangle",
    fill.by = "feature",
    color.palette = feature_colors,
    coordinate.space = "self",
    y.label.id = "ID",
    annotation.level = 0
  )

  n <- length(order)
  axis_breaks <- seq(100000, 1300000, 200000)
  axis_labels <- sprintf("%+d", (axis_breaks - 700000) / 1000)
  plt <- plt +
    ggplot2::geom_vline(xintercept = 700000, linewidth = 0.35, linetype = 3, color = "#404040") +
    ggplot2::scale_x_continuous(
      breaks = axis_breaks,
      labels = axis_labels,
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = c(1, 1400000), clip = "off") +
    ggplot2::scale_y_continuous(
      breaks = n:1,
      labels = unname(display_names[order]),
      expand = ggplot2::expansion(mult = c(0.035, 0.06))
    ) +
    ggplot2::labs(
      title = if (show_title) paste0("HSA19-orthologous ACRO1 locus - SVbyEye ", title_suffix) else NULL,
      subtitle = if (show_title) {
        "Progressive minimap2 alignments of 1.4-Mb orthologous windows; tracks oriented to human"
      } else {
        NULL
      },
      x = "Distance from inferred locus center (kb)",
      y = NULL
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title = "feature",
        nrow = if (legend_multiplier > 1) 2 else 1,
        byrow = TRUE,
        title.position = "left"
      ),
      fill_ggnewscale_1 = ggplot2::guide_legend(
        title = "Identity",
        nrow = 2,
        byrow = TRUE,
        title.position = "left"
      )
    ) +
    ggplot2::theme_minimal(base_size = 30) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35),
      axis.text.y = ggplot2::element_text(size = 27, color = "#24262B"),
      axis.text.x = ggplot2::element_text(size = 24),
      axis.title.x = ggplot2::element_text(size = 30, margin = ggplot2::margin(t = 24)),
      plot.title = if (show_title) {
        ggplot2::element_text(size = 45, face = "bold", color = "#16181D")
      } else {
        ggplot2::element_blank()
      },
      plot.subtitle = if (show_title) {
        ggplot2::element_text(size = 28.5, color = "#4B5563", margin = ggplot2::margin(b = 36))
      } else {
        ggplot2::element_blank()
      },
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 25.5 * legend_multiplier, face = "bold"),
      legend.text = ggplot2::element_text(size = 24 * legend_multiplier),
      plot.margin = ggplot2::margin(if (show_title) 54 else 24, 66, 36, 54)
    )

  for (extension in c("png", "pdf")) {
    ggplot2::ggsave(
      filename = file.path(figures, paste0(basename, ".", extension)),
      plot = plt,
      width = 33,
      height = height,
      units = "in",
      dpi = if (extension == "png") 160 else 300,
      bg = "white",
      limitsize = FALSE
    )
  }
  pdftocairo <- Sys.which("pdftocairo")
  if (nzchar(pdftocairo)) {
    font_cache <- file.path(work, ".fontconfig-cache")
    dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
    old_cache <- Sys.getenv("XDG_CACHE_HOME", unset = NA_character_)
    Sys.setenv(XDG_CACHE_HOME = font_cache)
    status <- system2(
      pdftocairo,
      args = c(
        "-svg",
        file.path(figures, paste0(basename, ".pdf")),
        file.path(figures, paste0(basename, ".svg"))
      )
    )
    if (is.na(old_cache)) Sys.unsetenv("XDG_CACHE_HOME") else Sys.setenv(XDG_CACHE_HOME = old_cache)
    if (status != 0) warning("pdftocairo could not create SVG for ", basename)
  } else {
    warning("pdftocairo not found; PNG and PDF were written, but SVG was skipped")
  }
  sips <- Sys.which("sips")
  if (nzchar(sips)) {
    status <- system2(
      sips,
      args = c(
        "-s", "format", "com.adobe.photoshop-image",
        file.path(figures, paste0(basename, ".png")),
        "--out", file.path(figures, paste0(basename, ".psd"))
      ),
      stdout = FALSE,
      stderr = FALSE
    )
    if (status != 0) warning("sips could not create PSD for ", basename)
  } else {
    warning("sips not found; PSD was skipped")
  }

  summary <- data.frame(
    plot = basename,
    SVbyEye_version = as.character(utils::packageVersion("SVbyEye")),
    PAF_records_read = nrow(SVbyEye::readPaf(paf.file = paf_file, include.paf.tags = FALSE)),
    PAF_records_plotted = nrow(paf),
    minimum_alignment_length = 3000,
    minimum_identity = 0.65,
    bin_size = 5000,
    sequence_count = length(order)
  )
  return(summary)
}

new_summaries <- list()
if (target %in% c("both", "all14")) {
  new_summaries[[length(new_summaries) + 1]] <- make_plot(
    "HSA19_all14_progressive.paf",
    all_order,
    "HSA19_ACRO1_synteny_SVbyEye_all14",
    "(all 14 assemblies)",
    24.6,
    show_title = FALSE,
    legend_multiplier = 2
  )
}

if (target %in% c("both", "representative9")) {
  new_summaries[[length(new_summaries) + 1]] <- make_plot(
    "HSA19_representative9_progressive.paf",
    representative_order,
    "HSA19_ACRO1_synteny_SVbyEye_representative9",
    "(representative panel)",
    17.6
  )
}

new_summary <- do.call(rbind, new_summaries)
summary_path <- file.path(work, "SVbyEye_plot_summary.tsv")
if (file.exists(summary_path)) {
  old_summary <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
  old_summary <- old_summary[!old_summary$plot %in% new_summary$plot, , drop = FALSE]
  new_summary <- rbind(old_summary, new_summary)
}
write.table(
  new_summary,
  summary_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Wrote SVbyEye figures to: ", figures)
