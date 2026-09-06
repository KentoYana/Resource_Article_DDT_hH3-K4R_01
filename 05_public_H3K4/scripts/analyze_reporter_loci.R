# analyze_reporter_loci.R
# Kento Yanagisawa
# This script analyzes public H3K4 ChIP-seq signal at candidate reporter loci.

# Make these packages and their associated functions
# available to use in this script
library("tikzDevice")
library("tidyverse")
library("patchwork")
library("here")
library("IRanges")
library("rtracklayer")
library("digest")
library("jsonlite")

# Clear R's brain
rm(list = ls())

# Set project root
here::i_am("05_public_H3K4/scripts/analyze_reporter_loci.R")

# =========================
# Functions
# =========================

# Parse optional command-line arguments used by the shell workflow
parse_arguments <- function(arguments, default_work_root, default_output_dir) {
  result <- list(
    work_root = default_work_root,
    output_dir = default_output_dir
  )

  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]

    if (argument %in% c("-h", "--help")) {
      cat(
        "Usage: analyze_reporter_loci.R [--work-root PATH] [--output-dir PATH]\n\n",
        "Defaults:\n",
        "  --work-root  ", default_work_root, "\n",
        "  --output-dir ", default_output_dir, "\n",
        sep = ""
      )
      quit(status = 0L)
    }

    if (!argument %in% c("--work-root", "--output-dir") ||
        index == length(arguments)) {
      stop("Unknown or incomplete argument: ", argument, call. = FALSE)
    }

    key <- argument %>%
      stringr::str_remove("^--") %>%
      stringr::str_replace_all("-", "_")
    result[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }

  result
}

# Extract one key from the semicolon-delimited GFF attribute column
extract_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  stringr::str_match(attributes, pattern)[, 2]
}

# Read protein-coding genes on the seven NC12 nuclear chromosomes
load_genes <- function(gff_path) {
  column_names <- c(
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "phase", "attributes"
  )

  gff <- readr::read_tsv(
    gff_path,
    comment = "#",
    col_names = column_names,
    col_types = readr::cols(
      chrom = readr::col_character(),
      source = readr::col_character(),
      feature = readr::col_character(),
      start = readr::col_integer(),
      end = readr::col_integer(),
      score = readr::col_character(),
      strand = readr::col_character(),
      phase = readr::col_character(),
      attributes = readr::col_character()
    ),
    progress = FALSE,
    show_col_types = FALSE
  )

  genes <- gff %>%
    dplyr::mutate(
      locus_tag = extract_attribute(attributes, "locus_tag"),
      gene_name = extract_attribute(attributes, "gene"),
      feature_name = extract_attribute(attributes, "Name"),
      gene_biotype = extract_attribute(attributes, "gene_biotype")
    ) %>%
    dplyr::filter(
      feature == "gene",
      gene_biotype == "protein_coding",
      stringr::str_detect(chrom, "^CM002(23[6-9]|24[0-2])\\.1$"),
      !is.na(locus_tag)
    ) %>%
    dplyr::transmute(
      locus_tag,
      symbol = dplyr::coalesce(gene_name, feature_name, locus_tag),
      chrom,
      start,
      end,
      strand,
      length_bp = end - start + 1L
    )

  if (anyDuplicated(genes$locus_tag)) {
    stop("Duplicate protein-coding gene features in GFF", call. = FALSE)
  }

  genes
}

# Return one gene row for a locus tag
gene_lookup <- function(genes, target_locus_tag) {
  gene <- genes %>%
    dplyr::filter(locus_tag == target_locus_tag)

  if (nrow(gene) != 1L) {
    stop(
      "Gene lookup did not return exactly one row: ",
      target_locus_tag,
      call. = FALSE
    )
  }

  gene
}

# Define promoter, gene-body, or flanking regions
metric_regions <- function(
  genes,
  metric,
  chromosome_lengths,
  window_flank_bp,
  promoter_upstream_bp,
  promoter_downstream_bp
) {
  if (metric == "gene_body") {
    regions <- genes %>%
      dplyr::transmute(locus_tag, chrom, start, end)
  } else if (metric == "gene_body_plus_minus_2kb") {
    regions <- genes %>%
      dplyr::transmute(
        locus_tag,
        chrom,
        start = start - window_flank_bp,
        end = end + window_flank_bp
      )
  } else if (metric == "promoter") {
    regions <- genes %>%
      dplyr::mutate(
        region_start = dplyr::if_else(
          strand == "+",
          start - promoter_upstream_bp,
          end - promoter_downstream_bp + 1L
        ),
        region_end = dplyr::if_else(
          strand == "+",
          start + promoter_downstream_bp - 1L,
          end + promoter_upstream_bp
        )
      ) %>%
      dplyr::transmute(
        locus_tag,
        chrom,
        start = region_start,
        end = region_end
      )
  } else {
    stop("Unknown metric: ", metric, call. = FALSE)
  }

  regions %>%
    dplyr::mutate(
      start = pmax(1L, as.integer(start)),
      end = pmin(
        as.integer(unname(chromosome_lengths[chrom])),
        as.integer(end)
      )
    )
}

# Apply one summary statistic to genomic intervals chromosome by chromosome
view_statistic <- function(signal, regions, statistic) {
  result <- rep(NA_real_, nrow(regions))

  for (chromosome in unique(regions$chrom)) {
    indices <- which(regions$chrom == chromosome)
    views <- IRanges::Views(
      signal[[chromosome]],
      start = regions$start[indices],
      end = regions$end[indices]
    )
    result[indices] <- statistic(views)
  }

  result
}

# Calculate mean CPM for each interval
region_means <- function(signal, regions) {
  view_statistic(
    signal,
    regions,
    function(views) IRanges::viewMeans(views, na.rm = TRUE)
  )
}

# Calculate an empirical midrank percentile and descending rank
empirical_midrank <- function(value, population) {
  finite <- population[is.finite(population)]

  if (!is.finite(value) || length(finite) == 0L) {
    return(c(
      percentile = NA_real_,
      rank_desc = NA_real_,
      population_n = length(finite)
    ))
  }

  less <- sum(finite < value)
  greater <- sum(finite > value)
  ties <- sum(finite == value)

  c(
    percentile = 100 * (less + 0.5 * ties) / length(finite),
    rank_desc = greater + (ties + 1) / 2,
    population_n = length(finite)
  )
}

# Summarize coverage, exact zeros, and signal in one reporter region
region_details <- function(signal, region) {
  values <- as.numeric(
    signal[[region$chrom]][region$start:region$end]
  )
  valid_values <- values[is.finite(values)]
  region_length <- length(values)
  covered <- length(valid_values)

  if (covered == 0L) {
    return(c(
      region_start_1based = region$start,
      region_end_1based = region$end,
      region_length_bp = region_length,
      covered_bases = 0,
      missing_bases = region_length,
      exact_zero_bases = 0,
      mean_cpm = NA_real_,
      max_cpm = NA_real_,
      nonzero_fraction = NA_real_
    ))
  }

  c(
    region_start_1based = region$start,
    region_end_1based = region$end,
    region_length_bp = region_length,
    covered_bases = covered,
    missing_bases = region_length - covered,
    exact_zero_bases = sum(valid_values == 0),
    mean_cpm = mean(valid_values),
    max_cpm = max(valid_values),
    nonzero_fraction = sum(valid_values > 0) / covered
  )
}

# Bin a strand-oriented profile around one gene
profile_for_gene <- function(signal, gene, window_flank_bp, profile_bin_bp) {
  region_start <- max(1L, gene$start - window_flank_bp)
  region_end <- min(
    length(signal[[gene$chrom]]),
    gene$end + window_flank_bp
  )
  positions <- seq.int(region_start, region_end)
  values <- as.numeric(signal[[gene$chrom]][region_start:region_end])

  if (gene$strand == "+") {
    relative_positions <- positions - gene$start
  } else {
    values <- rev(values)
    relative_positions <- rev(gene$end - positions)
  }

  tibble::tibble(
    relative_position = relative_positions,
    cpm = values,
    profile_bin = ceiling(seq_along(values) / profile_bin_bp)
  ) %>%
    dplyr::group_by(profile_bin) %>%
    dplyr::summarise(
      x = mean(relative_position),
      y = if (all(is.na(cpm))) 0 else mean(cpm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(x, y)
}

# Return one summary row
summary_lookup <- function(
  summary,
  sample,
  variant,
  target_locus_tag,
  region_metric
) {
  selected <- summary %>%
    dplyr::filter(
      sample_id == sample,
      track_variant == variant,
      locus_tag == target_locus_tag,
      metric == region_metric
    )

  if (nrow(selected) != 1L) {
    stop("Summary lookup did not return exactly one row", call. = FALSE)
  }

  selected
}

# Quantify reporter regions and collect profile data from every bigWig
collect_analysis_data <- function(
  genes,
  bigwig_root,
  targets,
  runs,
  variants,
  metrics,
  window_flank_bp,
  promoter_upstream_bp,
  promoter_downstream_bp,
  profile_bin_bp
) {
  summary_rows <- list()
  profile_rows <- list()
  summary_index <- 1L
  profile_index <- 1L

  target_genes <- targets %>%
    dplyr::select(locus_tag) %>%
    dplyr::left_join(genes, by = "locus_tag")

  for (variant in names(variants)) {
    for (run_index in seq_len(nrow(runs))) {
      run <- runs %>% dplyr::slice(run_index)
      bigwig_name <- paste0(run$sample_id, ".", variants[[variant]])
      bigwig_path <- file.path(bigwig_root, bigwig_name)

      if (!file.exists(bigwig_path)) {
        stop("Missing bigWig: ", bigwig_path, call. = FALSE)
      }

      signal <- rtracklayer::import(bigwig_path, as = "RleList")
      chromosome_lengths <- lengths(signal)

      if (!all(unique(genes$chrom) %in% names(signal))) {
        stop(
          "A nuclear chromosome is absent from ",
          bigwig_name,
          call. = FALSE
        )
      }

      genome_values <- list()
      target_values <- list()
      target_regions <- list()

      for (metric in metrics) {
        regions <- metric_regions(
          genes,
          metric,
          chromosome_lengths,
          window_flank_bp,
          promoter_upstream_bp,
          promoter_downstream_bp
        )
        values <- region_means(signal, regions)
        names(values) <- regions$locus_tag

        genome_values[[metric]] <- values
        target_values[[metric]] <- values[targets$locus_tag]
        target_regions[[metric]] <- metric_regions(
          target_genes,
          metric,
          chromosome_lengths,
          window_flank_bp,
          promoter_upstream_bp,
          promoter_downstream_bp
        )
      }

      for (target_index in seq_len(nrow(targets))) {
        target <- targets %>% dplyr::slice(target_index)
        gene <- gene_lookup(target_genes, target$locus_tag)

        profile_rows[[profile_index]] <- profile_for_gene(
          signal,
          gene,
          window_flank_bp,
          profile_bin_bp
        ) %>%
          dplyr::mutate(
            sample_id = run$sample_id,
            study = run$study,
            mark = run$mark,
            track_variant = variant,
            display_name = target$display_name,
            locus_tag = target$locus_tag,
            strand = gene$strand,
            gene_length_bp = gene$length_bp,
            .before = 1
          )
        profile_index <- profile_index + 1L

        for (metric in metrics) {
          value <- target_values[[metric]][[target$locus_tag]]
          genome_rank <- empirical_midrank(value, genome_values[[metric]])
          reporter_rank <- empirical_midrank(value, target_values[[metric]])
          target_region <- target_regions[[metric]] %>%
            dplyr::filter(locus_tag == target$locus_tag)
          details <- region_details(signal, target_region)

          summary_rows[[summary_index]] <- tibble::tibble(
            sample_id = run$sample_id,
            study = run$study,
            mark = run$mark,
            track_variant = variant,
            bigwig_file = bigwig_name,
            display_name = target$display_name,
            locus_tag = target$locus_tag,
            role = target$role,
            chrom = gene$chrom,
            gene_start_1based = gene$start,
            gene_end_1based = gene$end,
            strand = gene$strand,
            metric = metric,
            region_start_1based = details[["region_start_1based"]],
            region_end_1based = details[["region_end_1based"]],
            region_length_bp = details[["region_length_bp"]],
            covered_bases = details[["covered_bases"]],
            missing_bases = details[["missing_bases"]],
            exact_zero_bases = details[["exact_zero_bases"]],
            mean_cpm = details[["mean_cpm"]],
            max_cpm = details[["max_cpm"]],
            nonzero_fraction = details[["nonzero_fraction"]],
            genome_percentile_midrank = genome_rank[["percentile"]],
            genome_rank_desc = genome_rank[["rank_desc"]],
            genome_population_n = genome_rank[["population_n"]],
            reporter_set_rank_desc = reporter_rank[["rank_desc"]],
            reporter_set_n = reporter_rank[["population_n"]]
          )
          summary_index <- summary_index + 1L
        }
      }

      rm(signal)
      invisible(gc())
    }
  }

  list(
    summary = dplyr::bind_rows(summary_rows),
    profiles = dplyr::bind_rows(profile_rows)
  )
}

# Make one profile panel for a run and reporter locus
make_profile_panel <- function(
  profile_data,
  summary,
  run,
  target,
  gene,
  variant,
  y_max,
  show_title,
  show_y_axis,
  show_left_flank_label,
  show_right_flank_label,
  window_flank_bp,
  promoter_upstream_bp,
  promoter_downstream_bp,
  profile_color,
  promoter_color
) {
  body <- summary_lookup(
    summary,
    run$sample_id,
    variant,
    target$locus_tag,
    "gene_body"
  )
  promoter <- summary_lookup(
    summary,
    run$sample_id,
    variant,
    target$locus_tag,
    "promoter"
  )

  title_text <- if (show_title) {
    paste0(
      target$display_name,
      "\n",
      target$locus_tag,
      " (",
      gene$strand,
      ")"
    )
  } else {
    NULL
  }
  y_axis_text <- if (show_y_axis) {
    paste(run$study, run$mark, "CPM", sep = "\n")
  } else {
    NULL
  }
  annotation_text <- paste0(
    "body pct ", round(body$genome_percentile_midrank),
    "; promoter pct ", round(promoter$genome_percentile_midrank)
  )

  panel <- ggplot2::ggplot(profile_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::annotate(
      "rect",
      xmin = -promoter_upstream_bp,
      xmax = promoter_downstream_bp,
      ymin = 0,
      ymax = Inf,
      fill = promoter_color,
      alpha = 0.14
    ) +
    ggplot2::annotate(
      "rect",
      xmin = 0,
      xmax = gene$length_bp,
      ymin = 0,
      ymax = Inf,
      fill = "grey70",
      alpha = 0.25
    ) +
    ggplot2::geom_area(fill = profile_color, alpha = 0.45) +
    ggplot2::geom_line(color = profile_color, linewidth = 0.5) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey35",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.05,
      vjust = 1.25,
      size = 1.7,
      color = "grey25"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-window_flank_bp, gene$length_bp + window_flank_bp),
      breaks = c(
        -window_flank_bp,
        0,
        gene$length_bp,
        gene$length_bp + window_flank_bp
      ),
      labels = c(
        if (show_left_flank_label) "-2 kb" else "",
        "TSS",
        "TES",
        if (show_right_flank_label) "+2 kb" else ""
      ),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max),
      breaks = c(0, y_max),
      labels = function(values) format(round(values, 1), trim = TRUE),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(title = title_text, x = NULL, y = y_axis_text) +
    ggplot2::theme_bw(base_size = 6) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      axis.ticks = ggplot2::element_line(linewidth = 0.5),
      axis.text = ggplot2::element_text(size = 5, colour = "black"),
      axis.title.y = ggplot2::element_text(
        size = 6,
        margin = ggplot2::margin(r = 2)
      ),
      plot.title = ggplot2::element_text(
        size = 6,
        hjust = 0.5,
        lineheight = 0.9,
        colour = "black"
      ),
      plot.margin = ggplot2::margin(1, 1, 1, 1)
    )

  if (!show_y_axis) {
    panel <- panel +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  panel
}

# Assemble one mark-specific reporter-locus profile figure
make_profile_plot <- function(
  profiles,
  summary,
  genes,
  targets,
  runs,
  variant,
  selected_mark,
  window_flank_bp,
  promoter_upstream_bp,
  promoter_downstream_bp,
  plot_colors
) {
  profile_color <- plot_colors[["blue"]]
  mark_runs <- runs %>%
    dplyr::filter(mark == selected_mark)
  if (nrow(mark_runs) == 0L) {
    stop("No runs found for histone mark: ", selected_mark, call. = FALSE)
  }
  row_plots <- vector("list", nrow(mark_runs))

  for (run_index in seq_len(nrow(mark_runs))) {
    run <- mark_runs %>% dplyr::slice(run_index)
    run_profiles <- profiles %>%
      dplyr::filter(
        sample_id == run$sample_id,
        track_variant == variant
      )
    y_max <- max(run_profiles$y, na.rm = TRUE) * 1.05
    if (!is.finite(y_max) || y_max <= 0) {
      y_max <- 1
    }

    panels <- vector("list", nrow(targets))
    for (target_index in seq_len(nrow(targets))) {
      target <- targets %>% dplyr::slice(target_index)
      gene <- gene_lookup(genes, target$locus_tag)
      profile_data <- run_profiles %>%
        dplyr::filter(locus_tag == target$locus_tag)

      panels[[target_index]] <- make_profile_panel(
        profile_data,
        summary,
        run,
        target,
        gene,
        variant,
        y_max,
        show_title = run_index == 1L,
        show_y_axis = target_index == 1L,
        show_left_flank_label = target_index == 1L,
        show_right_flank_label = target_index == nrow(targets),
        window_flank_bp,
        promoter_upstream_bp,
        promoter_downstream_bp,
        profile_color,
        plot_colors[["orange"]]
      )
    }

    row_plots[[run_index]] <- patchwork::wrap_plots(panels, nrow = 1)
  }

  variant_label <- stringr::str_replace_all(variant, "_", " ")
  patchwork::wrap_plots(row_plots, ncol = 1) +
    patchwork::plot_annotation(
      title = paste0(
        "Candidate reporter loci: gene-oriented ",
        selected_mark,
        " profiles (",
        variant_label,
        ")"
      ),
      caption = paste0(
        "Each row has one shared y-axis across all seven loci; y-axes are not ",
        "shared between runs.\n",
        "Orange: promoter (-1 kb to +200 bp); grey: annotated gene body.\n",
        "Annotations are within-run genome-wide midrank percentiles."
      )
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 6, hjust = 0.5)
    )
}

# Make the percentile heatmap
make_percentile_plot <- function(
  summary,
  targets,
  runs,
  metrics,
  variant,
  plot_colors
) {
  metric_labels <- c(
    promoter = "Promoter (-1 kb/+200 bp)",
    gene_body = "Gene body",
    gene_body_plus_minus_2kb = "Gene body $\\pm$ 2 kb"
  )
  run_labels <- paste(runs$study, runs$mark, sep = "\n")

  plot_data <- summary %>%
    dplyr::filter(track_variant == variant) %>%
    dplyr::mutate(
      display_name = factor(display_name, levels = targets$display_name),
      run_label = factor(
        paste(study, mark, sep = "\n"),
        levels = rev(run_labels)
      ),
      metric_label = factor(
        dplyr::recode(metric, !!!metric_labels),
        levels = unname(metric_labels)
      ),
      percentile_label = round(genome_percentile_midrank),
      label_color = dplyr::if_else(
        genome_percentile_midrank < 18,
        "white",
        "black"
      )
    )

  variant_label <- stringr::str_replace_all(variant, "_", " ")
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = display_name,
      y = run_label,
      fill = genome_percentile_midrank
    )
  ) +
    ggplot2::geom_tile(color = "black", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = percentile_label, color = label_color),
      size = 2.4
    ) +
    ggplot2::facet_wrap(~metric_label, nrow = 1) +
    ggplot2::scale_fill_gradientn(
      colors = grDevices::colorRampPalette(
        unname(plot_colors),
        space = "Lab"
      )(100),
      limits = c(0, 100),
      name = "Within-run genome-wide\npercentile (midrank)"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        display = "rectangles",
        nbin = 20,
        direction = "horizontal"
      )
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = paste0(
        "Candidate reporter H3K4me signal percentiles (",
        variant_label,
        ")"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 6,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(size = 6, colour = "black"),
      axis.ticks = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(size = 8, face = "bold"),
      plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 6, colour = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key = ggplot2::element_blank(),
      legend.key.width = grid::unit(0.45, "cm"),
      legend.box.spacing = grid::unit(0, "cm")
    )
}

# Plot all figures first using the standard R graphics device
plot_standard_figures <- function(
  profile_plots,
  percentile_plots,
  figure_width_in,
  profile_heights_in
) {
  preview_path <- NULL

  if (!interactive()) {
    preview_path <- tempfile(
      pattern = "reporter_loci_preview_",
      fileext = ".pdf"
    )
    grDevices::pdf(
      preview_path,
      width = figure_width_in,
      height = max(profile_heights_in)
    )
    on.exit(
      {
        grDevices::dev.off()
        file.remove(preview_path)
      },
      add = TRUE
    )
  }

  for (variant in names(profile_plots)) {
    for (mark in names(profile_plots[[variant]])) {
      plot(profile_plots[[variant]][[mark]])
    }
  }
  plot(percentile_plots$nonduplicate)
  plot(percentile_plots$all_mapped)

  invisible(NULL)
}

# Switch to tikzDevice and draw the same plot as a tikzpicture fragment
write_tikz_plot <- function(plot_object, path, width, height) {
  figure_stem <- tools::file_path_sans_ext(basename(path))
  stale_sidecars <- list.files(dirname(path), full.names = TRUE) %>%
    purrr::keep(
      ~ stringr::str_starts(basename(.x), paste0(figure_stem, "_ras")) &&
        stringr::str_ends(basename(.x), ".png")
    )
  if (length(stale_sidecars) > 0L) {
    file.remove(stale_sidecars)
  }

  tikzDevice::tikz(
    path,
    width = width,
    height = height,
    standAlone = FALSE,
    timestamp = FALSE,
    sanitize = FALSE,
    lwdUnit = 72.27 / 96
  )
  plot(plot_object)
  grDevices::dev.off()

  generated_sidecars <- list.files(dirname(path), full.names = TRUE) %>%
    purrr::keep(
      ~ stringr::str_starts(basename(.x), paste0(figure_stem, "_ras")) &&
        stringr::str_ends(basename(.x), ".png")
    )
  if (length(generated_sidecars) > 0L) {
    stop("TikZ figure unexpectedly requires raster sidecars: ", basename(path),
         call. = FALSE)
  }
}

# Write reporter-locus coordinates
write_locus_table <- function(path, genes, targets) {
  locus_table <- targets %>%
    dplyr::left_join(genes, by = "locus_tag") %>%
    dplyr::transmute(
      display_name,
      locus_tag,
      role,
      chrom,
      start_1based = start,
      end_1based = end,
      strand,
      length_bp,
      gff_gene_name = symbol
    )

  readr::write_tsv(locus_table, path, na = "NA")
}

# Write the concise pan-2 interpretation table
write_analysis_summary <- function(path, summary, runs, metrics) {
  connection <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)

  writeLines(
    c(
      "# Reporter-locus H3K4me descriptive summary",
      "",
      "This table reports the nonduplicate tracks for `pan-2`. Percentiles are calculated separately for every run and region definition across all annotated protein-coding genes on the seven nuclear chromosomes.",
      "",
      "| Study | Mark | Region | Mean CPM | Nonzero fraction | Genome percentile | Reporter-set rank |",
      "|---|---|---|---:|---:|---:|---:|"
    ),
    connection
  )

  for (run_index in seq_len(nrow(runs))) {
    run <- runs %>% dplyr::slice(run_index)
    for (metric in metrics) {
      row <- summary_lookup(
        summary,
        run$sample_id,
        "nonduplicate",
        "NCU10048",
        metric
      )
      writeLines(
        sprintf(
          "| %s | %s | %s | %.4f | %.3f | %.1f | %.1f/7 |",
          run$study,
          run$mark,
          metric,
          row$mean_cpm,
          row$nonzero_fraction,
          row$genome_percentile_midrank,
          row$reporter_set_rank_desc
        ),
        connection
      )
    }
  }

  writeLines(
    c(
      "",
      "## Interpretation boundaries",
      "",
      "- These are descriptive CPM summaries, not statistical tests; no selected mark has biological replication within every study.",
      "- Nonzero coverage is not equivalent to enrichment because no matched input is used.",
      "- Studies and marks remain separate; absolute CPM values are not pooled or treated as directly exchangeable.",
      "- Exact zero bases are counted separately from missing bases in the long-form TSV.",
      "- `csr-1` is exploratory and is not treated as an established forward-mutation reporter."
    ),
    connection
  )
}

# Write the long-form summary with stable six-decimal numeric formatting
write_signal_summary <- function(path, summary) {
  six_decimal_columns <- c(
    "mean_cpm",
    "max_cpm",
    "nonzero_fraction",
    "genome_percentile_midrank",
    "genome_rank_desc",
    "reporter_set_rank_desc"
  )

  formatted <- summary %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(six_decimal_columns),
        ~ dplyr::if_else(
          is.na(.x),
          "NA",
          formatC(.x, format = "f", digits = 6L)
        )
      )
    )

  readr::write_tsv(formatted, path, na = "NA")
}

# Record software, inputs, and analysis choices
write_parameters <- function(
  path,
  script_path,
  gff_path,
  runs,
  targets,
  variants,
  profile_bin_bp,
  window_flank_bp,
  plot_colors,
  figure_width_in,
  profile_heights_in,
  percentile_height_in
) {
  input_bigwigs <- names(variants) %>%
    purrr::map(~ paste0(runs$sample_id, ".", variants[[.x]])) %>%
    unlist(use.names = FALSE)

  payload <- list(
    analysis_script = basename(script_path),
    R = as.character(getRversion()),
    tidyverse = as.character(packageVersion("tidyverse")),
    ggplot2 = as.character(packageVersion("ggplot2")),
    tikzDevice = as.character(packageVersion("tikzDevice")),
    patchwork = as.character(packageVersion("patchwork")),
    here = as.character(packageVersion("here")),
    rtracklayer = as.character(packageVersion("rtracklayer")),
    IRanges = as.character(packageVersion("IRanges")),
    jsonlite = as.character(packageVersion("jsonlite")),
    digest = as.character(packageVersion("digest")),
    reference_gff = basename(gff_path),
    reference_accession = "GCA_000182925.2",
    input_bigwigs = input_bigwigs,
    runs = purrr::transpose(runs),
    targets = purrr::transpose(targets),
    normalization = "CPM; inherited from the uniformly generated bigWig inputs",
    minimum_mapping_quality = 20L,
    profile_bin_bp = profile_bin_bp,
    window_flank_bp = window_flank_bp,
    plot_colors = as.list(plot_colors),
    plot_palette_reference = "https://okumuralab.org/~okumura/stat/colors.html",
    figure_width_in = figure_width_in,
    profile_heights_in = as.list(profile_heights_in),
    percentile_height_in = percentile_height_in,
    tikz_lwd_unit = 72.27 / 96,
    promoter_definition = "strand-aware TSS -1000 bp through +200 bp",
    gene_body_definition = "full NCBI GFF gene feature",
    ranking_population = "protein-coding genes on the seven NC12 nuclear chromosomes",
    percentile_method = "empirical midrank, calculated independently for each run, track variant, and region metric",
    exact_zero_treatment = "included as zero; counted separately in reporter_locus_signal_summary.tsv",
    missing_treatment = "excluded from regional means and counted separately; NA retained when a region has no covered bases",
    random_seed = NULL,
    pooling = "none",
    figure_format = "TikZ fragments generated from ggplot2 objects with tikzDevice"
  )

  jsonlite::write_json(
    payload,
    path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    na = "string"
  )
}

# Write SHA-256 checksums for all generated files
write_checksums <- function(output_dir) {
  checksum_path <- file.path(output_dir, "output_sha256.txt")
  output_files <- list.files(output_dir, full.names = TRUE) %>%
    sort() %>%
    purrr::keep(~ !file.info(.x)$isdir) %>%
    purrr::discard(~ basename(.x) == basename(checksum_path))
  digests <- output_files %>%
    purrr::map_chr(
      ~ digest::digest(.x, algo = "sha256", file = TRUE, serialize = FALSE)
    )

  writeLines(
    paste(digests, basename(output_files), sep = "  "),
    checksum_path
  )
}

# Run the complete reporter-locus analysis
run_analysis <- function(
  arguments,
  script_path,
  default_work_root,
  default_output_dir,
  targets,
  runs,
  variants,
  metrics,
  window_flank_bp,
  promoter_upstream_bp,
  promoter_downstream_bp,
  profile_bin_bp,
  plot_colors,
  figure_width_in,
  profile_heights_in,
  percentile_height_in
) {
  options(digits = 15)
  parsed_arguments <- parse_arguments(
    arguments,
    default_work_root,
    default_output_dir
  )
  work_root <- normalizePath(parsed_arguments$work_root, mustWork = TRUE)
  output_dir <- normalizePath(parsed_arguments$output_dir, mustWork = FALSE)
  gff_path <- file.path(
    work_root,
    "reference",
    "GCA_000182925.2_NC12_genomic.gff"
  )
  bigwig_root <- file.path(work_root, "bigwig")

  if (!file.exists(gff_path)) {
    stop("Missing reference GFF: ", gff_path, call. = FALSE)
  }
  if (!dir.exists(bigwig_root)) {
    stop("Missing bigWig directory: ", bigwig_root, call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Read genome annotation
  genes <- load_genes(gff_path)
  missing_targets <- setdiff(targets$locus_tag, genes$locus_tag)
  if (length(missing_targets) > 0L) {
    stop(
      "Target loci absent from GFF: ",
      paste(missing_targets, collapse = ", "),
      call. = FALSE
    )
  }

  # Quantify signal and prepare plotting data
  analysis_data <- collect_analysis_data(
    genes,
    bigwig_root,
    targets,
    runs,
    variants,
    metrics,
    window_flank_bp,
    promoter_upstream_bp,
    promoter_downstream_bp,
    profile_bin_bp
  )

  # Prepare figures
  marks <- unique(runs$mark)
  profile_plots <- names(variants) %>%
    purrr::set_names() %>%
    purrr::map(function(variant) {
      marks %>%
        purrr::set_names() %>%
        purrr::map(function(mark) {
          make_profile_plot(
            analysis_data$profiles,
            analysis_data$summary,
            genes,
            targets,
            runs,
            variant,
            mark,
            window_flank_bp,
            promoter_upstream_bp,
            promoter_downstream_bp,
            plot_colors
          )
        })
    })
  percentile_plots <- names(variants) %>%
    purrr::set_names() %>%
    purrr::map(
      ~ make_percentile_plot(
        analysis_data$summary,
        targets,
        runs,
        metrics,
        .x,
        plot_colors
      )
    )

  # Plot figures using the standard R graphics device
  plot_standard_figures(
    profile_plots,
    percentile_plots,
    figure_width_in,
    profile_heights_in
  )

  # =========================
  # Output results
  # =========================

  write_locus_table(
    file.path(output_dir, "reporter_loci.tsv"),
    genes,
    targets
  )
  write_signal_summary(
    file.path(output_dir, "reporter_locus_signal_summary.tsv"),
    analysis_data$summary
  )
  write_analysis_summary(
    file.path(output_dir, "analysis_summary.md"),
    analysis_data$summary,
    runs,
    metrics
  )

  obsolete_profile_files <- file.path(
    output_dir,
    paste0("reporter_locus_profiles.", names(variants), ".tex")
  )
  file.remove(obsolete_profile_files[file.exists(obsolete_profile_files)])

  for (variant in names(variants)) {
    for (mark in marks) {
      write_tikz_plot(
        profile_plots[[variant]][[mark]],
        file.path(
          output_dir,
          paste0("reporter_locus_profiles.", mark, ".", variant, ".tex")
        ),
        width = figure_width_in,
        height = profile_heights_in[[mark]]
      )
    }
    write_tikz_plot(
      percentile_plots[[variant]],
      file.path(
        output_dir,
        paste0("reporter_locus_percentiles.", variant, ".tex")
      ),
      width = figure_width_in,
      height = percentile_height_in
    )
  }

  write_parameters(
    file.path(output_dir, "analysis_parameters.json"),
    script_path,
    gff_path,
    runs,
    targets,
    variants,
    profile_bin_bp,
    window_flank_bp,
    plot_colors,
    figure_width_in,
    profile_heights_in,
    percentile_height_in
  )
  write_checksums(output_dir)

  # Print summaries to console
  cat(
    "Analyzed ", nrow(targets), " reporter loci across ", nrow(runs),
    " runs and ", length(variants), " track variants\n",
    sep = ""
  )
  cat(
    "Genome-wide ranking population: ", nrow(genes),
    " protein-coding genes\n",
    sep = ""
  )
  cat("Outputs: ", output_dir, "\n", sep = "")
}

# =========================
# Analysis settings
# =========================

# Define directories
script_path <- here::here(
  "05_public_H3K4",
  "scripts",
  "analyze_reporter_loci.R"
)
default_work_root <- Sys.getenv(
  "H3K4_WORK_ROOT",
  unset = "/Volumes/Garage/Re_analysis/260906_issue69_H3K4"
)
default_output_dir <- Sys.getenv(
  "H3K4_OUTPUT_DIR",
  unset = here::here("05_public_H3K4", "output", "reporter_loci")
)

# Define genomic regions and profile resolution
window_flank_bp <- 2000L
promoter_upstream_bp <- 1000L
promoter_downstream_bp <- 200L
profile_bin_bp <- 10L

# Define reporter loci
targets <- tibble::tribble(
  ~display_name, ~locus_tag, ~role,
  "pan-2", "NCU10048", "focal",
  "ad-3A", "NCU03166", "alternative",
  "ad-3B", "NCU03194", "alternative",
  "ad-8", "NCU09789", "alternative",
  "mtr", "NCU06619", "alternative",
  "his-3", "NCU03139", "alternative",
  "csr-1", "NCU00726", "exploratory"
)

# Define H3K4 ChIP-seq runs
runs <- tibble::tribble(
  ~sample_id, ~study, ~mark,
  "Ferraro2021_WT_H3K4me1", "Ferraro et al. 2021", "H3K4me1",
  "Ferraro2021_WT_H3K4me2", "Ferraro et al. 2021", "H3K4me2",
  "Sasaki2014_WT_H3K4me2", "Sasaki et al. 2014", "H3K4me2",
  "Ferraro2021_WT_H3K4me3", "Ferraro et al. 2021", "H3K4me3",
  "Storck2020_WT_H3K4me3", "Storck et al. 2020", "H3K4me3"
)

# Keep duplicate-retaining and nonduplicate tracks separate
variants <- c(
  nonduplicate = "q20.nonduplicate.cpm.bw",
  all_mapped = "q20.all_mapped.cpm.bw"
)
metrics <- c("promoter", "gene_body", "gene_body_plus_minus_2kb")

# Match the color and TikZ sizing conventions used by the other analyses
plot_colors <- c(
  blue = "#0068b7",
  midpoint = "#ffffff",
  orange = "#f39800"
)
figure_width_in <- 7.5
profile_heights_in <- c(
  H3K4me1 = 3.1,
  H3K4me2 = 4.5,
  H3K4me3 = 4.5
)
percentile_height_in <- 3.2

# =========================
# Run analysis
# =========================

runtime_arguments <- if (interactive()) {
  character()
} else {
  commandArgs(trailingOnly = TRUE)
}

run_analysis(
  runtime_arguments,
  script_path,
  default_work_root,
  default_output_dir,
  targets,
  runs,
  variants,
  metrics,
  window_flank_bp,
  promoter_upstream_bp,
  promoter_downstream_bp,
  profile_bin_bp,
  plot_colors,
  figure_width_in,
  profile_heights_in,
  percentile_height_in
)
