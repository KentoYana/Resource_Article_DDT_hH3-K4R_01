#!/usr/bin/env Rscript

# Quantify and plot H3K4 ChIP-seq signal at candidate reporter loci.
#
# The five selected runs remain separate. CPM values are never pooled across
# studies or compared across histone marks as if they were absolute quantities.

suppressPackageStartupMessages({
  library(IRanges)
  library(rtracklayer)
})

required_namespaces <- c("digest", "jsonlite")
missing_namespaces <- required_namespaces[
  !vapply(required_namespaces, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_namespaces) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_namespaces, collapse = ", "),
    call. = FALSE
  )
}

default_work_root <- "/Volumes/Garage/Re_analysis/260906_issue69_H3K4"
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Could not resolve the R script path", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_argument))
script_dir <- dirname(script_path)
default_output_dir <- file.path(dirname(script_dir), "output", "reporter_loci")

window_flank_bp <- 2000L
promoter_upstream_bp <- 1000L
promoter_downstream_bp <- 200L
profile_bin_bp <- 10L

targets <- data.frame(
  display_name = c("pan-2", "ad-3A", "ad-3B", "ad-8", "mtr", "his-3", "csr-1"),
  locus_tag = c(
    "NCU10048", "NCU03166", "NCU03194", "NCU09789", "NCU06619", "NCU03139", "NCU00726"
  ),
  role = c("focal", rep("alternative", 5L), "exploratory"),
  stringsAsFactors = FALSE
)

runs <- data.frame(
  sample_id = c(
    "Ferraro2021_WT_H3K4me1",
    "Ferraro2021_WT_H3K4me2",
    "Sasaki2014_WT_H3K4me2",
    "Ferraro2021_WT_H3K4me3",
    "Storck2020_WT_H3K4me3"
  ),
  study = c(
    "Ferraro et al. 2021",
    "Ferraro et al. 2021",
    "Sasaki et al. 2014",
    "Ferraro et al. 2021",
    "Storck et al. 2020"
  ),
  mark = c("H3K4me1", "H3K4me2", "H3K4me2", "H3K4me3", "H3K4me3"),
  stringsAsFactors = FALSE
)

variants <- c(
  nonduplicate = "q20.nonduplicate.cpm.bw",
  all_mapped = "q20.all_mapped.cpm.bw"
)

metrics <- c("promoter", "gene_body", "gene_body_plus_minus_2kb")

parse_arguments <- function(arguments) {
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
    if (!argument %in% c("--work-root", "--output-dir") || index == length(arguments)) {
      stop("Unknown or incomplete argument: ", argument, call. = FALSE)
    }
    key <- sub("^--", "", argument)
    key <- sub("-", "_", key, fixed = TRUE)
    result[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }
  result
}

extract_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)", key, "=([^;]+)")
  matches <- regexec(pattern, attributes, perl = TRUE)
  vapply(
    regmatches(attributes, matches),
    function(value) if (length(value) >= 2L) value[[2L]] else NA_character_,
    character(1)
  )
}

load_genes <- function(gff_path) {
  gff <- read.delim(
    gff_path,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "#",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  colnames(gff) <- c(
    "chrom", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"
  )
  locus_tag <- extract_attribute(gff$attributes, "locus_tag")
  gene_name <- extract_attribute(gff$attributes, "gene")
  feature_name <- extract_attribute(gff$attributes, "Name")
  gene_biotype <- extract_attribute(gff$attributes, "gene_biotype")
  keep <-
    gff$feature == "gene" &
    gene_biotype == "protein_coding" &
    grepl("^CM002(23[6-9]|24[0-2])\\.1$", gff$chrom) &
    !is.na(locus_tag)
  genes <- data.frame(
    locus_tag = locus_tag[keep],
    symbol = ifelse(
      !is.na(gene_name[keep]),
      gene_name[keep],
      ifelse(!is.na(feature_name[keep]), feature_name[keep], locus_tag[keep])
    ),
    chrom = gff$chrom[keep],
    start = as.integer(gff$start[keep]),
    end = as.integer(gff$end[keep]),
    strand = gff$strand[keep],
    stringsAsFactors = FALSE
  )
  genes$length_bp <- genes$end - genes$start + 1L
  if (anyDuplicated(genes$locus_tag)) {
    stop("Duplicate protein-coding gene features in GFF", call. = FALSE)
  }
  rownames(genes) <- genes$locus_tag
  genes
}

metric_regions <- function(genes, metric, chromosome_lengths) {
  if (metric == "gene_body") {
    region_start <- genes$start
    region_end <- genes$end
  } else if (metric == "gene_body_plus_minus_2kb") {
    region_start <- genes$start - window_flank_bp
    region_end <- genes$end + window_flank_bp
  } else if (metric == "promoter") {
    region_start <- ifelse(
      genes$strand == "+",
      genes$start - promoter_upstream_bp,
      genes$end - promoter_downstream_bp + 1L
    )
    region_end <- ifelse(
      genes$strand == "+",
      genes$start + promoter_downstream_bp - 1L,
      genes$end + promoter_upstream_bp
    )
  } else {
    stop("Unknown metric: ", metric, call. = FALSE)
  }
  region_start <- pmax(1L, as.integer(region_start))
  region_end <- pmin(
    as.integer(chromosome_lengths[genes$chrom]),
    as.integer(region_end)
  )
  data.frame(
    locus_tag = genes$locus_tag,
    chrom = genes$chrom,
    start = region_start,
    end = region_end,
    stringsAsFactors = FALSE
  )
}

view_statistic <- function(signal, regions, statistic) {
  result <- rep(NA_real_, nrow(regions))
  for (chromosome in unique(regions$chrom)) {
    indices <- which(regions$chrom == chromosome)
    views <- Views(
      signal[[chromosome]],
      start = regions$start[indices],
      end = regions$end[indices]
    )
    result[indices] <- statistic(views)
  }
  result
}

region_means <- function(signal, regions) {
  view_statistic(signal, regions, function(views) viewMeans(views, na.rm = TRUE))
}

empirical_midrank <- function(value, population) {
  finite <- population[is.finite(population)]
  if (!is.finite(value) || length(finite) == 0L) {
    return(c(percentile = NA_real_, rank_desc = NA_real_, population_n = length(finite)))
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

region_details <- function(signal, region) {
  values <- as.numeric(signal[[region$chrom]][region$start:region$end])
  valid <- is.finite(values)
  valid_values <- values[valid]
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

collect_summary <- function(genes, bigwig_root) {
  rows <- list()
  row_index <- 1L
  target_genes <- genes[targets$locus_tag, , drop = FALSE]

  for (variant in names(variants)) {
    for (run_index in seq_len(nrow(runs))) {
      run <- runs[run_index, , drop = FALSE]
      bigwig_name <- paste0(run$sample_id, ".", variants[[variant]])
      bigwig_path <- file.path(bigwig_root, bigwig_name)
      if (!file.exists(bigwig_path)) {
        stop("Missing bigWig: ", bigwig_path, call. = FALSE)
      }
      signal <- import(bigwig_path, as = "RleList")
      chromosome_lengths <- lengths(signal)
      if (!all(unique(genes$chrom) %in% names(signal))) {
        stop("A nuclear chromosome is absent from ", bigwig_name, call. = FALSE)
      }

      genome_values <- list()
      target_values <- list()
      target_regions <- list()
      for (metric in metrics) {
        regions <- metric_regions(genes, metric, chromosome_lengths)
        values <- region_means(signal, regions)
        names(values) <- regions$locus_tag
        genome_values[[metric]] <- values
        target_values[[metric]] <- values[targets$locus_tag]
        target_regions[[metric]] <- metric_regions(
          target_genes, metric, chromosome_lengths
        )
        rownames(target_regions[[metric]]) <- targets$locus_tag
      }

      for (target_index in seq_len(nrow(targets))) {
        target <- targets[target_index, , drop = FALSE]
        gene <- target_genes[target$locus_tag, , drop = FALSE]
        for (metric in metrics) {
          value <- target_values[[metric]][[target$locus_tag]]
          genome_rank <- empirical_midrank(value, genome_values[[metric]])
          reporter_rank <- empirical_midrank(value, target_values[[metric]])
          details <- region_details(
            signal,
            target_regions[[metric]][target$locus_tag, , drop = FALSE]
          )
          rows[[row_index]] <- data.frame(
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
            reporter_set_n = reporter_rank[["population_n"]],
            stringsAsFactors = FALSE
          )
          row_index <- row_index + 1L
        }
      }
      rm(signal)
      invisible(gc())
    }
  }
  do.call(rbind, rows)
}

summary_lookup <- function(summary, sample_id, variant, locus_tag, metric) {
  selected <-
    summary$sample_id == sample_id &
    summary$track_variant == variant &
    summary$locus_tag == locus_tag &
    summary$metric == metric
  if (sum(selected) != 1L) {
    stop("Summary lookup did not return exactly one row", call. = FALSE)
  }
  summary[selected, , drop = FALSE]
}

profile_for_gene <- function(signal, gene) {
  region_start <- max(1L, gene$start - window_flank_bp)
  region_end <- min(length(signal[[gene$chrom]]), gene$end + window_flank_bp)
  positions <- seq.int(region_start, region_end)
  values <- as.numeric(signal[[gene$chrom]][region_start:region_end])
  if (gene$strand == "+") {
    relative_positions <- positions - gene$start
  } else {
    values <- rev(values)
    relative_positions <- rev(gene$end - positions)
  }
  groups <- ceiling(seq_along(values) / profile_bin_bp)
  x <- as.numeric(tapply(relative_positions, groups, mean))
  y <- as.numeric(tapply(values, groups, mean, na.rm = TRUE))
  y[!is.finite(y)] <- 0
  data.frame(x = x, y = y)
}

tex_escape <- function(text) {
  text <- gsub("\\\\", "\\\\textbackslash{}", text)
  text <- gsub("([%&#_$])", "\\\\\\1", text, perl = TRUE)
  text
}

tex_number <- function(value, digits = 5L) {
  vapply(value, function(item) {
    if (!is.finite(item)) {
      return("0")
    }
    formatted <- formatC(item, format = "f", digits = digits)
    formatted <- sub("0+$", "", formatted)
    formatted <- sub("\\.$", "", formatted)
    if (formatted == "-0") "0" else formatted
  }, character(1), USE.NAMES = FALSE)
}

write_profile_tikz <- function(
  path, variant, genes, bigwig_root, summary
) {
  target_genes <- genes[targets$locus_tag, , drop = FALSE]
  profiles <- list()
  row_maxima <- numeric(nrow(runs))
  names(row_maxima) <- runs$sample_id

  for (run_index in seq_len(nrow(runs))) {
    run <- runs[run_index, , drop = FALSE]
    bigwig_path <- file.path(
      bigwig_root,
      paste0(run$sample_id, ".", variants[[variant]])
    )
    signal <- import(bigwig_path, as = "RleList")
    maximum <- 0
    for (target_index in seq_len(nrow(targets))) {
      target <- targets[target_index, , drop = FALSE]
      gene <- target_genes[target$locus_tag, , drop = FALSE]
      key <- paste(run$sample_id, target$locus_tag, sep = "|")
      profile <- profile_for_gene(signal, gene)
      profiles[[key]] <- profile
      maximum <- max(maximum, profile$y, na.rm = TRUE)
    }
    row_maxima[[run$sample_id]] <- max(1, maximum * 1.03)
    rm(signal)
    invisible(gc())
  }

  connection <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(
    c(
      "% Generated by analyze_reporter_loci.R. This file is a TikZ fragment.",
      "\\begin{tikzpicture}",
      "\\definecolor{signalblue}{HTML}{1546A0}",
      "\\pgfplotsset{compat=1.18}",
      "\\begin{groupplot}[",
      "  group style={group size=7 by 5, horizontal sep=0.42cm, vertical sep=0.55cm},",
      "  width=3.15cm, height=2.25cm, scale only axis,",
      "  axis line style={black!65, line width=0.3pt},",
      "  tick align=outside, tick style={black!65, line width=0.3pt},",
      "  tick label style={font=\\tiny}, title style={font=\\scriptsize, align=center},",
      "  label style={font=\\scriptsize}, clip=true",
      "]"
    ),
    connection
  )

  for (run_index in seq_len(nrow(runs))) {
    run <- runs[run_index, , drop = FALSE]
    y_max <- row_maxima[[run$sample_id]]
    for (target_index in seq_len(nrow(targets))) {
      target <- targets[target_index, , drop = FALSE]
      gene <- target_genes[target$locus_tag, , drop = FALSE]
      body <- summary_lookup(
        summary, run$sample_id, variant, target$locus_tag, "gene_body"
      )
      promoter <- summary_lookup(
        summary, run$sample_id, variant, target$locus_tag, "promoter"
      )
      options <- c(
        paste0("xmin=-", window_flank_bp),
        paste0("xmax=", gene$length_bp + window_flank_bp),
        "ymin=0",
        paste0("ymax=", tex_number(y_max, 3L)),
        paste0(
          "xtick={-", window_flank_bp, ",0,", gene$length_bp, ",",
          gene$length_bp + window_flank_bp, "}"
        ),
        "xticklabels={-2 kb,TSS,TES,+2 kb}"
      )
      if (run_index == 1L) {
        options <- c(
          options,
          paste0(
            "title={\\shortstack{", tex_escape(target$display_name), "\\\\",
            tex_escape(target$locus_tag), " (", gene$strand, ")}}"
          )
        )
      }
      if (target_index == 1L) {
        options <- c(
          options,
          paste0("ytick={0,", tex_number(y_max, 2L), "}"),
          paste0(
            "ylabel={\\shortstack{", tex_escape(run$study), "\\\\",
            tex_escape(run$mark), "\\\\CPM}}"
          )
        )
      } else {
        options <- c(options, "ytick=\\empty")
      }
      writeLines(
        paste0("\\nextgroupplot[", paste(options, collapse = ","), "]"),
        connection
      )
      writeLines(
        c(
          paste0(
            "\\path[fill=orange!14,draw=none] (axis cs:-", promoter_upstream_bp,
            ",0) rectangle (axis cs:", promoter_downstream_bp, ",",
            tex_number(y_max, 3L), ");"
          ),
          paste0(
            "\\path[fill=black!8,draw=none] (axis cs:0,0) rectangle (axis cs:",
            gene$length_bp, ",", tex_number(y_max, 3L), ");"
          ),
          "\\draw[black!55,dashed,line width=0.3pt] (axis cs:0,0) -- (axis cs:0,\\pgfkeysvalueof{/pgfplots/ymax});"
        ),
        connection
      )
      key <- paste(run$sample_id, target$locus_tag, sep = "|")
      profile <- profiles[[key]]
      coordinates <- paste0(
        "(", tex_number(profile$x, 2L), ",", tex_number(profile$y, 5L), ")"
      )
      area_coordinates <- c(
        paste0("(", tex_number(profile$x[[1L]], 2L), ",0)"),
        coordinates,
        paste0("(", tex_number(tail(profile$x, 1L), 2L), ",0)")
      )
      writeLines(
        paste0(
          "\\addplot[signalblue,line width=0.35pt,fill=signalblue!45,fill opacity=0.72] coordinates {",
          paste(area_coordinates, collapse = " "), "};"
        ),
        connection
      )
      writeLines(
        paste0(
          "\\node[anchor=north east,font=\\tiny,text=black!75] at (rel axis cs:0.98,0.97) {body pct ",
          round(body$genome_percentile_midrank), "; promoter pct ",
          round(promoter$genome_percentile_midrank), "};"
        ),
        connection
      )
    }
  }
  title_variant <- gsub("_", " ", variant, fixed = TRUE)
  writeLines(
    c(
      "\\end{groupplot}",
      paste0(
        "\\node[font=\\large\\bfseries,anchor=south] at ([yshift=0.85cm]group c4r1.north) ",
        "{Candidate reporter loci: gene-oriented H3K4me profiles (",
        tex_escape(title_variant), ")};"
      ),
      "\\node[font=\\scriptsize,align=center,text width=22cm,anchor=north] at ([yshift=-0.75cm]group c4r5.south) {Each row has one shared y-axis across all seven loci; y-axes are not shared between runs. Orange: promoter (-1 kb to +200 bp); gray: annotated gene body. Annotations are within-run genome-wide midrank percentiles.};",
      "\\end{tikzpicture}"
    ),
    connection
  )
}

heatmap_color <- function(percentile) {
  if (percentile <= 50) {
    paste0("heatlow!", round(100 - 2 * percentile), "!heatmid")
  } else {
    paste0("heatmid!", round(200 - 2 * percentile), "!heathigh")
  }
}

write_percentile_tikz <- function(path, variant, summary) {
  connection <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  title_variant <- gsub("_", " ", variant, fixed = TRUE)
  panel_titles <- c(
    promoter = "Promoter (-1 kb/+200 bp)",
    gene_body = "Gene body",
    gene_body_plus_minus_2kb = "Gene body $\\pm$ 2 kb"
  )
  panel_offsets <- c(0, 9.5, 19)
  writeLines(
    c(
      "% Generated by analyze_reporter_loci.R. This file is a TikZ fragment.",
      "\\begin{tikzpicture}[x=0.78cm,y=0.78cm]",
      "\\definecolor{heatlow}{HTML}{440154}",
      "\\definecolor{heatmid}{HTML}{21918C}",
      "\\definecolor{heathigh}{HTML}{FDE725}",
      paste0(
        "\\node[font=\\large\\bfseries] at (13,-0.1) {Candidate reporter H3K4me signal percentiles (",
        tex_escape(title_variant), ")};"
      )
    ),
    connection
  )
  for (metric_index in seq_along(metrics)) {
    metric <- metrics[[metric_index]]
    offset <- panel_offsets[[metric_index]]
    writeLines(
      paste0(
        "\\node[font=\\normalsize\\bfseries] at (", offset + 3.5,
        ",-1.25) {", panel_titles[[metric]], "};"
      ),
      connection
    )
    for (run_index in seq_len(nrow(runs))) {
      run <- runs[run_index, , drop = FALSE]
      for (target_index in seq_len(nrow(targets))) {
        target <- targets[target_index, , drop = FALSE]
        row <- summary_lookup(
          summary, run$sample_id, variant, target$locus_tag, metric
        )
        percentile <- row$genome_percentile_midrank
        x_left <- offset + target_index - 1L
        y_top <- -1.7 - (run_index - 1L)
        text_color <- if (percentile < 45) "white" else "black"
        writeLines(
          c(
            paste0(
              "\\fill[", heatmap_color(percentile), "] (", x_left, ",", y_top,
              ") rectangle ++(1,-1);"
            ),
            paste0(
              "\\node[font=\\scriptsize,text=", text_color, "] at (", x_left + 0.5,
              ",", y_top - 0.5, ") {", round(percentile), "};"
            )
          ),
          connection
        )
      }
      if (metric_index == 1L) {
        y_center <- -2.2 - (run_index - 1L)
        writeLines(
          paste0(
            "\\node[anchor=east,font=\\scriptsize,align=right] at (-0.2,", y_center,
            ") {", tex_escape(run$study), "\\\\", tex_escape(run$mark), "};"
          ),
          connection
        )
      }
    }
    for (target_index in seq_len(nrow(targets))) {
      x_center <- offset + target_index - 0.5
      writeLines(
        paste0(
          "\\node[anchor=east,rotate=45,font=\\scriptsize] at (", x_center,
          ",-6.95) {", tex_escape(targets$display_name[[target_index]]), "};"
        ),
        connection
      )
    }
    writeLines(
      paste0(
        "\\draw[black,line width=0.4pt] (", offset, ",-1.7) rectangle (",
        offset + 7, ",-6.7);"
      ),
      connection
    )
  }
  legend_x <- 27.4
  for (legend_index in 0:19) {
    percentile <- legend_index * 5
    y_bottom <- -6.7 + legend_index * 0.25
    writeLines(
      paste0(
        "\\fill[", heatmap_color(percentile), "] (", legend_x, ",", y_bottom,
        ") rectangle ++(0.35,0.25);"
      ),
      connection
    )
  }
  writeLines(
    c(
      paste0("\\draw[black,line width=0.4pt] (", legend_x, ",-6.7) rectangle ++(0.35,5);"),
      paste0("\\node[anchor=west,font=\\scriptsize] at (", legend_x + 0.45, ",-6.7) {0};"),
      paste0("\\node[anchor=west,font=\\scriptsize] at (", legend_x + 0.45, ",-4.2) {50};"),
      paste0("\\node[anchor=west,font=\\scriptsize] at (", legend_x + 0.45, ",-1.7) {100};"),
      paste0(
        "\\node[rotate=90,font=\\scriptsize] at (", legend_x + 1.35,
        ",-4.2) {Within-run genome-wide percentile (midrank)};"
      ),
      "\\end{tikzpicture}"
    ),
    connection
  )
}

write_locus_table <- function(path, genes) {
  target_genes <- genes[targets$locus_tag, , drop = FALSE]
  table <- data.frame(
    display_name = targets$display_name,
    locus_tag = targets$locus_tag,
    role = targets$role,
    chrom = target_genes$chrom,
    start_1based = target_genes$start,
    end_1based = target_genes$end,
    strand = target_genes$strand,
    length_bp = target_genes$length_bp,
    gff_gene_name = target_genes$symbol,
    stringsAsFactors = FALSE
  )
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_analysis_summary <- function(path, summary) {
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
    run <- runs[run_index, , drop = FALSE]
    for (metric in metrics) {
      row <- summary_lookup(
        summary, run$sample_id, "nonduplicate", "NCU10048", metric
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

write_signal_summary <- function(path, summary) {
  formatted <- summary
  six_decimal_columns <- c(
    "mean_cpm",
    "max_cpm",
    "nonzero_fraction",
    "genome_percentile_midrank",
    "genome_rank_desc",
    "reporter_set_rank_desc"
  )
  for (column in six_decimal_columns) {
    formatted[[column]] <- ifelse(
      is.na(formatted[[column]]),
      "NA",
      formatC(formatted[[column]], format = "f", digits = 6L)
    )
  }
  write.table(
    formatted,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

write_parameters <- function(path, gff_path) {
  input_bigwigs <- unlist(lapply(
    names(variants),
    function(variant) paste0(runs$sample_id, ".", variants[[variant]])
  ))
  payload <- list(
    analysis_script = basename(script_path),
    R = as.character(getRversion()),
    rtracklayer = as.character(packageVersion("rtracklayer")),
    IRanges = as.character(packageVersion("IRanges")),
    jsonlite = as.character(packageVersion("jsonlite")),
    digest = as.character(packageVersion("digest")),
    reference_gff = basename(gff_path),
    reference_accession = "GCA_000182925.2",
    input_bigwigs = input_bigwigs,
    runs = lapply(seq_len(nrow(runs)), function(index) {
      as.list(runs[index, , drop = FALSE])
    }),
    targets = lapply(seq_len(nrow(targets)), function(index) {
      as.list(targets[index, , drop = FALSE])
    }),
    normalization = "CPM; inherited from the uniformly generated bigWig inputs",
    minimum_mapping_quality = 20L,
    profile_bin_bp = profile_bin_bp,
    window_flank_bp = window_flank_bp,
    promoter_definition = "strand-aware TSS -1000 bp through +200 bp",
    gene_body_definition = "full NCBI GFF gene feature",
    ranking_population = "protein-coding genes on the seven NC12 nuclear chromosomes",
    percentile_method = "empirical midrank, calculated independently for each run, track variant, and region metric",
    exact_zero_treatment = "included as zero; counted separately in reporter_locus_signal_summary.tsv",
    missing_treatment = "excluded from regional means and counted separately; NA retained when a region has no covered bases",
    random_seed = NULL,
    pooling = "none",
    figure_format = "TikZ fragments containing tikzpicture environments"
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

write_checksums <- function(output_dir) {
  checksum_path <- file.path(output_dir, "output_sha256.txt")
  output_files <- sort(list.files(output_dir, full.names = TRUE))
  output_files <- output_files[
    file.info(output_files)$isdir == FALSE & basename(output_files) != basename(checksum_path)
  ]
  digests <- vapply(
    output_files,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
  writeLines(paste(digests, basename(output_files), sep = "  "), checksum_path)
}

main <- function() {
  options(digits = 15)
  arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
  work_root <- normalizePath(arguments$work_root, mustWork = TRUE)
  output_dir <- normalizePath(
    arguments$output_dir,
    mustWork = FALSE
  )
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

  genes <- load_genes(gff_path)
  missing_targets <- setdiff(targets$locus_tag, genes$locus_tag)
  if (length(missing_targets) > 0L) {
    stop(
      "Target loci absent from GFF: ",
      paste(missing_targets, collapse = ", "),
      call. = FALSE
    )
  }

  write_locus_table(file.path(output_dir, "reporter_loci.tsv"), genes)
  summary <- collect_summary(genes, bigwig_root)
  write_signal_summary(
    file.path(output_dir, "reporter_locus_signal_summary.tsv"),
    summary
  )
  for (variant in names(variants)) {
    write_profile_tikz(
      file.path(output_dir, paste0("reporter_locus_profiles.", variant, ".tex")),
      variant,
      genes,
      bigwig_root,
      summary
    )
    write_percentile_tikz(
      file.path(output_dir, paste0("reporter_locus_percentiles.", variant, ".tex")),
      variant,
      summary
    )
  }
  write_analysis_summary(
    file.path(output_dir, "analysis_summary.md"),
    summary
  )
  write_parameters(
    file.path(output_dir, "analysis_parameters.json"),
    gff_path
  )
  write_checksums(output_dir)
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

main()
