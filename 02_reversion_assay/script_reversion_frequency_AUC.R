# script_reversion_frequency_AUC.R
# Kento Yanagisawa
# This script analyzes reversion frequency using baseline-corrected AUC.

# Make these packages and their associated functions
# available to use in this script
library("tikzDevice")
library('RColorBrewer')
library("tidyverse")
library("ggfortify")
library("patchwork")
library("here")
library("ggh4x")

# Clear R's brain
rm(list = ls())

# Set project root
here::i_am("02_reversion_assay/script_reversion_frequency_AUC.R")

# Define directories
script_dir  <- here("02_reversion_assay")
dataset_dir <- here("02_reversion_assay", "dataset")
output_dir  <- here("02_reversion_assay", "output")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read raw data
raw_csv <- read.csv(here("02_reversion_assay", "dataset", "reversion_data.csv"))

raw_csv <- raw_csv %>%
  dplyr::filter(!is.na(Survival), !is.na(Revertant))

# Read strain information
strain_list <- read.csv(here("02_reversion_assay", "dataset", "strain_list.csv"))

genotype_list <- read.csv(here("02_reversion_assay", "dataset", "genotype_list.csv"))

raw_csv <- raw_csv %>%
  dplyr::left_join(strain_list, by = "Strain")

raw_csv$genotype <- factor(raw_csv$genotype, levels = genotype_list$genotype)

# Summarize technical replicates within each Date × Strain × UV_dose
expData_raw <- raw_csv %>%
  dplyr::group_by(Date, Strain, sibling, allele, genotype, UV_dose) %>%
  dplyr::summarise(
    Suv_samples = sum(!is.na(Survival)),
    Suv_ave     = mean(Survival, na.rm = TRUE),
    Suv_se      = sd(Survival, na.rm = TRUE) / sqrt(sum(!is.na(Survival))),
    Rev_samples = sum(!is.na(Revertant)),
    Rev_ave     = mean(Revertant, na.rm = TRUE),
    Rev_se      = sd(Revertant, na.rm = TRUE) / sqrt(sum(!is.na(Revertant))),
    .groups = "drop"
  )

# Calculate reversion frequency
propagate_error_division <- function(A, A_err, B, B_err) {
  sqrt((A_err / A)^2 + (B_err / B)^2) * (A / B) * (10^7 / 10^4)
}

expData <- expData_raw %>%
  dplyr::mutate(
    Rev_freq_raw = (Rev_ave / (Suv_ave * 10^4)) * 10^7,
    Rev_freq_se  = propagate_error_division(
      Rev_ave, Rev_se,
      Suv_ave, Suv_se
    )
  )

# Sanity check: confirm one row per Date × Strain × UV_dose
dup_check <- expData %>%
  dplyr::count(Date, Strain, sibling, allele, genotype, UV_dose, name = "n_rows") %>%
  dplyr::filter(n_rows != 1)

if (nrow(dup_check) > 0) {
  stop(
    "Duplicated grouping structure detected in Date × Strain × UV_dose.\n",
    "Please inspect object 'dup_check'."
  )
}

# 6. Plotting data for dose-response curves mean of biological replicates across Date
plotData <- expData %>%
  dplyr::group_by(sibling, allele, genotype, UV_dose) %>%
  dplyr::summarise(
    Reversion_frequency    = mean(Rev_freq_raw, na.rm = TRUE),
    Reversion_frequency_se = sd(Rev_freq_raw, na.rm = TRUE) / sqrt(sum(!is.na(Rev_freq_raw))),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    seMax = Reversion_frequency + Reversion_frequency_se,
    seMin = Reversion_frequency - Reversion_frequency_se
  )

# Trapezoidal AUC function
calc_auc_trapz <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# Calculate baseline-corrected AUC for each Date × Strain
required_doses <- c(0, 60, 90, 120)

aucData <- expData %>%
  dplyr::group_by(Date, Strain, sibling, allele, genotype) %>%
  dplyr::filter(all(required_doses %in% UV_dose)) %>%
  dplyr::arrange(UV_dose, .by_group = TRUE) %>%
  dplyr::mutate(
    baseline_0J = Rev_freq_raw[UV_dose == 0][1],
    Rev_freq_bc = Rev_freq_raw - baseline_0J
  ) %>%
  dplyr::filter(UV_dose > 0) %>%
  dplyr::summarise(
    baseline_0J = first(baseline_0J),
    AUC_bc      = calc_auc_trapz(UV_dose, Rev_freq_bc),
    .groups = "drop"
  ) %>%
  dplyr::ungroup()

# Analyze each sibling pair without requiring Date matching
# Welch's t-test is used by default
analyze_auc_pair_unpaired <- function(auc_data, sibling_id) {

  pair_data <- auc_data %>%
    dplyr::filter(sibling == sibling_id)

  genotypes_in_pair <- unique(as.character(pair_data$genotype))

  if (!("wild-type" %in% genotypes_in_pair)) {
    stop(paste("wild-type is missing in sibling pair:", sibling_id))
  }

  mutant_name <- setdiff(genotypes_in_pair, "wild-type")

  if (length(mutant_name) != 1) {
    stop(paste("Sibling pair", sibling_id, "must contain exactly one mutant genotype."))
  }

  mutant_name <- mutant_name[1]

  pair_data$genotype <- droplevels(pair_data$genotype)

  wt_n  <- sum(pair_data$genotype == "wild-type")
  mut_n <- sum(pair_data$genotype == mutant_name)

  if (wt_n < 2 || mut_n < 2) {
    warning(paste("Too few observations for Welch's t-test in", sibling_id))
    test_res <- NULL
  } else {
    test_res <- t.test(AUC_bc ~ genotype, data = pair_data, var.equal = FALSE)
  }

  # Safe labels for pivoted column names
  pair_data_safe <- pair_data %>%
    dplyr::mutate(
      genotype_safe = dplyr::case_when(
        genotype == "wild-type" ~ "wild_type",
        genotype == "hH3-K4R"   ~ "hH3_K4R",
        genotype == "set-1 KO"  ~ "set_1_KO",
        TRUE ~ gsub("[^A-Za-z0-9]+", "_", as.character(genotype))
      )
    )

  mutant_safe <- unique(pair_data_safe$genotype_safe[pair_data_safe$genotype == mutant_name])

  summary_table <- pair_data_safe %>%
    dplyr::group_by(sibling, allele, genotype, genotype_safe) %>%
    dplyr::summarise(
      mean_AUC = mean(AUC_bc, na.rm = TRUE),
      se_AUC   = sd(AUC_bc, na.rm = TRUE) / sqrt(sum(!is.na(AUC_bc))),
      n        = sum(!is.na(AUC_bc)),
      .groups = "drop"
    )

  wide_summary <- summary_table %>%
    dplyr::select(sibling, allele, genotype_safe, mean_AUC, se_AUC, n) %>%
    tidyr::pivot_wider(
      names_from = genotype_safe,
      values_from = c(mean_AUC, se_AUC, n)
    )

  result_summary <- tibble(
    sibling    = sibling_id,
    allele     = unique(pair_data$allele)[1],
    mutant     = mutant_name,
    wt_n       = wt_n,
    mut_n      = mut_n,
    wt_mean    = wide_summary$mean_AUC_wild_type,
    mut_mean   = wide_summary[[paste0("mean_AUC_", mutant_safe)]],
    delta_mean = wide_summary[[paste0("mean_AUC_", mutant_safe)]] - wide_summary$mean_AUC_wild_type,
    p_value    = if (is.null(test_res)) NA_real_ else test_res$p.value
  ) %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        is.na(p_value)  ~ "n.d.",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "n.s."
      )
    )

  return(list(
    sibling        = sibling_id,
    mutant         = mutant_name,
    pair_data      = pair_data,
    summary_table  = summary_table,
    result_summary = result_summary,
    t_test         = test_res
  ))
}

sibling_ids <- unique(strain_list$sibling)

pair_results <- lapply(sibling_ids, function(x) analyze_auc_pair_unpaired(aucData, x))
names(pair_results) <- sibling_ids

pair_summary <- dplyr::bind_rows(lapply(pair_results, function(x) x$result_summary))

# Helper table for plotting AUC bars
aucPlotData_pair <- aucData %>%
  dplyr::group_by(sibling, allele, genotype) %>%
  dplyr::summarise(
    AUC_bc_mean = mean(AUC_bc, na.rm = TRUE),
    AUC_bc_se   = sd(AUC_bc, na.rm = TRUE) / sqrt(sum(!is.na(AUC_bc))),
    n           = sum(!is.na(AUC_bc)),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    pair_summary %>% dplyr::select(sibling, significance),
    by = "sibling"
  ) %>%
  dplyr::mutate(
    ymin = AUC_bc_mean - AUC_bc_se,
    ymax = AUC_bc_mean + AUC_bc_se
  )

# Plot dose-response curves in a faceted layout
plot_reversion_facet <- function(DataForPlot, pair_summary_df) {

  FilteredData <- DataForPlot %>%
    dplyr::filter(UV_dose > 0) %>%
    dplyr::mutate(
      genotype_label = dplyr::case_when(
        sibling %in% c("Op46", "Op47") ~ "hH3-K4R",
        sibling %in% c("Op55", "Op51") ~ "set-1 KO"
      ),
      allele = factor(allele, levels = c("B36", "OGW1")),
      genotype_label = factor(genotype_label, levels = c("hH3-K4R", "set-1 KO"))
    )

  pval_df <- pair_summary_df %>%
    dplyr::mutate(
      genotype_label = dplyr::case_when(
        sibling %in% c("Op46", "Op47") ~ "hH3-K4R",
        sibling %in% c("Op55", "Op51") ~ "set-1 KO"
      ),
      label = dplyr::case_when(
        is.na(p_value)  ~ "Welch's t-test: n.d.",
        p_value < 0.001 ~ "p < 0.001",
        TRUE            ~ paste0("p = ", signif(p_value, 3))
      ),
      allele = factor(allele, levels = c("B36", "OGW1")),
      genotype_label = factor(genotype_label, levels = c("hH3-K4R", "set-1 KO"))
    ) %>%
    dplyr::select(sibling, allele, genotype_label, label)

  anno_pos <- FilteredData %>%
    dplyr::group_by(sibling, allele, genotype_label) %>%
    dplyr::summarise(
      x = min(UV_dose, na.rm = TRUE),
      y = 90,
      .groups = "drop"
    ) %>%
    dplyr::left_join(pval_df, by = c("sibling", "allele", "genotype_label")) %>%
    dplyr::filter(!is.na(label))

  ggplot(
    FilteredData,
    aes(
      x = UV_dose,
      y = Reversion_frequency,
      color = genotype,
      shape = genotype,
      group = genotype
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = Reversion_frequency - Reversion_frequency_se,
        ymax = Reversion_frequency + Reversion_frequency_se
      ),
      width = 6,
      linewidth = 1
    ) +
    geom_text(
      data = anno_pos,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.5
    ) +
    scale_y_log10(
      limits = c(
        max(min(FilteredData$seMin[FilteredData$seMin > 0], na.rm = TRUE), 1e-3),
        max(FilteredData$seMax, na.rm = TRUE) + 25
      )
    ) +
    ggh4x::facet_nested(~ allele + genotype_label) +
    theme_bw(base_size = 10) +
    xlab("UV dose (J/m2)") +
    ylab("UV-induced reversion frequency\nat the pan-2 locus (pan+ / 10 7 survivors)") +
    scale_shape_manual(values = c(15, 0, 16, 1)) +
    scale_color_manual(values = c("#02010C", "#0068b7", "#f39800", "#009944")) +
    theme(
      axis.text = element_text(size = 10, colour = "black"),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 1.2),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(0.0, "cm"),
      legend.box.spacing = unit(0.0, "cm"),
      aspect.ratio = (1 + sqrt(5)) / 2
    )
}

# Plot 0 J data in a faceted layout
plot_reversion_0J_facet <- function(DataForPlot) {

  FilteredData <- DataForPlot %>%
    dplyr::filter(UV_dose == 0) %>%
    dplyr::mutate(
      genotype_label = dplyr::case_when(
        sibling %in% c("Op46", "Op47") ~ "hH3-K4R",
        sibling %in% c("Op55", "Op51") ~ "set-1 KO"
      ),
      allele = factor(allele, levels = c("B36", "OGW1")),
      genotype_label = factor(genotype_label, levels = c("hH3-K4R", "set-1 KO"))
    )

  ggplot(
    FilteredData,
    aes(
      x = genotype,
      y = Reversion_frequency,
      color = genotype,
      fill = genotype
    )
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = Reversion_frequency - Reversion_frequency_se,
        ymax = Reversion_frequency + Reversion_frequency_se
      ),
      width = 0.25,
      linewidth = 1
    ) +
    ggh4x::facet_nested(~ allele + genotype_label) +
    theme_bw(base_size = 10) +
    xlab("") +
    ylab("Spontaneous reversion frequency\nat the pan-2 locus (pan+ / 10 7 survivors)") +
    scale_color_manual(values = c("#02010C", "#0068b7", "#f39800", "#009944")) +
    scale_fill_manual(values = c("#02010C", "#0068b7", "#f39800", "#009944")) +
    theme(
      axis.text = element_text(size = 10, colour = "black"),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 1.2),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(0.0, "cm"),
      legend.box.spacing = unit(0.0, "cm"),
      aspect.ratio = (1 + sqrt(5)) / 2
    )
}

# 13. Plot AUC in a faceted layout
plot_auc_facet <- function(DataForPlot) {

  FilteredData <- DataForPlot %>%
    dplyr::mutate(
      genotype_label = dplyr::case_when(
        sibling %in% c("Op46", "Op47") ~ "hH3-K4R",
        sibling %in% c("Op55", "Op51") ~ "set-1 KO"
      ),
      allele = factor(allele, levels = c("B36", "OGW1")),
      genotype_label = factor(genotype_label, levels = c("hH3-K4R", "set-1 KO"))
    )

  anno_pos <- FilteredData %>%
    dplyr::group_by(sibling, allele, genotype_label) %>%
    dplyr::summarise(
      x = 1.5,
      y = max(ymax, na.rm = TRUE) + 0.05 * max(abs(AUC_bc_mean), na.rm = TRUE),
      label = unique(significance)[1],
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(label))

  ggplot(
    FilteredData,
    aes(
      x = genotype,
      y = AUC_bc_mean,
      color = genotype,
      fill = genotype
    )
  ) +
    geom_col(width = 0.7, alpha = 0.8) +
    geom_errorbar(
      aes(
        ymin = ymin,
        ymax = ymax
      ),
      width = 0.25,
      linewidth = 1
    ) +
    geom_text(
      data = anno_pos,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 5
    ) +
    ggh4x::facet_nested(~ allele + genotype_label) +
    theme_bw(base_size = 10) +
    xlab("") +
    ylab("Baseline-corrected AUC of reversion frequency\n(pan+ / 10 7 survivors × J/m 2)") +
    scale_color_manual(values = c("#02010C", "#0068b7", "#f39800", "#009944")) +
    scale_fill_manual(values = c("#02010C", "#0068b7", "#f39800", "#009944")) +
    theme(
      axis.text = element_text(size = 10, colour = "black"),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 1.2),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(0.0, "cm"),
      legend.box.spacing = unit(0.0, "cm"),
      aspect.ratio = (1 + sqrt(5)) / 2
    )
}

# Make faceted figure
fig_reversion <- plot_reversion_facet(plotData, pair_summary)
plot(fig_reversion)

# fig_0J <- plot_reversion_0J_facet(plotData)
# plot(fig_0J)
#
# fig_auc <- plot_auc_facet(aucPlotData_pair)
# plot(fig_auc)

# Output results

# Output processed data
write.csv(
  raw_csv,
  here("02_reversion_assay", "output", "raw_csv_filtered.csv"),
  row.names = FALSE
)

write.csv(
  expData_raw,
  here("02_reversion_assay", "output", "expData_raw.csv"),
  row.names = FALSE
)

write.csv(
  expData,
  here("02_reversion_assay", "output", "expData.csv"),
  row.names = FALSE
)

write.csv(
  plotData,
  here("02_reversion_assay", "output", "plotData.csv"),
  row.names = FALSE
)

write.csv(
  aucData,
  here("02_reversion_assay", "output", "aucData.csv"),
  row.names = FALSE
)

write.csv(
  pair_summary,
  here("02_reversion_assay", "output", "pair_summary.csv"),
  row.names = FALSE
)

write.csv(
  aucPlotData_pair,
  here("02_reversion_assay", "output", "aucPlotData_pair.csv"),
  row.names = FALSE
)

# Output pair-wise data and summary tables
for (sid in names(pair_results)) {

  write.csv(
    pair_results[[sid]]$pair_data,
    here("02_reversion_assay", "output", paste0("pair_data_", sid, ".csv")),
    row.names = FALSE
  )

  write.csv(
    pair_results[[sid]]$summary_table,
    here("02_reversion_assay", "output", paste0("pair_summary_table_", sid, ".csv")),
    row.names = FALSE
  )
}

# Output Welch's t-test results as text
sink(here("02_reversion_assay", "output", "welch_t_tests.txt"))

for (sid in names(pair_results)) {
  cat("========================================\n")
  cat("Sibling:", sid, "\n")
  cat("Mutant:", pair_results[[sid]]$mutant, "\n")
  cat("========================================\n")

  if (is.null(pair_results[[sid]]$t_test)) {
    cat("Welch's t-test was not performed because of too few observations.\n\n")
  } else {
    print(pair_results[[sid]]$t_test)
    cat("\n")
  }
}

sink()

# Output figure as TikZ
tikz(
  here("02_reversion_assay", "output", "reversion_frequency_AUC.tex"),
  width = 7.5,
  lwdUnit = 72.27 / 96
)
plot(fig_reversion)
dev.off()

# Optional: output additional figures as TikZ if generated
# fig_0J <- plot_reversion_0J_facet(plotData)
# tikz(
#   here("02_reversion_assay", "output", "reversion_frequency_0J.tex"),
#   width = 7.5,
#   lwdUnit = 72.27 / 96
# )
# plot(fig_0J)
# dev.off()
#
# fig_auc <- plot_auc_facet(aucPlotData_pair)
# tikz(
#   here("02_reversion_assay", "output", "reversion_frequency_AUC_bar.tex"),
#   width = 7.5,
#   lwdUnit = 72.27 / 96
# )
# plot(fig_auc)
# dev.off()

# Print summaries to console
aucData
pair_summary

pair_results$Op47$t_test
pair_results$Op51$t_test
pair_results$Op46$t_test
pair_results$Op55$t_test

pair_results$Op47$pair_data
pair_results$Op51$pair_data
pair_results$Op46$pair_data
pair_results$Op55$pair_data