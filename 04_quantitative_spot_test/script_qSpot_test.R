# script_qSpot_test.R
# Kento Yanagisawa
# This script is for the analysis of semi-quantitative spot-test(4*6 matrix)

# make these packages and their associated functions
# available to use in this script
library("tikzDevice")
library('RColorBrewer')
library("tidyverse")
library("ggfortify")
library("patchwork")
library("here")
library('betareg')
library('emmeans')

# clear R's brain
rm(list = ls())
par(mfrow = c(2, 2), ask = FALSE)

# Set project root
here::i_am("04_quantitative_spot_test/script_qSpot_test.R")

# Define directories
script_dir  <- here("04_quantitative_spot_test")
dataset_dir <- here("04_quantitative_spot_test", "dataset")
output_dir  <- here("04_quantitative_spot_test", "output")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to perform likelihood ratio test
likelihood_ratio_test <- function(model1, model2) {
  # Get log-likelihood for each model
  logLik_model1 <- logLik(model1)
  logLik_model2 <- logLik(model2)
  # Calculate the likelihood ratio statistic
  LR_statistic <- 2 * (logLik_model1 - logLik_model2)
  # Calculate the difference in degrees of freedom
  df_diff <- attr(logLik_model1, "df") - attr(logLik_model2, "df")
  # Calculate the p-value
  p_value <- pchisq(LR_statistic, df = df_diff, lower.tail = FALSE)
  # Print the results
  cat("Likelihood Ratio Statistic:", LR_statistic, "\n")
  cat("Degrees of Freedom Difference:", df_diff, "\n")
  cat("p-value:", p_value, "\n")
  # Return the results as a list
  return(list(
    LR_statistic = as.numeric(LR_statistic),
    df_diff = as.numeric(df_diff),
    p_value = as.numeric(p_value)
  ))
}

# palette("Okabe-Ito")
colsBlack = brewer.pal(6, "Greys")
colsOrange = brewer.pal(6, "Oranges")
colsBlue = brewer.pal(6, "Blues")
colsPurple = brewer.pal(5, "Purples")

# Read experiment list
exp_list <- read_csv(here("04_quantitative_spot_test", "dataset", "exp_list.csv"))

# Function to analyze one target
analyze_qspot_target <- function(target_info) {

  target_name <- target_info$target
  target_dose <- target_info$dose

  cat("\n")
  cat("========================================\n")
  cat("Analyzing target:", target_name, "\n")
  cat("Target suspension:", target_dose, "\n")
  cat("========================================\n")

  # Define target-specific output directory
  target_output_dir <- here(
    "04_quantitative_spot_test",
    "output",
    target_name
  )

  if (!dir.exists(target_output_dir)) {
    dir.create(target_output_dir, recursive = TRUE)
  }

  # set analyzing data set
  spotpath <- paste('qSpot-', target_name, sep = '')

  # read dataset
  csv_list <- list.files(
    path = here("04_quantitative_spot_test", "dataset", spotpath),
    recursive = TRUE,
    pattern = '.*_results.csv',
    full.names = TRUE
  )

  if (length(csv_list) == 0) {
    warning(paste("No result CSV files found for target:", target_name))
    return(NULL)
  }

  raw_csv <- read_csv(csv_list, id = 'csv_path')

  strain_list <- read_csv(
    here("04_quantitative_spot_test", "dataset", spotpath, "strain_list.csv")
  )

  raw_csv <- raw_csv %>%
    mutate(exp_ID = str_extract(csv_path, '/results_.*/UV_')) %>%
    mutate(exp_ID = str_remove_all(exp_ID, '/results_|/UV_')) %>%
    mutate(csv_path = str_remove_all(csv_path, '.*data_availability')) %>%
    separate(Image, c('exp_condition', 'suspension', 'dose', 'unit', 'photo_ID'), sep = '_') %>%
    mutate(dose = as.numeric(dose)) %>%
    mutate(suspension = as.numeric(suspension)) %>%
    mutate(Conidia = 2 * 10^suspension / 5^(column - 1)) %>%
    left_join(strain_list, by = 'row') %>%
    filter(suspension == target_dose) %>%
    na.omit()

  raw_csv$strain <- factor(raw_csv$strain, levels = strain_list$strain)

  # calculate mean of each duplicate spot
  raw_csv_summarize <- raw_csv %>%
    group_by(suspension, dose, strain, column, exp_ID, Conidia) %>%
    summarize(
      Colony_mean = mean(Colony, na.rm = TRUE),
      .groups = "drop"
    )

  raw_csv_summarize <- mutate(
    raw_csv_summarize,
    Colony_mean = round(Colony_mean, digit = 2)
  )

  # get control value
  raw_csv_control_raw <- filter(raw_csv_summarize, dose == 0)

  raw_csv_control <- data.frame(
    suspension = raw_csv_control_raw$suspension,
    strain = raw_csv_control_raw$strain,
    column = raw_csv_control_raw$column,
    exp_ID = raw_csv_control_raw$exp_ID,
    Colony_control = raw_csv_control_raw$Colony_mean,
    Conidia_control = max(raw_csv_control_raw$Conidia)
  )

  raw_csv_summarize <- left_join(
    raw_csv_summarize,
    raw_csv_control,
    by = c("suspension", "strain", "column", "exp_ID")
  )

  # A zero UV-control mean cannot define a normalized spot-coverage ratio.
  # Retain the excluded observations for auditing, but remove them before
  # calculating ratios and the sample size for the endpoint adjustment.
  excluded_zero_control <- raw_csv_summarize %>%
    filter(Colony_control == 0)

  if (any(!is.finite(raw_csv_summarize$Colony_control)) ||
      any(raw_csv_summarize$Colony_control < 0)) {
    stop("Missing, non-finite, or negative UV-control coverage for ", target_name)
  }

  ratio_data <- raw_csv_summarize %>%
    filter(Colony_control != 0)

  # calculate ratio of Colony_mean to Control
  expData_raw <- data.frame(
    exp_ID = ratio_data$exp_ID,
    suspension = ratio_data$suspension,
    dose = ratio_data$dose,
    strain = ratio_data$strain,
    column = ratio_data$column,
    Spot_ratio = ratio_data$Colony_mean / ratio_data$Colony_control,
    Conidia_ratio = ratio_data$Conidia / ratio_data$Conidia_control
  )

  n_observations <- nrow(expData_raw)
  if (n_observations < 2 || any(!is.finite(expData_raw$Spot_ratio)) ||
      any(expData_raw$Spot_ratio < 0)) {
    stop("Invalid normalized spot-coverage ratios for ", target_name)
  }

  expData_raw <- expData_raw %>%
    mutate(Spot_ratio_scaled = ifelse(Spot_ratio > 1, 1, Spot_ratio)) %>%
    mutate(Spot_ratio_scaled = (Spot_ratio_scaled * (n_observations - 1) + 0.5) / n_observations) %>%
    mutate(Conidia_ratio_scaled = (Conidia_ratio * (n_observations - 1) + 0.5) / n_observations)

  spot_beta <- betareg(
    Spot_ratio_scaled ~ dose * strain + Conidia_ratio_scaled + exp_ID,
    data = expData_raw,
    link = 'logit'
  )

  spot_beta_reduced <- betareg(
    Spot_ratio_scaled ~ dose + Conidia_ratio_scaled + exp_ID,
    data = expData_raw,
    link = 'logit'
  )

  lrt_result <- likelihood_ratio_test(spot_beta, spot_beta_reduced)
  lrt_result_df <- data.frame(
    target = target_name,
    LR_statistic = lrt_result$LR_statistic,
    df_diff = lrt_result$df_diff,
    p_value = lrt_result$p_value
  )

  UV_dose <- as.vector(unique(expData_raw$dose, nmax = 4))

  dose_comparisons <- emmeans(
    spot_beta,
    ~ dose * strain,
    type = 'response',
    at = list(dose = UV_dose)
  )

  slope_comparisons <- emtrends(
    spot_beta,
    specs = 'strain',
    type = 'response',
    var = 'dose'
  )

  slope_pairwise <- contrast(
    slope_comparisons,
    method = "pairwise",
    type = 'response',
    adjust = 'BH'
  )

  plot_prediction <- as.data.frame(dose_comparisons)

  end_points <- plot_prediction %>%
    group_by(strain) %>%
    filter(dose == max(dose)) %>%
    ungroup()

  g <- ggplot(
    plot_prediction,
    aes(
      x = dose,
      y = emmean * 100,
      color = strain,
      fill = strain,
      shape = strain
    )) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(
      # aes(
      #   ymin = asymp.LCL*100,
      #   ymax = asymp.UCL*100,
      #   width = max(dose)/16
      # ),
      aes(
        ymin = (emmean - SE) * 100,
        ymax = (emmean + SE) * 100,
        width = max(dose) / 16
      ),
      linewidth = 1
    ) +
    geom_text(
      data = end_points,
      aes(
        x = max(dose) * 1.1,
        y = emmean * 100,
        label = rep('xx', nrow(end_points))
      ),
      # hjust = -1,  # 右端に少し離して表示
      show.legend = FALSE,
      size = 4
    ) +
    scale_y_log10(limits = c(NA, 100)) +
    theme_bw(
      base_size = 10
    ) +
    # geom_segment(
    #   aes(
    #     x = 1,
    #     xend = 3,
    #     y = 100,
    #     yend = 100
    #   ),
    #   color = 'red',
    #   linewidth = 0.8
    # ) +
    xlab('UV dose (unit{Jpersquaremeter})') +
    ylab('Estimated spot-coverage ratio (unit{percent})') +
    scale_shape_manual(values = c(15, 0, 16, 1)) +
    scale_color_manual(values = c("#02010C", "#009944", "#0068b7", "#f39800")) +
    theme(
      axis.text = element_text(size = 10, colour = "black"),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 3),
      legend.position = 'bottom',
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(0.0, 'cm'),
      legend.box.spacing = unit(0.0, 'cm'),
      aspect.ratio = 1
      # aspect.ratio = (1+sqrt(5)) / 2
      # aspect.ratio = 4/6
    ) +
    guides(shape = guide_legend(nrow = 2, byrow = TRUE))

  plot(g)

  # =========================
  # Output results
  # =========================

  write.csv(
    raw_csv,
    file = here("04_quantitative_spot_test", "output", target_name, "raw_csv.csv"),
    row.names = FALSE
  )

  write.csv(
    raw_csv_summarize,
    file = here("04_quantitative_spot_test", "output", target_name, "raw_csv_summarize.csv"),
    row.names = FALSE
  )

  write.csv(
    raw_csv_control,
    file = here("04_quantitative_spot_test", "output", target_name, "raw_csv_control.csv"),
    row.names = FALSE
  )

  write.csv(
    expData_raw,
    file = here("04_quantitative_spot_test", "output", target_name, "expData_raw.csv"),
    row.names = FALSE
  )

  write.csv(
    excluded_zero_control,
    file = here("04_quantitative_spot_test", "output", target_name, "excluded_zero_control.csv"),
    row.names = FALSE
  )

  write.csv(
    plot_prediction,
    file = here("04_quantitative_spot_test", "output", target_name, "plot_prediction.csv"),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(summary(dose_comparisons)),
    file = here("04_quantitative_spot_test", "output", target_name, "dose_comparisons.csv"),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(summary(slope_comparisons)),
    file = here("04_quantitative_spot_test", "output", target_name, "slope_comparisons.csv"),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(summary(slope_pairwise)),
    file = here("04_quantitative_spot_test", "output", target_name, "slope_pairwise.csv"),
    row.names = FALSE
  )

  write.csv(
    lrt_result_df,
    file = here("04_quantitative_spot_test", "output", target_name, "likelihood_ratio_test.csv"),
    row.names = FALSE
  )

  sink(here("04_quantitative_spot_test", "output", target_name, "spot_beta_summary.txt"))
  print(summary(spot_beta))
  sink()

  sink(here("04_quantitative_spot_test", "output", target_name, "spot_beta_reduced_summary.txt"))
  print(summary(spot_beta_reduced))
  sink()

  sink(here("04_quantitative_spot_test", "output", target_name, "dose_comparisons_summary.txt"))
  print(summary(dose_comparisons))
  sink()

  sink(here("04_quantitative_spot_test", "output", target_name, "slope_comparisons_summary.txt"))
  print(summary(slope_comparisons))
  sink()

  sink(here("04_quantitative_spot_test", "output", target_name, "slope_pairwise_summary.txt"))
  print(summary(slope_pairwise))
  sink()

  sink(here("04_quantitative_spot_test", "output", target_name, "likelihood_ratio_test.txt"))
  cat("Target:", target_name, "\n")
  cat("Likelihood Ratio Statistic:", lrt_result$LR_statistic, "\n")
  cat("Degrees of Freedom Difference:", lrt_result$df_diff, "\n")
  cat("p-value:", lrt_result$p_value, "\n")
  sink()

  tikz_file <- here(
    "04_quantitative_spot_test",
    "output",
    target_name,
    paste0("qSpot_", target_name, ".tex")
  )

  tikz(
    tikz_file,
    width = 2.75,
    lwdUnit = 72.27 / 96
  )

  tryCatch(
    {
      plot(g)
    },
    finally = {
      dev.off()
    }
  )

  # Return result object
  return(list(
    target = target_name,
    raw_csv = raw_csv,
    raw_csv_summarize = raw_csv_summarize,
    raw_csv_control = raw_csv_control,
    expData_raw = expData_raw,
    spot_beta = spot_beta,
    spot_beta_reduced = spot_beta_reduced,
    lrt_result = lrt_result_df,
    dose_comparisons = dose_comparisons,
    slope_comparisons = slope_comparisons,
    slope_pairwise = slope_pairwise,
    plot_prediction = plot_prediction,
    plot = g
  ))
}

# =========================
# Run all targets
# =========================

qspot_results <- list()

for (i in seq_len(nrow(exp_list))) {

  target_info <- exp_list[i, ]
  target_name <- target_info$target

  qspot_results[[target_name]] <- analyze_qspot_target(target_info)
}

# =========================
# Output combined summaries
# =========================

lrt_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      return(x$lrt_result)
    }
  })
)

write.csv(
  lrt_all,
  file = here("04_quantitative_spot_test", "output", "likelihood_ratio_test_all_targets.csv"),
  row.names = FALSE
)

slope_pairwise_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      as.data.frame(summary(x$slope_pairwise)) %>%
        mutate(target = x$target, .before = 1)
    }
  })
)

write.csv(
  slope_pairwise_all,
  file = here("04_quantitative_spot_test", "output", "slope_pairwise_all_targets.csv"),
  row.names = FALSE
)

slope_comparisons_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      as.data.frame(summary(x$slope_comparisons)) %>%
        mutate(target = x$target, .before = 1)
    }
  })
)

write.csv(
  slope_comparisons_all,
  file = here("04_quantitative_spot_test", "output", "slope_comparisons_all_targets.csv"),
  row.names = FALSE
)

dose_comparisons_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      as.data.frame(summary(x$dose_comparisons)) %>%
        mutate(target = x$target, .before = 1)
    }
  })
)

write.csv(
  dose_comparisons_all,
  file = here("04_quantitative_spot_test", "output", "dose_comparisons_all_targets.csv"),
  row.names = FALSE
)

# Print combined summaries to console
lrt_all
slope_pairwise_all
