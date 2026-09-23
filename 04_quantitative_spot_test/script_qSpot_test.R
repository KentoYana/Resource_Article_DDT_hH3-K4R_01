# script_qSpot_test.R
# Kento Yanagisawa
# This script is for the analysis of semi-quantitative spot-test(4*6 matrix)

# make these packages and their associated functions
# available to use in this script
library("tikzDevice")
library("tidyverse")
library("here")
library('betareg')

# clear R's brain
rm(list = ls())

# Set project root
here::i_am("04_quantitative_spot_test/script_qSpot_test.R")

# Define directories
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

# Create a compact-letter display directly from the Holm-adjusted pairwise
# tests. Maximal cliques of mutually non-significant genotypes share a letter;
# overlapping cliques naturally produce labels such as "ab".
compact_auc_letters <- function(strains, estimates, pairwise, alpha = 0.05) {
  n <- length(strains)
  nonsignificant <- diag(TRUE, n)
  dimnames(nonsignificant) <- list(strains, strains)
  for (i in seq_len(nrow(pairwise))) {
    first <- match(pairwise$strain_1[i], strains)
    second <- match(pairwise$strain_2[i], strains)
    is_nonsignificant <- pairwise$p.value.Holm[i] >= alpha
    nonsignificant[first, second] <- is_nonsignificant
    nonsignificant[second, first] <- is_nonsignificant
  }

  candidate_cliques <- lapply(seq_len(2^n - 1), function(mask) {
    which(as.logical(intToBits(mask)[seq_len(n)]))
  })
  candidate_cliques <- Filter(function(indices) {
    all(nonsignificant[indices, indices, drop = FALSE])
  }, candidate_cliques)
  maximal_cliques <- Filter(function(indices) {
    !any(vapply(candidate_cliques, function(other) {
      length(other) > length(indices) && all(indices %in% other)
    }, logical(1)))
  }, candidate_cliques)
  maximal_cliques <- maximal_cliques[order(vapply(maximal_cliques, function(indices) {
    -max(estimates[indices])
  }, numeric(1)))]

  letter_codes <- letters[seq_along(maximal_cliques)]
  labels <- vapply(seq_len(n), function(i) {
    paste0(letter_codes[vapply(maximal_cliques, function(indices) i %in% indices, logical(1))], collapse = "")
  }, character(1))
  tibble(strain = strains, auc_group = labels)
}

# Read experiment list
exp_list <- read_csv(here("04_quantitative_spot_test", "dataset", "exp_list.csv"))

# Function to analyze one target
analyze_qspot_target <- function(target_info) {

  target_name <- target_info$target
  target_image_suspension <- target_info$image_suspension

  cat("\n")
  cat("========================================\n")
  cat("Analyzing target:", target_name, "\n")
  cat("Target image-suspension code:", target_image_suspension, "\n")
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

  required_strain_columns <- c("row", "strain", "initial_conidia_per_ml")
  if (!all(required_strain_columns %in% names(strain_list)) ||
      any(!is.finite(strain_list$initial_conidia_per_ml)) ||
      any(strain_list$initial_conidia_per_ml <= 0)) {
    stop("Missing or invalid strain-specific conidial concentration for ", target_name)
  }

  raw_csv <- raw_csv %>%
    mutate(exp_ID = str_extract(csv_path, '/results_.*/UV_')) %>%
    mutate(exp_ID = str_remove_all(exp_ID, '/results_|/UV_')) %>%
    mutate(
      csv_path = str_remove(
        csv_path,
        fixed(paste0(normalizePath(here()), "/"))
      )
    ) %>%
    separate(Image, c('exp_condition', 'image_suspension', 'dose', 'unit', 'photo_ID'), sep = '_') %>%
    mutate(dose = as.numeric(dose)) %>%
    mutate(image_suspension = as.numeric(image_suspension)) %>%
    left_join(strain_list, by = 'row') %>%
    filter(image_suspension == target_image_suspension) %>%
    mutate(Conidia = initial_conidia_per_ml / 5^(column - 1)) %>%
    na.omit()

  raw_csv$strain <- factor(raw_csv$strain, levels = strain_list$strain)

  # calculate mean of each duplicate spot
  raw_csv_summarize <- raw_csv %>%
    group_by(image_suspension, dose, strain, initial_conidia_per_ml, column, exp_ID, Conidia) %>%
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
    image_suspension = raw_csv_control_raw$image_suspension,
    strain = raw_csv_control_raw$strain,
    initial_conidia_per_ml = raw_csv_control_raw$initial_conidia_per_ml,
    column = raw_csv_control_raw$column,
    exp_ID = raw_csv_control_raw$exp_ID,
    Colony_control = raw_csv_control_raw$Colony_mean,
    Conidia_control = max(raw_csv_control_raw$Conidia)
  )

  raw_csv_summarize <- left_join(
    raw_csv_summarize,
    raw_csv_control,
    by = c("image_suspension", "strain", "initial_conidia_per_ml", "column", "exp_ID")
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
    image_suspension = ratio_data$image_suspension,
    dose = ratio_data$dose,
    strain = ratio_data$strain,
    initial_conidia_per_ml = ratio_data$initial_conidia_per_ml,
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
    mutate(
      Conidia_ratio_scaled = (Conidia_ratio * (n_observations - 1) + 0.5) / n_observations,
      strain = factor(strain, levels = strain_list$strain),
      exp_ID = factor(exp_ID)
    )

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

  ## Model diagnostic plots (active graphics device only; not exported via tikzDevice)
  plot(spot_beta, which = 1:4, ask = FALSE)
  plot(spot_beta_reduced, which = 1:4, ask = FALSE)

  lrt_result <- likelihood_ratio_test(spot_beta, spot_beta_reduced)
  lrt_result_df <- data.frame(
    target = target_name,
    LR_statistic = lrt_result$LR_statistic,
    df_diff = lrt_result$df_diff,
    p_value = lrt_result$p_value
  )

  # Generate response-scale curves from the beta regression model. Predictions
  # use the pooled mean conidial covariate. Within each experiment-specific
  # predicted curve, divide by its predicted 0-J value; then average experiments
  # with equal weight. This makes every displayed curve exactly one at 0 J.
  mean_coefficients <- coef(spot_beta, model = "mean")
  mean_covariance <- vcov(spot_beta)[
    names(mean_coefficients), names(mean_coefficients), drop = FALSE
  ]
  observed_doses <- sort(unique(expData_raw$dose))
  dense_doses <- sort(unique(c(
    seq(min(observed_doses), max(observed_doses), length.out = 2001),
    observed_doses
  )))
  exp_levels <- levels(expData_raw$exp_ID)
  conidia_reference <- mean(expData_raw$Conidia_ratio_scaled)

  trapezoid_weights <- function(x) {
    dx <- diff(x)
    c(dx[1] / 2, (head(dx, -1) + tail(dx, -1)) / 2, tail(dx, 1) / 2)
  }
  auc_weights <- trapezoid_weights(dense_doses)

  design_for_strain <- function(strain_value) {
    newdata <- expand.grid(
      dose = dense_doses,
      strain = strain_value,
      exp_ID = exp_levels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        strain = factor(strain, levels = strain_list$strain),
        exp_ID = factor(exp_ID, levels = exp_levels),
        Conidia_ratio_scaled = conidia_reference
      )
    design <- model.matrix(
      ~ dose * strain + Conidia_ratio_scaled + exp_ID,
      data = newdata
    )
    design[, names(mean_coefficients), drop = FALSE]
  }

  normalized_curve <- function(coefficients, design) {
    predictions <- matrix(
      plogis(as.vector(design %*% coefficients)),
      nrow = length(dense_doses)
    )
    normalized_by_experiment <- sweep(
      predictions, 2, predictions[1, ], FUN = "/"
    )
    rowMeans(normalized_by_experiment)
  }

  curve_and_gradient <- function(design) {
    estimate <- normalized_curve(mean_coefficients, design)
    gradient <- vapply(seq_along(mean_coefficients), function(j) {
      step <- 1e-5 * max(1, abs(mean_coefficients[j]))
      upper <- mean_coefficients
      lower <- mean_coefficients
      upper[j] <- upper[j] + step
      lower[j] <- lower[j] - step
      (normalized_curve(upper, design) - normalized_curve(lower, design)) / (2 * step)
    }, numeric(length(dense_doses)))
    list(estimate = estimate, gradient = gradient)
  }

  strain_results <- lapply(strain_list$strain, function(strain_value) {
    result <- curve_and_gradient(design_for_strain(strain_value))
    curve_variance <- rowSums((result$gradient %*% mean_covariance) * result$gradient)
    curve_se <- sqrt(pmax(curve_variance, 0))
    auc_estimate <- sum(auc_weights * result$estimate)
    auc_gradient <- colSums(result$gradient * auc_weights)
    tibble(
      strain = strain_value,
      dose = dense_doses,
      response = result$estimate,
      SE = curve_se,
      lower.CL = result$estimate - qnorm(0.975) * curve_se,
      upper.CL = result$estimate + qnorm(0.975) * curve_se
    ) %>%
      mutate(
        response = if_else(dose == 0, 1, response),
        SE = if_else(dose == 0, 0, SE),
        lower.CL = if_else(dose == 0, 1, lower.CL),
        upper.CL = if_else(dose == 0, 1, upper.CL),
        strain = factor(strain, levels = strain_list$strain)
      ) %>%
      list(
        curve = .,
        auc = tibble(
          target = target_name,
          strain = strain_value,
          min_dose = min(dense_doses),
          max_dose = max(dense_doses),
          reference_Conidia_ratio_scaled = conidia_reference,
          raw_AUC = auc_estimate,
          SE = sqrt(drop(auc_gradient %*% mean_covariance %*% auc_gradient)),
          auc_gradient = list(auc_gradient)
        )
      )
  })

  plot_prediction <- bind_rows(lapply(strain_results, `[[`, "curve"))
  auc_summary <- bind_rows(lapply(strain_results, `[[`, "auc")) %>%
    mutate(
      lower.CL = raw_AUC - qnorm(0.975) * SE,
      upper.CL = raw_AUC + qnorm(0.975) * SE
    )

  auc_gradients <- do.call(rbind, auc_summary$auc_gradient)
  auc_covariance <- auc_gradients %*% mean_covariance %*% t(auc_gradients)
  strain_names <- as.character(auc_summary$strain)
  pair_indices <- combn(seq_along(strain_names), 2)
  auc_pairwise <- bind_rows(lapply(seq_len(ncol(pair_indices)), function(i) {
    first <- pair_indices[1, i]
    second <- pair_indices[2, i]
    estimate <- auc_summary$raw_AUC[first] - auc_summary$raw_AUC[second]
    variance <- auc_covariance[first, first] + auc_covariance[second, second] -
      2 * auc_covariance[first, second]
    SE <- sqrt(max(variance, 0))
    tibble(
      target = target_name,
      strain_1 = strain_names[first],
      strain_2 = strain_names[second],
      estimate = estimate,
      SE = SE,
      z.ratio = estimate / SE,
      lower.CL = estimate - qnorm(0.975) * SE,
      upper.CL = estimate + qnorm(0.975) * SE,
      p.value.raw = 2 * pnorm(-abs(z.ratio))
    )
  })) %>%
    mutate(
      p.value.Holm = p.adjust(p.value.raw, method = "holm"),
      significance.Holm = case_when(
        p.value.Holm < 0.001 ~ "***",
        p.value.Holm < 0.01 ~ "**",
        p.value.Holm < 0.05 ~ "*",
        TRUE ~ "n.s."
      )
    )

  if (length(strain_names) != 4) {
    stop("AUC interaction analysis requires four ordered genotypes for ", target_name)
  }
  interaction_weights <- c(1, -1, -1, 1)
  interaction_estimate <- sum(interaction_weights * auc_summary$raw_AUC)
  interaction_gradient <- colSums(auc_gradients * interaction_weights)
  interaction_se <- sqrt(drop(
    interaction_gradient %*% mean_covariance %*% interaction_gradient
  ))
  auc_interaction <- tibble(
    target = target_name,
    contrast = "(double - background) - (hH3-K4R - wild type)",
    wild_type = strain_names[1],
    hH3_K4R = strain_names[2],
    background = strain_names[3],
    double_mutant = strain_names[4],
    estimate = interaction_estimate,
    SE = interaction_se,
    z.ratio = interaction_estimate / interaction_se,
    lower.CL = interaction_estimate - qnorm(0.975) * interaction_se,
    upper.CL = interaction_estimate + qnorm(0.975) * interaction_se,
    p.value.raw = 2 * pnorm(-abs(z.ratio))
  )

  auc_letters <- compact_auc_letters(
    strain_names,
    auc_summary$raw_AUC,
    auc_pairwise
  )
  auc_summary_output <- auc_summary %>%
    select(-auc_gradient) %>%
    left_join(auc_letters, by = "strain")
  point_prediction <- plot_prediction %>% filter(dose %in% observed_doses)
  strain_plot_labels <- switch(
    target_name,
    "mus-9" = c("wild type", "\\textit{hH3-K4R}", "\\textit{mus-9}", "\\textit{mus-9 hH3-K4R}"),
    "uvs-2" = c("wild type", "\\textit{hH3-K4R}", "\\textit{$\\Delta$uvs-2}", "\\textit{$\\Delta$uvs-2 hH3-K4R}"),
    "mei-3" = c("wild type", "\\textit{hH3-K4R}", "\\textit{$\\Delta$mei-3}", "\\textit{$\\Delta$mei-3 hH3-K4R}"),
    "mus-11" = c("wild type", "\\textit{hH3-K4R}", "\\textit{$\\Delta$mus-11}", "\\textit{$\\Delta$mus-11 hH3-K4R}"),
    "recQ" = c("wild type", "\\textit{hH3-K4R}", "\\textit{qde-3-RIP $\\Delta$recQ2}", "\\textit{qde-3-RIP $\\Delta$recQ2 hH3-K4R}"),
    "mus-26-polh" = c("wild type", "\\textit{hH3-K4R}", "\\textit{$\\Delta$mus-26 $\\Delta$polh}", "\\textit{$\\Delta$mus-26 $\\Delta$polh hH3-K4R}")
  )
  display_doses <- sort(unique(c(
    dense_doses[seq(1, length(dense_doses), length.out = 201) %>% round()],
    observed_doses
  )))
  plot_prediction_display <- plot_prediction %>% filter(dose %in% display_doses)
  endpoint_letters <- plot_prediction %>%
    group_by(strain) %>%
    filter(dose == max(dose)) %>%
    ungroup() %>%
    left_join(auc_letters, by = "strain")

  g <- ggplot(
    plot_prediction_display,
    aes(x = dose, y = response, color = strain, fill = strain, shape = strain)
  ) +
    geom_ribbon(
      aes(ymin = pmax(lower.CL, 0), ymax = upper.CL),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 1) +
    geom_point(data = point_prediction, size = 2.5) +
    geom_text(
      data = endpoint_letters,
      aes(x = max(dense_doses) * 1.1, y = response, label = auc_group),
      show.legend = FALSE,
      size = 4
    ) +
    coord_cartesian(ylim = c(0, 1.06)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    theme_bw(base_size = 10) +
    xlab('UV dose (unit{Jpersquaremeter})') +
    ylab('Predicted response') +
    scale_shape_manual(values = c(15, 0, 16, 1), labels = strain_plot_labels) +
    scale_color_manual(
      values = c("#02010C", "#009944", "#0068b7", "#f39800"),
      labels = strain_plot_labels
    ) +
    scale_fill_manual(
      values = c("#02010C", "#009944", "#0068b7", "#f39800"),
      labels = strain_plot_labels
    ) +
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
    guides(shape = guide_legend(nrow = 2, byrow = FALSE))

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
    auc_summary_output,
    file = here("04_quantitative_spot_test", "output", target_name, "auc_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    auc_pairwise,
    file = here("04_quantitative_spot_test", "output", target_name, "auc_pairwise.csv"),
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

  sink(here("04_quantitative_spot_test", "output", target_name, "auc_analysis_summary.txt"))
  cat("Beta-regression response curves standardized at the pooled mean conidial covariate.\n")
  cat("Each experiment-specific curve was normalized to its predicted 0-J value before equal averaging.\n")
  cat("Raw AUC was integrated over the target-specific dose range by the trapezoidal rule.\n")
  cat("Delta-method uncertainty includes coefficient covariance and 0-J normalization.\n")
  cat("Holm correction was applied across all six pairwise strain contrasts within this target.\n\n")
  print(auc_summary_output)
  cat("\n")
  print(auc_pairwise)
  sink()

  sink(here("04_quantitative_spot_test", "output", target_name, "likelihood_ratio_test.txt"))
  cat("Target:", target_name, "\n")
  cat("Likelihood Ratio Statistic:", lrt_result$LR_statistic, "\n")
  cat("Degrees of Freedom Difference:", lrt_result$df_diff, "\n")
  cat("p-value:", lrt_result$p_value, "\n")
  sink()

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
    auc_summary = auc_summary_output,
    auc_pairwise = auc_pairwise,
    auc_interaction = auc_interaction,
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

auc_pairwise_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      x$auc_pairwise
    }
  })
)

write.csv(
  auc_pairwise_all,
  file = here("04_quantitative_spot_test", "output", "auc_pairwise_all_targets.csv"),
  row.names = FALSE
)

auc_summary_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      x$auc_summary
    }
  })
)

auc_interaction_all <- bind_rows(
  lapply(qspot_results, function(x) {
    if (is.null(x)) NULL else x$auc_interaction
  })
) %>%
  mutate(
    p.value.Holm = p.adjust(p.value.raw, method = "holm"),
    significance.Holm = case_when(
      p.value.Holm < 0.001 ~ "***",
      p.value.Holm < 0.01 ~ "**",
      p.value.Holm < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
  )

write.csv(
  auc_interaction_all,
  file = here("04_quantitative_spot_test", "output", "auc_interaction_all_targets.csv"),
  row.names = FALSE
)

format_interaction_p <- function(p_value) {
  if (p_value < 0.0001) {
    "p < 0.0001"
  } else if (p_value < 0.1) {
    paste0("p = ", formatC(p_value, format = "f", digits = 4))
  } else {
    paste0("p = ", formatC(p_value, format = "f", digits = 3))
  }
}

for (target_name in names(qspot_results)) {
  result <- qspot_results[[target_name]]
  if (is.null(result)) next
  interaction_row <- auc_interaction_all %>% filter(target == target_name)
  write.csv(
    interaction_row,
    file = here("04_quantitative_spot_test", "output", target_name, "auc_interaction.csv"),
    row.names = FALSE
  )
  sink(
    here("04_quantitative_spot_test", "output", target_name, "auc_analysis_summary.txt"),
    append = TRUE
  )
  cat("\nFormal AUC interaction contrast; Holm correction across the six target backgrounds.\n")
  print(interaction_row)
  sink()

  p_label <- format_interaction_p(interaction_row$p.value.Holm)
  max_dose <- max(result$plot_prediction$dose)
  annotated_plot <- result$plot +
    annotate(
      "text", x = max_dose * 0.04, y = 0.1875,
      label = "Genetic interaction:", hjust = 0, size = 3.5
    ) +
    annotate(
      "text", x = max_dose * 0.04, y = 0.0625,
      label = p_label, hjust = 0, size = 3.5
    )
  tikz_file <- here(
    "04_quantitative_spot_test", "output", target_name,
    paste0("qSpot_", target_name, ".tex")
  )
  tikz(tikz_file, width = 2.75, height = 7, lwdUnit = 72.27 / 96)
  tryCatch(plot(annotated_plot), finally = dev.off())
}

write.csv(
  auc_summary_all,
  file = here("04_quantitative_spot_test", "output", "auc_summary_all_targets.csv"),
  row.names = FALSE
)

primary_contrasts <- tibble(
  target = c("mus-9", "uvs-2", "mei-3", "mus-11", "recQ", "mus-26-polh"),
  background = c(
    "mus-9(FK129)", "uvs-2 KO", "mei-3 KO", "mus-11 KO",
    "qde-3-RIP recQ2 KO", "mus-26 polh dKO"
  ),
  plus_hH3_K4R = c(
    "mus-9(FK129) hH3-K4R", "uvs-2 KO hH3-K4R",
    "mei-3 KO hH3-K4R", "mus-11 KO hH3-K4R",
    "qde-3-RIP recQ2 KO hH3-K4R", "mus-26 polh dKO hH3-K4R"
  )
)

auc_primary_contrasts <- auc_pairwise_all %>%
  inner_join(primary_contrasts, by = "target") %>%
  filter(strain_1 == background, strain_2 == plus_hH3_K4R) %>%
  select(
    target, background, plus_hH3_K4R, estimate, SE, z.ratio,
    lower.CL, upper.CL, p.value.raw, p.value.Holm, significance.Holm
  )

write.csv(auc_primary_contrasts,
  file = here("04_quantitative_spot_test", "output", "auc_primary_contrasts.csv"),
  row.names = FALSE)

writeLines(
  capture.output(sessionInfo()),
  here("04_quantitative_spot_test", "output", "sessionInfo.txt")
)

# Print combined summaries to console
lrt_all
auc_primary_contrasts
