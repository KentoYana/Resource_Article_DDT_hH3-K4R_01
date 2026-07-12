# script_killing_test.R
# Kento Yanagisawa
# This script is for the analysis of killing test

# Make these packages and their associated functions
# available to use in this script
library("tikzDevice")
library('RColorBrewer')
library("tidyverse")
library("ggfortify")
library("patchwork")
library("here")
library("betareg")
library("emmeans")


# Clear R's brain
rm(list = ls())

# Set project root
here::i_am("01_killing_test/script_killing_test.R")

# Define directories
script_dir  <- here("01_killing_test")
dataset_dir <- here("01_killing_test", "dataset")
output_dir  <- here("01_killing_test", "output")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read dataset
raw_csv <- read.csv(here("01_killing_test", "dataset", "killing_data.csv"))

# Correct data
raw_csv <- raw_csv %>%
  mutate(Survival_fold = Survival / Fold)

# Summarize results of technical replicates
expData_raw <- raw_csv %>%
  group_by(Date, Strain, UV_dose) %>%
  summarise(
    Suv_samples = n(),
    Suv_ave = mean(Survival_fold, na.rm = TRUE),
    Suv_se = sd(Survival_fold, na.rm = TRUE) / sqrt(sum(!is.na(Survival_fold))),
    .groups = "drop"
  )

# Standardize using 0 J UV controls and calculate survival rate
Control_0J <- expData_raw %>%
  filter(UV_dose == 0) %>%
  reframe(
    Date,
    Strain,
    Standard_ave = Suv_ave,
    Standard_se = Suv_se
  )

expData_raw <- left_join(
  expData_raw,
  Control_0J,
  by = c("Date", "Strain")
)

# Calculate survival rate
expData_raw <- expData_raw %>%
  mutate(Suv_rate = (Suv_ave * 100) / Standard_ave)

# Calculate propagation of error
propagate_error_division <- function(A, A_err, B, B_err) {
  sqrt((A_err / A)^2 + (B_err / B)^2) * (A / B) * 100
}

expData_raw <- expData_raw %>%
  mutate(
    Suv_rate_se = propagate_error_division(
      Suv_ave,
      Suv_se,
      Standard_ave,
      Standard_se
    )
  )

# Input genotype
strain_list <- read.csv(here("01_killing_test", "dataset", "strain_list.csv"))

expData <- left_join(expData_raw, strain_list, by = "Strain")

expData$genotype <- factor(expData$genotype, levels = strain_list$genotype)

# Summarize results of biological replicates
plotData <- expData %>%
  group_by(genotype, UV_dose) %>%
  summarize(
    Survival = mean(Suv_rate, na.rm = TRUE),
    Survival_se = sqrt(sum(Suv_rate_se^2) / sum(!is.na(Suv_rate_se))),
    .groups = "drop"
  )

# Prepare beta-regression data
expData_beta <- expData %>%
  reframe(
    genotype,
    UV_dose,
    Survival = Suv_rate / 100
  )

expData_beta$Survival <- ifelse(expData_beta$Survival > 1, 1, expData_beta$Survival)

expData_beta$Survival <- (
  expData_beta$Survival * (nrow(expData_beta) - 1) + 0.5
) / nrow(expData_beta)

# Beta regression
killing_beta <- betareg(
  Survival ~ UV_dose * genotype,
  data = expData_beta
)

UV_list <- as.vector(unique(expData_beta$UV_dose, nmax = 4))

emmeans_results <- emmeans(
  killing_beta,
  ~ genotype * UV_dose,
  type = "response",
  at = list(UV_dose = UV_list),
  adjust = "holm"
)

contrast_results <- contrast(
  emmeans_results,
  method = "pairwise",
  by = "UV_dose",
  adjust = "holm"
)

UV_dose_pvalues <- joint_tests(
  killing_beta,
  at = list(UV_dose = UV_list),
  by = "UV_dose"
)

slope_comparisons <- emtrends(
  killing_beta,
  specs = "genotype",
  type = "response",
  var = "UV_dose",
  adjust = "holm"
)

slope_pairwise <- contrast(
  slope_comparisons,
  method = "pairwise",
  type = "response",
  adjust = "holm"
)

end_points <- plotData %>%
  group_by(genotype) %>%
  filter(UV_dose == max(UV_dose)) %>%
  ungroup()

# Plot
g <- ggplot(
  plotData,
  aes(
    x = UV_dose,
    y = Survival,
    color = genotype,
    shape = genotype
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_log10() +
  geom_text(
    data = end_points,
    aes(
      x = max(UV_dose) * 1.075,
      y = Survival,
      label = c("xx", "xx", "xx")
    ),
    show.legend = FALSE,
    size = 4
  ) +
  geom_errorbar(
    aes(
      ymin = Survival - Survival_se,
      ymax = Survival + Survival_se,
      width = max(UV_dose) / 16
    ),
    linewidth = 1
  ) +
  theme_bw(base_size = 10) +
  xlab("UV dose (unit{Jpersquaremeter})") +
  ylab("Survival rate(unit{percent})") +
  scale_shape_manual(values = c(15, 0, 16, 1)) +
  scale_color_manual(
    values = c(
      "#02010C",
      "#0068b7",
      "#f39800",
      "#009944"
    )
  ) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 3),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(0.0, "cm"),
    legend.box.spacing = unit(0.0, "cm"),
    aspect.ratio = 1
  )

plot(g)

# Output results

# Output processed data
write.csv(
  raw_csv,
  here("01_killing_test", "output", "raw_csv_corrected.csv"),
  row.names = FALSE
)

write.csv(
  expData_raw,
  here("01_killing_test", "output", "expData_raw.csv"),
  row.names = FALSE
)

write.csv(
  expData,
  here("01_killing_test", "output", "expData.csv"),
  row.names = FALSE
)

write.csv(
  plotData,
  here("01_killing_test", "output", "plotData.csv"),
  row.names = FALSE
)

write.csv(
  expData_beta,
  here("01_killing_test", "output", "expData_beta.csv"),
  row.names = FALSE
)

# Output statistical results as CSV
write.csv(
  as.data.frame(summary(emmeans_results)),
  here("01_killing_test", "output", "emmeans_results.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(summary(contrast_results)),
  here("01_killing_test", "output", "contrast_results.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(UV_dose_pvalues),
  here("01_killing_test", "output", "UV_dose_pvalues.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(summary(slope_comparisons)),
  here("01_killing_test", "output", "slope_comparisons.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(summary(slope_pairwise)),
  here("01_killing_test", "output", "slope_pairwise.csv"),
  row.names = FALSE
)

# Output model and statistical summaries as text
sink(here("01_killing_test", "output", "killing_beta_summary.txt"))
print(summary(killing_beta))
sink()

sink(here("01_killing_test", "output", "emmeans_summary.txt"))
print(summary(emmeans_results))
sink()

sink(here("01_killing_test", "output", "contrast_summary.txt"))
print(summary(contrast_results))
sink()

sink(here("01_killing_test", "output", "slope_comparisons_summary.txt"))
print(summary(slope_comparisons))
sink()

sink(here("01_killing_test", "output", "slope_pairwise_summary.txt"))
print(summary(slope_pairwise))
sink()

# Output figure as TikZ
tikz(
  here("01_killing_test", "output", "killing_uv.tex"),
  width = 3.25,
  height = 3.25,
  lwdUnit = 72.27 / 96
)
plot(g)
dev.off()

# Print summaries to console
summary(killing_beta)
summary(emmeans_results)
summary(contrast_results)
summary(slope_comparisons)
summary(slope_pairwise)