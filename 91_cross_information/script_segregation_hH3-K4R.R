# script_segregation_hH3-K4R.R
# Kento Yanagisawa
# This script is for the analysis of cross information

# make these packages and their associated functions
# available to use in this script
library("tidyverse")
library("here")

# clear R's brain
rm(list = ls())

# Set project root
here::i_am("91_cross_information/script_segregation_hH3-K4R.R")

# Define directories
output_dir  <- here("91_cross_information", "output")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read dataset
raw_csv <- read_csv(
  here("91_cross_information", "dataset", "cross_data.csv")
)

# Calculate germination and segregation ratios
expData <- raw_csv %>%
  mutate(germination_ratio = germinated_spores / isolated_spores) %>%
  mutate(segregation_ratio = BA_resistant / germinated_spores)

# Summarize results
expData_summary <- expData %>%
  summarise(
    isolated_spore_n = sum(isolated_spores),
    germinated_spores_n = sum(germinated_spores),
    germination_ratio_n = sum(!is.na(germination_ratio)),
    germination_ratio_mean = mean(germination_ratio, na.rm = TRUE),
    germination_ratio_sd = sd(germination_ratio, na.rm = TRUE),
    germination_ratio_se = germination_ratio_sd / sqrt(germination_ratio_n),
    segregation_ratio_n = sum(!is.na(segregation_ratio)),
    segregation_ratio_mean = mean(segregation_ratio, na.rm = TRUE),
    segregation_ratio_sd = sd(segregation_ratio, na.rm = TRUE),
    segregation_ratio_se = segregation_ratio_sd / sqrt(segregation_ratio_n)
  )

# =========================
# Output results
# =========================

write.csv(
  x = raw_csv,
  file = here("91_cross_information", "output", "raw_csv.csv"),
  row.names = FALSE
)

write.csv(
  x = expData,
  file = here("91_cross_information", "output", "expData.csv"),
  row.names = FALSE
)

write.csv(
  x = expData_summary,
  file = here("91_cross_information", "output", "expData_summary.csv"),
  row.names = FALSE
)

# Print summaries to console
expData_summary
