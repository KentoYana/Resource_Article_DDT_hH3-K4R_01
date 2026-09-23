# script_killing_test.R
# Kento Yanagisawa
# This script analyzes UV killing-test data using raw area under the curve (AUC).

library("tikzDevice")
library("tidyverse")
library("here")

rm(list = ls())
here::i_am("01_killing_test/script_killing_test.R")

output_dir <- here("01_killing_test", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

calc_auc_trapz <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

raw_csv <- read.csv(here("01_killing_test", "dataset", "killing_data.csv")) %>%
  mutate(Survival_fold = Survival / Fold)

# Technical replicates are averaged before normalization. Each biological
# replicate is one experimental date.
expData_raw <- raw_csv %>%
  group_by(Date, Strain, UV_dose) %>%
  summarise(
    Suv_samples = sum(!is.na(Survival_fold)),
    Suv_ave = mean(Survival_fold, na.rm = TRUE),
    Suv_se = sd(Survival_fold, na.rm = TRUE) / sqrt(Suv_samples),
    .groups = "drop"
  )

Control_0J <- expData_raw %>%
  filter(UV_dose == 0) %>%
  transmute(Date, Strain, Standard_ave = Suv_ave)

expData_raw <- expData_raw %>%
  left_join(Control_0J, by = c("Date", "Strain")) %>%
  mutate(Suv_rate = Suv_ave / Standard_ave)

strain_list <- read.csv(here("01_killing_test", "dataset", "strain_list.csv"))
expData <- expData_raw %>%
  left_join(strain_list, by = "Strain") %>%
  mutate(genotype = factor(genotype, levels = unique(strain_list$genotype)))

required_doses <- c(0, 100, 200, 400)
aucData <- expData %>%
  group_by(Date, genotype) %>%
  filter(all(required_doses %in% UV_dose)) %>%
  arrange(UV_dose, .by_group = TRUE) %>%
  summarise(
    min_dose = min(UV_dose),
    max_dose = max(UV_dose),
    AUC = calc_auc_trapz(UV_dose, Suv_rate),
    .groups = "drop"
  )

if (n_distinct(aucData$Date) < 2 || any(table(aucData$Date) != 3)) {
  stop("Paired AUC analysis requires every retained date to contain all three genotypes.")
}

aucSummary <- aucData %>%
  group_by(genotype) %>%
  summarise(
    n = sum(!is.na(AUC)),
    mean_AUC = mean(AUC, na.rm = TRUE),
    sd_AUC = sd(AUC, na.rm = TRUE),
    se_AUC = sd_AUC / sqrt(n),
    .groups = "drop"
  )

auc_wide <- aucData %>%
  select(Date, genotype, AUC) %>%
  pivot_wider(names_from = genotype, values_from = AUC)
genotypes <- levels(expData$genotype)
pair_indices <- combn(seq_along(genotypes), 2)

aucPairwise <- bind_rows(lapply(seq_len(ncol(pair_indices)), function(i) {
  genotype_1 <- genotypes[pair_indices[1, i]]
  genotype_2 <- genotypes[pair_indices[2, i]]
  difference <- auc_wide[[genotype_1]] - auc_wide[[genotype_2]]
  test_result <- t.test(auc_wide[[genotype_1]], auc_wide[[genotype_2]], paired = TRUE)
  tibble(
    genotype_1 = genotype_1,
    genotype_2 = genotype_2,
    n_pairs = sum(complete.cases(auc_wide[[genotype_1]], auc_wide[[genotype_2]])),
    estimate = mean(difference, na.rm = TRUE),
    SE = sd(difference, na.rm = TRUE) / sqrt(sum(!is.na(difference))),
    df = unname(test_result$parameter),
    lower.CL = test_result$conf.int[1],
    upper.CL = test_result$conf.int[2],
    p.value.raw = test_result$p.value
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

# Means and SEs are calculated across biological-replicate normalized ratios.
# At 0 J every replicate equals one, so its SE is exactly zero.
plotData <- expData %>%
  group_by(genotype, UV_dose) %>%
  summarise(
    n = sum(!is.na(Suv_rate)),
    Survival = mean(Suv_rate, na.rm = TRUE) * 100,
    Survival_se = sd(Suv_rate, na.rm = TRUE) / sqrt(n) * 100,
    .groups = "drop"
  )

auc_labels <- tibble(
  genotype = factor(genotypes, levels = genotypes),
  auc_group = c("a", "b", "a")
) %>%
  left_join(
    plotData %>% group_by(genotype) %>% filter(UV_dose == max(UV_dose)) %>% ungroup(),
    by = "genotype"
  )

g <- ggplot(plotData, aes(x = UV_dose, y = Survival, color = genotype, shape = genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = Survival - Survival_se, ymax = Survival + Survival_se),
    width = max(plotData$UV_dose) / 16,
    linewidth = 1
  ) +
  geom_text(
    data = auc_labels,
    aes(x = max(plotData$UV_dose) * 1.075, y = Survival, label = auc_group),
    show.legend = FALSE,
    size = 4
  ) +
  scale_y_log10() +
  theme_bw(base_size = 10) +
  xlab("UV dose (unit{Jpersquaremeter})") +
  ylab("Survival rate(unit{percent})") +
  scale_shape_manual(values = c(15, 0, 16)) +
  scale_color_manual(values = c("#02010C", "#0068b7", "#f39800")) +
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

write.csv(raw_csv, here("01_killing_test", "output", "raw_csv_corrected.csv"), row.names = FALSE)
write.csv(expData_raw, here("01_killing_test", "output", "expData_raw.csv"), row.names = FALSE)
write.csv(expData, here("01_killing_test", "output", "expData.csv"), row.names = FALSE)
write.csv(plotData, here("01_killing_test", "output", "plotData.csv"), row.names = FALSE)
write.csv(aucData, here("01_killing_test", "output", "aucData.csv"), row.names = FALSE)
write.csv(aucSummary, here("01_killing_test", "output", "auc_summary.csv"), row.names = FALSE)
write.csv(aucPairwise, here("01_killing_test", "output", "auc_pairwise_paired_t_tests.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), here("01_killing_test", "output", "sessionInfo.txt"))

sink(here("01_killing_test", "output", "auc_pairwise_paired_t_tests.txt"))
cat("Raw AUC of survival normalized to 0 J within Date x Strain.\n")
cat("Integration range: 0--400 J/m2; trapezoidal rule.\n")
cat("Paired two-sided t-tests by experimental date; Holm correction across three comparisons.\n\n")
print(aucSummary)
cat("\n")
print(aucPairwise)
sink()

tikz(
  here("01_killing_test", "output", "killing_uv.tex"),
  width = 3.25,
  height = 7,
  lwdUnit = 72.27 / 96
)
plot(g)
dev.off()

aucSummary
aucPairwise
