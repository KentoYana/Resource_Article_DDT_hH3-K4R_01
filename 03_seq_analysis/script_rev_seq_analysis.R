# script_rev_analysis.R
# Kento Yanagisawa
# This script is for the analysis of revertant's sequences

# make these packages and their associated functions
# available to use in this script
library("tikzDevice")
library('RColorBrewer')
library("tidyverse")
library("ggfortify")
library("patchwork")
library("here")
library('MASS')
library('seqinr')

# clear R's brain
rm(list = ls())

# Set project root
here::i_am("03_seq_analysis/script_rev_seq_analysis.R")

# Define directories
script_dir  <- here("03_seq_analysis")
dataset_dir <- here("03_seq_analysis", "dataset")
output_dir  <- here("03_seq_analysis", "output")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# palette('Okabe-Ito')
colsBlack  <- brewer.pal(6, 'Greys')
colsOrange <- brewer.pal(6, 'Oranges')
colsBlue   <- brewer.pal(6, 'Blues')
colsPurple <- brewer.pal(5, 'Purples')

# Read dataset
raw_csv <- read.csv(here("03_seq_analysis", "dataset", "reversion_sequence.csv"))

# Input pan-2 sequence
pan2_raw <- read.fasta(
  file = here("03_seq_analysis", "dataset", "pan-2_genomic.fasta"),
  as.string = TRUE,
  forceDNAtolower = FALSE
)

pan2 <- data.frame(
  allele = names(pan2_raw),
  pan2_sequence = as.character(pan2_raw)
)

pan2 <- mutate(
  pan2,
  pan2_sequence = ifelse(
    allele == 'OGW1',
    paste(
      substring(pan2_sequence, 1, 55),
      substring(pan2_sequence, 56, 680),
      substring(pan2_sequence, 681, nchar(pan2_sequence)),
      sep = '-'
    ),
    pan2_sequence
  )
)

# Input strain information
strain_list <- read.csv(here("03_seq_analysis", "dataset", "strain_list.csv"))
strain_list <- left_join(strain_list, pan2, by = 'allele')

# Tidy up raw_csv
expData_raw <- raw_csv[seq(1, 7)]
expData_raw <- left_join(expData_raw, strain_list, by = 'strain')
expData_raw <- na.omit(expData_raw)

# Generate Revertant's information and sequences
expData <- expData_raw %>%
  mutate(rev_sequence = paste(
    substring(pan2_sequence, 1, mut_start),
    mutation,
    substring(pan2_sequence, mut_finish, nchar(pan2_sequence)),
    sep = ''
  )) %>%
  mutate(rev_sequence_extract = paste(
    substring(pan2_sequence, mut_start - 7, mut_start),
    mut_start,
    '(',
    tolower(
      ifelse(
        mutation == '-',
        tolower(substring(pan2_sequence, mut_start + 1, mut_finish - 1)),
        mutation
      )
    ),
    ')',
    substring(pan2_sequence, mut_finish, mut_finish),
    mut_finish,
    substring(pan2_sequence, mut_finish + 1, mut_finish + 7),
    ' (',
    genotype,
    ')',
    sep = ''
  )) %>%
  mutate(pan2_sequence = gsub('-', '', pan2_sequence)) %>%
  mutate(rev_sequence = gsub('-', '', rev_sequence)) %>%
  mutate(mut_size = nchar(rev_sequence) - nchar(pan2_sequence)) %>%
  mutate(indel_size = abs(mut_size)) %>%
  mutate(mut_type = ifelse(
    mut_size == 0,
    'Substitution',
    ifelse(mut_size < 0, 'Deletion', 'Insertion')
  )) %>%
  mutate(delins_chk = case_when(
    mut_type == 'Substitution' ~ FALSE,
    nchar(gsub('-', '', mutation)) != indel_size &
      nchar(gsub('-', '', mutation)) > 0 ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(mut_before8_loc = mut_start - 7) %>%
  mutate(mut_before8 = substring(pan2_sequence, mut_start - 7, mut_start)) %>%
  mutate(mut_after8 = substring(pan2_sequence, mut_finish, mut_finish + 7)) %>%
  mutate(mut_after8_loc = mut_finish + 7)

indelData <- expData %>%
  filter(mut_type != 'Substitution') %>%
  mutate(strain_legend = factor(strain_legend, levels = c("wild-type", "hH3-K4R")))

indel_sample <- indelData %>%
  group_by(mut_size) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  dplyr::select(mut_type, delins_chk, mut_size, rev_sequence_extract)

indelData_summary <- indelData %>%
  group_by(genotype) %>%
  summarise(
    indel_median = median(indel_size),
    indel_max = max(indel_size),
    indel_min = min(indel_size),
    .groups = "drop"
  )

indelData_summary_table <- indelData %>%
  dplyr::count(mut_size, genotype, name = 'n') %>%
  pivot_wider(
    names_from  = genotype,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join(indel_sample, by = 'mut_size') %>%
  arrange(desc(mut_type), -mut_size) %>%
  relocate(mut_type, .before = mut_size)

## Negative binomial GLM
m_nb <- glm.nb(
  indel_size ~ genotype,
  data = indelData
)

## Null model
m_nb_null <- glm.nb(
  indel_size ~ 1,
  data = indelData
)

## Model summary
summary(m_nb)

## Likelihood ratio test
lrt_nb <- anova(
  m_nb_null,
  m_nb,
  test = "Chisq"
)

## Extract p-value
p_lrt <- lrt_nb$`Pr(Chi)`[2]

## Format p-value label
p_label <- paste0("p = ", formatC(p_lrt, digits = 3, format = "f"))

indelData_pos_count <- indelData %>%
  dplyr::count(
    strain_legend,
    mut_start,
    mut_finish,
    name = "n"
  ) %>%
  group_by(strain_legend) %>%
  mutate(
    prop = n / sum(n),
    prop_percent = 100 * prop
  ) %>%
  ungroup()

irs_single <- read.csv(here("03_seq_analysis", "output", "ir_pairs.csv"))

ir_guides <- irs_single %>%
  transmute(
    IR_ID,
    x1 = left_start,
    x2 = left_end,
    y1 = right_start,
    y2 = right_end
  )

ir_arms <- bind_rows(
  irs_single %>%
    transmute(IR_ID, xmin = left_start, xmax = left_end),
  irs_single %>%
    transmute(IR_ID, xmin = right_start, xmax = right_end)
)

ir_arms_by_genotype <- crossing(
  strain_legend = levels(indelData$strain_legend),
  ir_arms
) %>%
  mutate(
    strain_legend = factor(
      strain_legend,
      levels = levels(indelData$strain_legend)
    )
  )

indel_plot_data <- indelData %>%
  arrange(
    strain_legend,
    indel_size,
    mut_start,
    mut_finish,
    mut_type,
    expID,
    seqID,
    colony
  ) %>%
  group_by(strain_legend) %>%
  mutate(event_row = rev(row_number())) %>%
  ungroup()

median_data <- indel_plot_data %>%
  group_by(strain_legend) %>%
  summarise(median_indel_size = median(indel_size), .groups = "drop")

event_cols <- c("Insertion" = "#0068B7", "Deletion" = "#F39800")
median_col <- "#009944"
ir_cols <- c("IR1" = "#DDEBF7", "IR2" = "#FCE4D6", "IR3" = "#E2F0D9")

indel_theme <- theme_bw(base_size = 9) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.2),
    strip.background = element_rect(fill = "#F2F2F2", colour = "black", linewidth = 1.2),
    strip.text.y.left = element_text(angle = 90, face = "plain", size = 8),
    plot.title = element_text(face = "bold", size = 10),
    plot.title.position = "plot",
    plot.margin = margin(5.5, 0, 5.5, 0)
  )

h <- ggplot() +
  geom_rect(
    data = ir_arms_by_genotype,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = IR_ID),
    alpha = 0.58,
    colour = NA
  ) +
  geom_vline(
    xintercept = 56,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "black"
  ) +
  geom_segment(
    data = indel_plot_data,
    aes(
      x = mut_start,
      xend = mut_finish,
      y = event_row,
      yend = event_row,
      colour = mut_type
    ),
    linewidth = 0.9,
    lineend = "round"
  ) +
  geom_point(
    data = indel_plot_data,
    aes(x = mut_start, y = event_row, colour = mut_type),
    shape = 21,
    fill = "white",
    stroke = 0.55,
    size = 1.15
  ) +
  geom_point(
    data = indel_plot_data,
    aes(x = mut_finish, y = event_row, colour = mut_type),
    shape = 21,
    fill = "white",
    stroke = 0.55,
    size = 1.15
  ) +
  facet_grid(
    strain_legend ~ .,
    scales = "free_y",
    switch = "y",
    labeller = as_labeller(c(
      "wild-type" = "wild type",
      "hH3-K4R" = "\\textit{hH3-K4R}"
    ))
  ) +
  scale_x_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, 10),
    expand = expansion(mult = c(0.005, 0.015))
  ) +
  scale_colour_manual(values = event_cols, drop = FALSE) +
  scale_fill_manual(values = ir_cols, drop = FALSE) +
  labs(
    x = "Indel-boundary position relative to the \\textit{pan-2} start codon (bp)",
    y = NULL,
    fill = "Predicted IR pair"
  ) +
  guides(
    colour = "none",
    fill = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  indel_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(5.5, 3, 5.5, 0)
  )

plot(h)

g <- ggplot(indel_plot_data) +
  geom_segment(
    aes(x = 0, xend = indel_size, y = event_row, yend = event_row),
    colour = "#333333",
    linewidth = 0.9,
    lineend = "round"
  ) +
  geom_point(
    aes(x = indel_size, y = event_row),
    colour = "#333333",
    size = 1.15
  ) +
  geom_vline(
    data = median_data,
    aes(xintercept = median_indel_size),
    linewidth = 0.55,
    colour = median_col
  ) +
  geom_label(
    data = median_data,
    aes(
      x = 34,
      y = Inf,
      label = paste0("median = ", median_indel_size, " bp")
    ),
    hjust = 1,
    vjust = 1.8,
    size = 3.5,
    linewidth = 0,
    label.padding = unit(0.08, "lines"),
    fill = "white",
    colour = median_col
  ) +
  facet_grid(strain_legend ~ ., scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 36),
    breaks = c(0, 5, 10, 20, 30, 35)
  ) +
  labs(
    title = paste0("Genotype effect: ", p_label),
    x = "Absolute indel size (bp)",
    y = NULL
  ) +
  indel_theme +
  theme(
    aspect.ratio = (1 + sqrt(5)) / 2,
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "plain", size = 10, hjust = 0.68),
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(5.5, 0, 5.5, 3)
  )

plot(g)

merge_plot <- (h | g) +
  plot_layout(widths = c(3.6, 1), guides = "collect") &
  theme(legend.position = "bottom")

plot(merge_plot)

# Output results

# Output processed data
write.csv(
  raw_csv,
  here("03_seq_analysis", "output", "raw_csv.csv"),
  row.names = FALSE
)

write.csv(
  expData_raw,
  here("03_seq_analysis", "output", "expData_raw.csv"),
  row.names = FALSE
)

write.csv(
  expData,
  here("03_seq_analysis", "output", "expData.csv"),
  row.names = FALSE
)

write.csv(
  indelData,
  here("03_seq_analysis", "output", "indelData.csv"),
  row.names = FALSE
)

write.csv(
  indel_sample,
  here("03_seq_analysis", "output", "indel_sample.csv"),
  row.names = FALSE
)

write.csv(
  indelData_summary,
  here("03_seq_analysis", "output", "indelData_summary.csv"),
  row.names = FALSE
)

write.csv(
  indelData_summary_table,
  here("03_seq_analysis", "output", "indelData_summary_table.csv"),
  row.names = FALSE
)

write.csv(
  indelData_pos_count,
  here("03_seq_analysis", "output", "indelData_pos_count.csv"),
  row.names = FALSE
)

write.csv(
  irs_single,
  here("03_seq_analysis", "output", "ir_pairs_used.csv"),
  row.names = FALSE
)

write.csv(
  ir_guides,
  here("03_seq_analysis", "output", "ir_guides.csv"),
  row.names = FALSE
)

# Output negative binomial GLM results
write.csv(
  as.data.frame(coef(summary(m_nb))),
  here("03_seq_analysis", "output", "negative_binomial_GLM_coefficients.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(lrt_nb),
  here("03_seq_analysis", "output", "negative_binomial_GLM_LRT.csv"),
  row.names = FALSE
)

# Output model summaries as text
sink(here("03_seq_analysis", "output", "negative_binomial_GLM_summary.txt"))
print(summary(m_nb))
sink()

sink(here("03_seq_analysis", "output", "negative_binomial_GLM_LRT.txt"))
print(lrt_nb)
cat("\n")
cat("Formatted p-value label used in plot:\n")
cat(p_label, "\n")
sink()

# Output FASTA files
for (i in strain_list$strain) {
  expData_temp <- expData %>%
    filter(strain == i) %>%
    mutate(rownum = row_number()) %>%
    mutate(dataID = paste('rev', rownum, sep = ''))

  fasta_list <- as.list(expData_temp$rev_sequence)
  names(fasta_list) <- as.vector(expData_temp$dataID)

  fasta_filename <- here(
    "03_seq_analysis",
    "output",
    paste0("expData_", i, ".fasta")
  )

  write.fasta(
    sequences = fasta_list,
    names = names(fasta_list),
    file.out = fasta_filename
  )
}

# Output figures as TikZ
tikz(
  here("03_seq_analysis", "output", "indel-size.tex"),
  width = 3.5,
  lwdUnit = 72.27 / 96
)
plot(g)
dev.off()

tikz(
  here("03_seq_analysis", "output", "indel-dist.tex"),
  height = 3.5,
  lwdUnit = 72.27 / 96
)
plot(h)
dev.off()

tikz(
  here("03_seq_analysis", "output", "indel-plot.tex"),
  width = 8.25,
  height = 6.6,
  lwdUnit = 72.27 / 96
)
plot(merge_plot)
dev.off()

# Print summaries to console
summary(m_nb)
lrt_nb
p_label
