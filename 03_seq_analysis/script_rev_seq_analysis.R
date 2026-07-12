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
p_label <- paste0("p = ", signif(p_lrt, 3))

g <- ggplot(
  indelData,
  aes(
    x = strain_legend,
    y = indel_size,
    color = strain_legend
  )
) +
  geom_violin(
    quantiles = 0.5,
    quantile.linetype = "solid"
  ) +
  geom_point(size = 1, position = position_jitter(width = .15, height = 0)) +
  annotate(
    "text",
    x = 2,
    y = 30,
    label = p_label,
    size = 3.5
  ) +
  theme_bw(base_size = 10) +
  ylab('Absolute indel size (unit{bp})') +
  xlab('') +
  scale_color_manual(values = c("#02010C", "#0068b7")) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 3),
    legend.position = "none",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(0.0, "cm"),
    legend.box.spacing = unit(0.0, "cm"),
    aspect.ratio = 1
  )

plot(g)

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
  dplyr::transmute(
    IR_ID,
    x1 = left_start,
    x2 = left_end,
    y1 = right_start,
    y2 = right_end
  )

facet_df <- tidyr::crossing(
  strain_legend = unique(indelData_pos_count$strain_legend)
)

ir_lines_all <- tidyr::crossing(ir_guides, facet_df)

h <- ggplot(
  indelData_pos_count,
  aes(
    x = mut_start,
    y = mut_finish,
    color = strain_legend
  )
) +
  geom_rect(
    data = ir_lines_all,
    aes(
      xmin = x1,
      xmax = x2,
      ymin = -Inf,
      ymax = Inf,
      fill = IR_ID
    ),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_rect(
    data = ir_lines_all,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = y1,
      ymax = y2,
      fill = IR_ID
    ),
    inherit.aes = FALSE,
    alpha = 0.5
  ) +
  geom_point(
    aes(size = prop_percent),
    shape = 1
  ) +
  scale_size_continuous(
    range = c(1.5, 7.0),
    name = "Revertants (%)"
  ) +
  geom_vline(xintercept = 56, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 56, linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 0.5) +
  theme_bw(base_size = 10) +
  scale_x_continuous(limits = c(0, 90)) +
  scale_y_continuous(limits = c(0, 90)) +
  scale_color_manual(values = c("#02010C", "#0068b7")) +
  scale_fill_brewer(palette = "Set2", name = "IR pair") +
  ylab("3 prime -side boundaries \n relative to the start codon (bp)") +
  xlab("5 prime -side boundaries relative to the start codon (bp)") +
  facet_grid(~ strain_legend) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 3),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(0.0, "cm"),
    legend.box.spacing = unit(0.0, "cm"),
    aspect.ratio = 1
  )

plot(h)

merge_plot <- g + h +
  plot_layout(widths = c(1, 2))

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
  lwdUnit = 72.27 / 96
)
plot(merge_plot)
dev.off()

# Print summaries to console
summary(m_nb)
lrt_nb
p_label