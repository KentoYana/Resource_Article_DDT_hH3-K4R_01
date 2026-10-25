#!/usr/bin/env Rscript
# Adapted from script230525.R by Kento Yanagisawa.
# Reproduce the descriptive race-tube growth curves in Figure S1e.
# Usage: Rscript script_race_tube.R [output-directory]
suppressPackageStartupMessages({
  library(tikzDevice)
  library(tidyverse)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript --vanilla.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript script_race_tube.R [output-directory]")
output_dir <- if (length(args)) args[[1]] else file.path(script_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir)
input_file <- file.path(script_dir, "dataset", "racetube210111.csv")

# read dataset
expData_raw <- read.csv(input_file)

# Check the supplied dataset rather than silently excluding unexpected records.
stopifnot(nrow(expData_raw) == 234L,
          all(c("Time", "Strain", "Length") %in% names(expData_raw)),
          !anyDuplicated(expData_raw[c("Time", "Strain")]),
          sum(is.na(expData_raw$Length)) == 2L)


# Length is already the measured cumulative distance; calculate elapsed time only.
expData_raw <- mutate(expData_raw, DataTime = as.POSIXct(Time, format="%Y-%m-%d %H:%M:%S", tz = "Japan"))
t0 = min(expData_raw$DataTime)
expData_raw <- mutate(expData_raw, growthTime = as.numeric(DataTime-t0, units="hours"))

# Retain measured distances without baseline subtraction or normalization.
expData <- data.frame(
  Strain = expData_raw$Strain,
  Time = expData_raw$growthTime,
  Hyphal_growth = expData_raw$Length
)
expData <- arrange(expData, Time)

# Input genotype
expData <- mutate(expData,
                  genotype = str_replace_all(
                    Strain,
                    c(
                      "C1-T10-28a-1" = "wild-type",
                      "C1-T10-28a-2" = "wild-type",
                      "C1-T10-28a-3" = "wild-type",
                      "KNT-Op3-87a-1" = "hH3-K4R",
                      "KNT-Op3-87a-2" = "hH3-K4R",
                      "KNT-Op3-87a-3" = "hH3-K4R",
                      "FGSC15828-1" = "set-1 KO",
                      "FGSC15828-2" = "set-1 KO",
                      "FGSC15828-3" = "set-1 KO"
                    )
                  )
)

expData_mean <- expData %>%
  group_by(genotype, Time) %>%
  summarise(
    Hyphal_mean = mean(Hyphal_growth, na.rm = TRUE),
    Hyphal_se = sd(Hyphal_growth, na.rm = TRUE)/sqrt(length(Hyphal_growth[!is.na(Hyphal_growth)])),
    n = sum(!is.na(Hyphal_growth)),
    n_missing = sum(is.na(Hyphal_growth)),
    .groups = "drop"
  )

expData_mean$genotype <- factor(expData_mean$genotype,
                                levels = c(
                                  "wild-type",
                                  "hH3-K4R",
                                  "set-1 KO"
                                )
)

g = ggplot(expData_mean,
           aes(x = Time,
               y = Hyphal_mean,
               shape = genotype,
               color = genotype,
               group = genotype
           )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = Hyphal_mean - Hyphal_se,
      ymax = Hyphal_mean + Hyphal_se,
      width = max(Time)/32
    ),
    linewidth = 1
  )+
  theme_bw(
    base_size = 10
  )+
  xlab("Time course (hours)") +
  ylab("Hyphal growth (mm)")+
  scale_shape_manual(values = c(15, 0, 16))+
  scale_color_manual(values = c("#02010C",
                                "#0068b7",
                                "#f39800"
  )
  )+
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 3),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(0.0, 'cm'),
    legend.box.spacing = unit(0.0, 'cm'),
    aspect.ratio = 2/(1+sqrt(5))
  )

expData_velocity <- expData %>%
  filter(Time == max(expData$Time)) %>%
  group_by(genotype) %>%
  summarize(
    velocity = mean(Hyphal_growth / Time),
    velocity_se = sd(Hyphal_growth / Time)/sqrt(3)
  )


# Full precision is retained in memory and in the exported numerical tables.
write.csv(expData, file.path(output_dir, "observations.csv"), row.names = FALSE, na = "NA")
write.csv(expData_mean, file.path(output_dir, "growth_summary.csv"), row.names = FALSE)
write.csv(expData_velocity, file.path(output_dir, "endpoint_velocity.csv"), row.names = FALSE)

# A separate resource plot; this does not overwrite the manuscript presentation copy.
ggsave(file.path(output_dir, "race_tube.pdf"), plot = g, width = 3, height = 3,
       units = "in", device = grDevices::pdf, useDingbats = FALSE)
tikz_file <- file.path(output_dir, "race_tube.tex")
tikz(tikz_file, width = 3, height = 3, lwdUnit = 72.27/96)
print(g)
invisible(dev.off())
# The creation timestamp is not analytical provenance; keep reruns deterministic.
tex <- readLines(tikz_file, warn = FALSE)
tex[1] <- paste0("% Generated by script_race_tube.R with tikzDevice ", packageVersion("tikzDevice"))
writeLines(tex, tikz_file)
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(output_dir, "R-session.txt"))
message("Wrote race-tube outputs to ", output_dir)
