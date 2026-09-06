# check_package_versions_and_citations.R

packages <- c(
  # meta-package
  "tidyverse",

  # tidyverse packages used in the scripts
  "dplyr",
  "ggplot2",
  "tidyr",
  "readr",
  "tibble",
  "stringr",

  # other packages used in the scripts
  "tikzDevice",
  "betareg",
  "emmeans",
  "RColorBrewer",
  "MASS",
  "seqinr",
  "patchwork",
  "ggh4x",
  "here",

  # public H3K4 reporter-locus analysis
  "IRanges",
  "rtracklayer",
  "digest",
  "jsonlite"
)

# ---- Package versions ----

pkg_versions <- data.frame(
  Package = packages,
  Installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  Version = vapply(
    packages,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(pkg_versions, row.names = FALSE)

write.csv(pkg_versions, "package_versions.csv", row.names = FALSE)


# ---- BibTeX entries ----

bib_file <- "package_citations.bib"

con <- file(bib_file, open = "w", encoding = "UTF-8")

writeLines(
  c(
    "% BibTeX entries for R and R packages used in the analysis scripts",
    paste0("% Generated on: ", Sys.Date()),
    ""
  ),
  con
)

# R itself
writeLines("% ---- R ----", con)
writeLines(toBibtex(citation()), con)
writeLines("", con)

# Packages
for (pkg in packages) {
  writeLines(paste0("% ---- ", pkg, " ----"), con)

  if (requireNamespace(pkg, quietly = TRUE)) {
    bib <- tryCatch(
      toBibtex(citation(pkg)),
      error = function(e) {
        paste0("% Could not generate BibTeX entry for ", pkg, ": ", e$message)
      }
    )

    writeLines(bib, con)
  } else {
    writeLines(
      paste0("% Package not installed: ", pkg),
      con
    )
  }

  writeLines("", con)
}

close(con)

message("Wrote package versions to: package_versions.csv")
message("Wrote BibTeX entries to: ", bib_file)
