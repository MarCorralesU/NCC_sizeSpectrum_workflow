#!/usr/bin/env Rscript
# =============================================================================
# R requirements for NCC-SPECTRA_sizeSpectrum_workflow
# =============================================================================
# This script checks for (and optionally installs) every third-party R package
# imported by the *.R scripts in this repository.
#
# Usage:
#   Rscript R_requirements.R            # install anything that is missing
#   Rscript R_requirements.R --check    # only report, do not install
#   Rscript R_requirements.R --help
#
# It is safe to re-run: already-installed packages are left untouched.
#
# Base R packages (grid, stats, utils, tools, grDevices, parallel) ship with R
# and are intentionally NOT listed here.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  cat(
    "Usage: Rscript R_requirements.R [--check]\n\n",
    "  --check   Report missing packages without installing them.\n",
    "  --help    Show this message.\n",
    sep = ""
  )
  quit(status = 0)
}

check_only <- "--check" %in% args

# --- Minimum R version -------------------------------------------------------
# The scripts use the native pipe `|>` and the `\\(x)` lambda shorthand,
# both introduced in R 4.1.0.
required_r_version <- "4.1.0"
if (getRversion() < required_r_version) {
  stop(
    "R >= ", required_r_version, " is required (found ", getRversion(), ").\n",
    "The scripts use the native pipe `|>` and `\\(x)` lambda syntax.",
    call. = FALSE
  )
}

# --- Required CRAN packages --------------------------------------------------
# Grouped for readability; `tidyverse` pulls in dplyr, tidyr, readr, purrr,
# stringr, forcats and tibble, but they are listed explicitly for clarity.
cran_packages <- c(
  # Data wrangling / I/O
  "tidyverse", "dplyr", "tidyr", "readr", "purrr", "stringr", "forcats",
  "tibble", "data.table", "reshape2", "plyr", "vroom", "LaF", "chunked",
  "fuzzyjoin", "fs", "here", "magrittr", "rlang", "vctrs", "zoo", "hms",
  "lubridate", "parsedate", "runner", "slider", "suncalc", "NLP",
  # Spatial / interpolation / kriging
  "sf", "terra", "raster", "sp", "gstat", "kriging", "akima", "fields", "MBA",
  "FNN", "geosphere", "metR", "marmap", "ggOceanMaps", "rnaturalearth",
  "rnaturalearthdata", "maps",
  # Statistics / modelling
  "lme4", "lmerTest", "emmeans", "multcomp", "multcompView", "performance",
  "spdep", "vegan", "caret", "ranger", "randomForest", "iml", "broom",
  # Plotting / reporting
  "ggplot2", "cowplot", "patchwork", "colorspace", "colorRamps",
  "RColorBrewer", "viridis", "ggthemes", "ggnewscale", "gt", "knitr",
  "rmarkdown",
  # Utilities / CLI
  "optparse", "scales"
)

# --- CRAN mirror -------------------------------------------------------------
repos <- getOption("repos")
if (is.null(repos) || length(repos) == 0 ||
    identical(unname(repos["CRAN"]), "@CRAN@") || is.na(repos["CRAN"])) {
  repos <- c(CRAN = "https://cloud.r-project.org")
}

# --- Helper ------------------------------------------------------------------
is_installed <- function(pkg) {
  isTRUE(requireNamespace(pkg, quietly = TRUE))
}

missing <- cran_packages[!vapply(cran_packages, is_installed, logical(1))]

cat("R version:", as.character(getRversion()), "\n")
cat("Required third-party packages:", length(cran_packages), "\n")
cat("Already installed:", length(cran_packages) - length(missing), "\n")
cat("Missing:", length(missing), "\n")

if (length(missing) == 0) {
  cat("\nAll R requirements are satisfied.\n")
  quit(status = 0)
}

cat("\nMissing packages:\n")
cat(paste0("  - ", missing, collapse = "\n"), "\n")

if (check_only) {
  cat("\nRun without --check to install them.\n")
  quit(status = 1)
}

# --- Install missing packages ------------------------------------------------
# Use multiple cores for compilation where possible.
options(Ncpus = max(1L, parallel::detectCores() - 1L))

cat("\nInstalling from CRAN (", repos[["CRAN"]], ") ...\n", sep = "")
install.packages(missing, repos = repos, dependencies = TRUE)

# --- Verify ------------------------------------------------------------------
still_missing <- cran_packages[!vapply(cran_packages, is_installed, logical(1))]

if (length(still_missing) > 0) {
  warning(
    "Some packages could not be installed:\n  ",
    paste(still_missing, collapse = ", "),
    "\nCheck the compiler/system libraries (e.g. GDAL, GEOS, PROJ for sf/terra).",
    call. = FALSE
  )
  quit(status = 1)
}

cat("\nAll R requirements are now satisfied.\n")
