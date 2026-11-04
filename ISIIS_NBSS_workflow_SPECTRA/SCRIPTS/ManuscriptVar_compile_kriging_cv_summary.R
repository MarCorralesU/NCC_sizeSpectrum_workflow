#!/usr/bin/env Rscript
# Compile LOOCV summaries from the SPECTRA kriging workflow.
# Usage: Rscript compile_kriging_cv_summary.R [project_root] [interp_tag] [output_csv]
# interp_tag: krige-only (default) or interp_and_krige_dist10

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[[1]] else getwd()
interp_tag <- if (length(args) >= 2) args[[2]] else "krige-only"
if (!interp_tag %in% c("krige-only", "interp_and_krige_dist10")) {
  stop("interp_tag must be 'krige-only' or 'interp_and_krige_dist10'")
}
output_csv <- if (length(args) >= 3) args[[3]] else file.path(base_dir, "figures", interp_tag, "kriging_cv_summary_compiled.csv")

cv_root <- file.path(base_dir, "figures", interp_tag, "size_spectra_coeffs")
if (!dir.exists(cv_root)) stop("CV directory not found: ", cv_root)

variables_keep <- c(
  "slope", "intercept", "log10_total_biovolume_um3_m3",
  "Temperature", "Salinity", "Oxygen_ml_l", "pH", "chl_a_ug_l"
)
columns_keep <- c("Dataset", "Cruise", "Transect_ID", "variable", "ME", "MAE", "RMSE", "R2", "R2_pred")

files <- list.files(cv_root, pattern = "^kriging_cv_summary\\.csv$", recursive = TRUE, full.names = TRUE)
files <- files[grepl("[/\\\\]diagnostic_figsAndTables[/\\\\]kriging_cv_summary\\.csv$", files)]
if (!length(files)) stop("No kriging_cv_summary.csv files found under: ", cv_root)

read_one <- function(path) {
  relative <- substring(normalizePath(path, winslash = "/"), nchar(normalizePath(cv_root, winslash = "/")) + 2)
  parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
  # instrument / cruise / transect / diagnostic_figsAndTables / kriging_cv_summary.csv
  if (length(parts) != 5 || !parts[[4]] == "diagnostic_figsAndTables") {
    warning("Skipping unexpected directory layout: ", path)
    return(NULL)
  }
  instrument <- parts[[1]]
  if (!instrument %in% c("ISIIS", "LISST-ISIIS")) return(NULL)
  dataset <- if (instrument == "ISIIS") "Plankton" else "Seston"

  dat <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("variable", "ME", "MAE", "RMSE", "R2", "R2_pred")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    warning("Skipping file missing columns [", paste(missing, collapse = ", "), "]: ", path)
    return(NULL)
  }
  if ("Cruise" %in% names(dat) && any(!is.na(dat$Cruise) & dat$Cruise != parts[[2]])) {
    warning("Cruise mismatch between CSV and path: ", path)
  }
  if ("Transect_ID" %in% names(dat) && any(!is.na(dat$Transect_ID) & dat$Transect_ID != parts[[3]])) {
    warning("Transect mismatch between CSV and path: ", path)
  }

  dat %>%
    dplyr::filter(variable %in% variables_keep) %>%
    dplyr::mutate(Dataset = dataset, Cruise = parts[[2]], Transect_ID = parts[[3]]) %>%
    dplyr::select(dplyr::all_of(columns_keep))
}

compiled <- purrr::map_dfr(files, read_one) %>%
  dplyr::mutate(
    .dataset_order = match(Dataset, c("Seston", "Plankton")),
    .cruise_order = match(Cruise, c("w22", "s22", "w23", "s23")),
    .variable_order = match(variable, variables_keep)
  ) %>%
  dplyr::arrange(.dataset_order, .cruise_order, Transect_ID, .variable_order) %>%
  dplyr::select(-dplyr::starts_with(".dataset_order"), -dplyr::starts_with(".cruise_order"), -dplyr::starts_with(".variable_order"))

if (!nrow(compiled)) stop("No matching ISIIS/LISST-ISIIS records with requested variables.")
key <- compiled %>% dplyr::count(Dataset, Cruise, Transect_ID, variable) %>% dplyr::filter(n > 1)
if (nrow(key)) warning("Duplicate dataset/cruise/transect/variable combinations detected: ", nrow(key), ". Check source files.")

# Round validation statistics to four significant figures in the exported table.
# Identifiers and variable names remain unchanged.
compiled <- compiled %>%
  dplyr::mutate(
    dplyr::across(
      c(ME, MAE, RMSE, R2, R2_pred),
      ~ ifelse(
        abs(.x) < 0.01 & .x != 0,
        formatC(.x, format = "e", digits = 3),
        format(signif(.x, digits = 4), scientific = FALSE, trim = TRUE)
      )
    )
  )

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(compiled, output_csv)
message("Files found: ", length(files))
message("Rows written: ", nrow(compiled))
message("Output: ", normalizePath(output_csv, winslash = "/"))
