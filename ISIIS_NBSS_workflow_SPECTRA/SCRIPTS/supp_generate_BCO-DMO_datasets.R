#!/usr/bin/env Rscript

# supp_generate_BCO-DMO_datasets.R
#
# Purpose:
#   Generate four merged datasets for BCO-DMO submission from SPECTRA workflow outputs.
#
# Outputs:
#   data/BCO-DMO_submission/BCO_DMO_LISST_size_spectra_1m_depthBin.csv
#   data/BCO-DMO_submission/BCO_DMO_ISIIS_binned_abundance_chunk_1.csv
#   data/BCO-DMO_submission/BCO_DMO_ISIIS_NBSS_by_taxa_allParticles_5-10m-bins_thresholded.csv
#   data/BCO-DMO_submission/BCO_DMO_ISIIS_NBSS_by_community_allParticles_5-10m-bins_thresholded.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

base_dir <- "/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA"
output_dir <- file.path(base_dir, "data", "BCO-DMO_submission")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_paths <- list(
  LISST_size_spectra_1m_depthBin = list(
    files = file.path(base_dir, "data", "LISST_data", "LISST_size_spectra", "1m_depthBin"),
    pattern = "\\.csv$",
    recursive = FALSE,
    output = file.path(output_dir, "LISST_size_spectra_1m_depthBin.csv")
  ),
  ISIIS_binned_abundance_chunk_1 = list(
    files = file.path(base_dir, "data", "binned_abundance"),
    pattern = "_binned_abundance_chunk_1\\.csv$",
    recursive = FALSE,
    output = file.path(output_dir, "ISIIS_binned_abundance.csv")
  ),
  ISIIS_NBSS_by_taxa = list(
    files = file.path(base_dir, "data", "size_spectra_allParticles_5-10m-bins_thresholded", "contour"),
    pattern = "ISIIS_NBSS_by_taxa\\.csv$",
    recursive = TRUE,
    output = file.path(output_dir, "ISIIS_NBSS_by_taxa_5m-bins_thresholded.csv")
  ),
  ISIIS_NBSS_seston = list(
    files = file.path(base_dir, "data", "size_spectra_allParticles_5-10m-bins_thresholded", "contour"),
    pattern = "ISIIS_NBSS_by_community\\.csv$",
    recursive = TRUE,
    output = file.path(output_dir, "ISIIS_NBSS_seston_5m-bins_thresholded.csv")
  ),
  ISIIS_NBSS_plankton = list(
    files = file.path(base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded", "contour"),
    pattern = "ISIIS_NBSS_by_community\\.csv$",
    recursive = TRUE,
    output = file.path(output_dir, "ISIIS_NBSS_plankton_5m-bins_thresholded.csv")
  )
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

read_one_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::mutate(source_file = basename(path), .before = 1)
}

read_and_bind_csvs <- function(dir_path, pattern, recursive = FALSE) {
  file_paths <- list.files(
    path = dir_path,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive
  )

  if (length(file_paths) == 0) {
    stop("No files found in: ", dir_path, " with pattern: ", pattern)
  }

  message("Reading ", length(file_paths), " files from ", dir_path)

  purrr::map_dfr(file_paths, read_one_csv)
}

remove_OB_if_present <- function(df) {
  if ("Transect_ID" %in% names(df)) {
    df <- df %>%
      dplyr::filter(.data$Transect_ID != "OB" | is.na(.data$Transect_ID))
  }

  df
}

renumber_group_id_if_present <- function(df) {
  if ("group_id" %in% names(df)) {
    df <- df %>%
      dplyr::rename(group_id_old = .data$group_id) %>%
      dplyr::mutate(
        group_id = cumsum(
          dplyr::coalesce(.data$group_id_old != dplyr::lag(.data$group_id_old), TRUE)
        )
      ) %>%
      dplyr::relocate(.data$group_id, .data$group_id_old)
  }

  df
}

add_UTC_time_if_Pacific_present <- function(df) {
  
  time_cols <- names(df)[stringr::str_detect(names(df), "time_Pacific")]
  
  if (length(time_cols) > 0) {
    
    time_col <- time_cols[1]
    
    df <- df %>%
      dplyr::mutate(
        time_UTC = format(
          lubridate::with_tz(
            lubridate::force_tz(
              lubridate::ymd_hms(.data[[time_col]], tz = "UTC"),
              tzone = "America/Los_Angeles"
            ),
            tzone = "UTC"
          ),
          "%Y-%m-%dT%H:%M:%SZ"
        )
      ) %>%
      dplyr::relocate(.data$time_UTC, .after = dplyr::all_of(time_col))
  }
  
  df
}

add_flag_if_column_exists <- function(df, flag_name, value_col, flag_expression) {
  if (value_col %in% names(df)) {
    df <- df %>%
      dplyr::mutate("{flag_name}" := as.integer({{ flag_expression }}))
  }

  df
}

get_col_by_pattern <- function(df, pattern) {
  matches <- names(df)[stringr::str_detect(names(df), regex(pattern, ignore_case = TRUE))]
  if (length(matches) == 0) return(NULL)
  matches[1]
}

add_QC_flags <- function(df) {
  
  depth_col            <- get_col_by_pattern(df, "Depth|depth_mid")
  pressure_col         <- get_col_by_pattern(df, "Pressure")
  salinity_col         <- get_col_by_pattern(df, "Salinity")
  conductivity_col     <- get_col_by_pattern(df, "Conductivity")
  temperature_col      <- get_col_by_pattern(df, "Temperature")
  horizontal_speed_col <- get_col_by_pattern(df, "Horizontal_Speed")
  vertical_speed_col   <- get_col_by_pattern(df, "Vertical_Speed")
  par_col              <- get_col_by_pattern(df, "(^|_)PAR($|_)")
  ph_col               <- get_col_by_pattern(df, "(^|_)pH($|_)")
  oxygen_col           <- get_col_by_pattern(df, "Oxygen")
  chl_a_col            <- get_col_by_pattern(df, "chl_a")
  
  if (!is.null(depth_col)) {
    df <- df %>%
      dplyr::mutate(
        Depth_flag = as.integer(
          is.na(.data[[depth_col]]) |
            .data[[depth_col]] < 1 |
            .data[[depth_col]] > 110
        )
      )
  }
  
  if (!is.null(pressure_col)) {
    df <- df %>%
      dplyr::mutate(
        Pressure_flag = as.integer(
          is.na(.data[[pressure_col]]) |
            .data[[pressure_col]] < 1 |
            .data[[pressure_col]] > 110
        )
      )
  }
  
  if (!is.null(salinity_col)) {
    df <- df %>%
      dplyr::mutate(
        Salinity_flag = as.integer(
          is.na(.data[[salinity_col]]) |
            .data[[salinity_col]] < 20
        )
      )
  }
  
  if (!is.null(conductivity_col)) {
    df <- df %>%
      dplyr::mutate(
        Conductivity_flag = as.integer(
          is.na(.data[[conductivity_col]]) |
            .data[[conductivity_col]] < 3
        )
      )
  }
  
  if (!is.null(temperature_col)) {
    df <- df %>%
      dplyr::mutate(
        Temperature_flag = as.integer(
          is.na(.data[[temperature_col]])
        )
      )
  }
  
  if (!is.null(horizontal_speed_col)) {
    df <- df %>%
      dplyr::mutate(
        Horizontal_Speed_flag = as.integer(
          is.na(.data[[horizontal_speed_col]]) |
            .data[[horizontal_speed_col]] < -500 |
            .data[[horizontal_speed_col]] > 500
        )
      )
  }
  
  if (!is.null(vertical_speed_col)) {
    df <- df %>%
      dplyr::mutate(
        Vertical_Speed_flag = as.integer(
          is.na(.data[[vertical_speed_col]]) |
            .data[[vertical_speed_col]] < 0 |
            .data[[vertical_speed_col]] > 4000
        )
      )
  }
  
  if (!is.null(par_col)) {
    df <- df %>%
      dplyr::mutate(
        PAR_flag = as.integer(
          is.na(.data[[par_col]]) |
            .data[[par_col]] < 0 |
            .data[[par_col]] > 0.1
        )
      )
  }
  
  if (!is.null(ph_col)) {
    df <- df %>%
      dplyr::mutate(
        pH_flag = as.integer(
          is.na(.data[[ph_col]]) |
            .data[[ph_col]] < 7
        )
      )
  }
  
  if (!is.null(oxygen_col)) {
    df <- df %>%
      dplyr::mutate(
        Oxygen_flag = as.integer(
          is.na(.data[[oxygen_col]])
        )
      )
  }
  
  if (!is.null(chl_a_col)) {
    df <- df %>%
      dplyr::mutate(
        Chl_a_flag = as.integer(
          is.na(.data[[chl_a_col]]) |
            .data[[chl_a_col]] < 0 |
            .data[[chl_a_col]] > 45
        )
      )
  }
  
  flag_cols <- names(df)[stringr::str_detect(names(df), "_flag$")]
  flag_cols <- setdiff(flag_cols, "Flag_totals")
  
  if (length(flag_cols) > 0) {
    df <- df %>%
      dplyr::mutate(
        Flag_totals = rowSums(
          dplyr::across(dplyr::all_of(flag_cols)),
          na.rm = TRUE
        )
      )
  }
  
  df
}

process_dataset <- function(dataset_name, config) {
  message("\nProcessing: ", dataset_name)

  df <- read_and_bind_csvs(
    dir_path = config$files,
    pattern = config$pattern,
    recursive = config$recursive
  ) %>%
    remove_OB_if_present() %>%
    renumber_group_id_if_present() %>%
    add_UTC_time_if_Pacific_present() %>%
    add_QC_flags()
  
  # Drop helper columns from all four BCO-DMO outputs
  df <- df %>%
    dplyr::select(-dplyr::any_of(c("group_id_old", "source_file")))
  
  # Drop LISST abundance column from allParticles NBSS outputs only
  if (dataset_name %in% c("ISIIS_NBSS_by_taxa", "ISIIS_NBSS_seston", "ISIIS_NBSS_plankton")) {
    df <- df %>%
      dplyr::select(-dplyr::any_of("ind_L"))
  }
  
  # Remove duplicate spatial and temporal columns from NBSS products
  if (dataset_name %in% c(
    "ISIIS_NBSS_by_taxa",
    "ISIIS_NBSS_seston",
    "ISIIS_NBSS_plankton"
  )) {
    
    df <- df %>%
      dplyr::select(
        -dplyr::any_of(c("time", "date", "lat", "lon"))
      )
    
  }

  readr::write_csv(df, config$output, na = "")

  message("Rows written: ", nrow(df))
  message("Columns written: ", ncol(df))
  message("Output: ", config$output)

  invisible(df)
}

# -----------------------------------------------------------------------------
# Generate datasets
# -----------------------------------------------------------------------------

purrr::iwalk(input_paths, ~ process_dataset(.y, .x))

message("\nDone. BCO-DMO datasets saved to: ", output_dir)
