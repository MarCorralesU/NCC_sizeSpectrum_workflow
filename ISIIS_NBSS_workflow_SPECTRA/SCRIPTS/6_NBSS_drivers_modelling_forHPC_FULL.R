#!/usr/bin/env Rscript

# ============================================================================
# 6_NBSS_drivers_modelling_forHPC_FULL.R
#
# Full HPC workflow for:
#   1. Loading NBSS and binned-abundance data for MEZCAL, MBON_ISIIS-2,
#      MBON_ISIIS-3, and SPECTRA.
#   2. Pairing NBSS samples with abundance samples.
#   3. Collapsing taxonomic abundance columns into functional groups.
#   4. Assigning shelf location.
#   5. Building a pooled analysis data set.
#   6. Hellinger-transforming community composition and running PCA.
#   7. Testing ISIIS-version differences in community composition using the
#      COMPLETE data set (no subsampling).
#   8. Running a full-data, memory-efficient Euclidean PERMANOVA-equivalent
#      permutation test via RDA and a full-data Euclidean centroid PERMDISP
#      permutation test without constructing an n x n distance matrix.
#   9. Saving figures/tables and rendering a self-contained R Markdown report.
#  10. Fitting the random-forest slope model using instrument-standardized
#      functional-group abundances.
#
# IMPORTANT STATISTICAL NOTE:
# ISIIS_version is confounded with project/time/space because the projects do
# not have overlapping spatiotemporal coverage. Therefore the ISIIS-version
# PERMANOVA should be interpreted as a preliminary diagnostic of systematic
# compositional separation, not as a clean causal instrument-calibration test.
# ============================================================================

# ----------------------------- USER SETTINGS ---------------------------------

root_dir <- "/nfs4/FW_HMSC/Cowen_Lab/corralma"
analysis_dir <- file.path(root_dir,  "NBSS_drivers_HPC")
figure_dir <- file.path(analysis_dir, "figures")
table_dir <- file.path(analysis_dir, "tables")
report_dir <- file.path(analysis_dir, "report")
pairing_dir <- file.path(analysis_dir, "NBSS_abundance_pairings")

N_PERM <- 9999
RANDOM_SEED <- 88
PERM_CORES <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(PERM_CORES) || PERM_CORES < 1L) PERM_CORES <- 1L

RUN_PAIRING_DIAGNOSTICS <- TRUE
RUN_RF_SLOPE <- TRUE
RENDER_RMARKDOWN <- TRUE

# Matching tolerances
LAT_TOL <- 0.005
LON_TOL <- 0.005
TIME_TOL_MIN <- 30

# ------------------------------- PACKAGES ------------------------------------

setwd(root_dir)

library(tidyverse)
library(caret)
library(ranger)
library(iml)
library(vegan)
library(scales)

# These are retained because they are used by sourced workflow functions.
library(ggthemes)
library(colorRamps)

# --------------------------- OUTPUT DIRECTORIES -------------------------------

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pairing_dir, recursive = TRUE, showWarnings = FALSE)

message("Analysis directory: ", analysis_dir)
message("Permutation count: ", N_PERM)
message("Permutation cores: ", PERM_CORES)

# --------------------------- SOURCE WORKFLOW CODE ------------------------------

source(file.path(root_dir, "ISIIS_NBSS_workflow_SPECTRA", "SCRIPTS", "functions", "ISIIS_PSS_binning_funcs.R"))
source(file.path(root_dir, "ISIIS_NBSS_workflow_SPECTRA", "SCRIPTS", "functions", "model_functions.R"))

# ------------------------------ PROJECT PATHS ---------------------------------

project_dirs <- c(
  MEZCAL = file.path(root_dir, "ISIIS_NBSS_workflow_MEZCAL"),
  `MBON_ISIIS-2` = file.path(root_dir, "ISIIS_NBSS_workflow_MBON_ISIIS-2"),
  `MBON_ISIIS-3` = file.path(root_dir, "ISIIS_NBSS_workflow_MBON_ISIIS-3"),
  SPECTRA = file.path(root_dir, "ISIIS_NBSS_workflow_SPECTRA")
)

env_vars <- c("Temperature", "Salinity", "Oxygen_ml_l", "pH", "PAR", "chl_a_ug_l")
tow_id_candidates <- c("tow_id", "tow_ID", "Tow_ID", "TOW_ID")
shelf_break_path <- file.path(root_dir, "ISIIS_NBSS_workflow_SPECTRA", "ancillary", "shelf_break_200.csv")

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Standardize the tow-ID column name when one is present.
standardize_tow_id <- function(df) {

  tow_col <- intersect(tow_id_candidates, names(df))

  if (length(tow_col) > 1L) stop("More than one possible tow-ID column was found: ", paste(tow_col, collapse = ", "))

  if (length(tow_col) == 1L && tow_col != "tow_id") df <- df %>% dplyr::rename(tow_id = dplyr::all_of(tow_col))

  df
}

# Correct string-based shelf-location function.
assign_shelf_location <- function(df, shelf_break_path, lon_col_name = "Longitude") {

  shelf_ref <- readr::read_csv(shelf_break_path, show_col_types = FALSE) %>%
    dplyr::mutate(Transect_ID = as.character(Transect_ID), Long_shelf_break_200 = as.numeric(Long_shelf_break_200)) %>%
    dplyr::select(Transect_ID, Long_shelf_break_200)

  df %>%
    dplyr::left_join(shelf_ref, by = "Transect_ID") %>%
    dplyr::mutate(
      shelf_loc_200 = dplyr::case_when(
        is.na(.data[[lon_col_name]]) ~ NA_character_,
        is.na(Long_shelf_break_200) ~ NA_character_,
        .data[[lon_col_name]] < Long_shelf_break_200 ~ "off_shelf",
        .data[[lon_col_name]] >= Long_shelf_break_200 ~ "on_shelf"
      )
    ) %>%
    dplyr::select(-Long_shelf_break_200)
}

# Collapse abundance taxonomic columns into mutually interpretable functional groups.
# Diatoms are explicitly removed from the phytoplankton bucket so they are not double-counted.
collapse_taxa_cols <- function(df, start_name = "appendicularian", end_name = "tunicate_salp", ignore_case = FALSE, detritus_handling = c("keep", "remove")) {

  detritus_handling <- match.arg(detritus_handling)
  nm <- names(df)

  if (ignore_case) {
    nm_l <- tolower(nm)
    s <- match(tolower(start_name), nm_l)
    e <- match(tolower(end_name), nm_l)
  } else {
    s <- match(start_name, nm)
    e <- match(end_name, nm)
  }

  if (is.na(s) || is.na(e) || s > e) stop("Couldn't locate a valid taxa block between '", start_name, "' and '", end_name, "'.")

  data_cols <- nm[s:e]

  patterns <- list(
    detritus = "detritus|art[ie]fact",
    protist = "protist",
    phytoplankton = "phytoplankton",
    diatom = "diatom",
    copepod = "copepod",
    crustacean = "crustacean",
    predatory_gelatinous = "(hydromedusae|chaetognath|siphonophore)",
    grazer_gelatinous = "(appendicularian|tunicate|pyrosome)",
    fish = "fish"
  )

  idx <- lapply(patterns, function(p) grepl(p, data_cols, ignore.case = TRUE))
  idx$phytoplankton <- idx$phytoplankton & !idx$diatom

  detritus_cols <- data_cols[idx$detritus]
  matched_non_det <- Reduce(`|`, idx[c("protist", "phytoplankton", "diatom", "copepod", "crustacean", "predatory_gelatinous", "grazer_gelatinous", "fish")])
  matched_cols <- data_cols[matched_non_det]
  other_cols <- setdiff(data_cols, c(matched_cols, detritus_cols))

  rs <- function(cols) {
    if (!length(cols)) return(rep(0, nrow(df)))
    m <- lapply(df[, cols, drop = FALSE], function(z) suppressWarnings(as.numeric(z)))
    rowSums(as.data.frame(m), na.rm = TRUE)
  }

  grouped <- tibble::tibble(
    protist = rs(data_cols[idx$protist]),
    phytoplankton = rs(data_cols[idx$phytoplankton]),
    diatom = rs(data_cols[idx$diatom]),
    copepod = rs(data_cols[idx$copepod]),
    crustacean = rs(data_cols[idx$crustacean]),
    predatory_gelatinous = rs(data_cols[idx$predatory_gelatinous]),
    grazer_gelatinous = rs(data_cols[idx$grazer_gelatinous]),
    fish = rs(data_cols[idx$fish]),
    other = rs(other_cols)
  )

  keep_left <- if (s > 1) df[, 1:(s - 1), drop = FALSE] else df[, 0, drop = FALSE]
  keep_right <- if (e < ncol(df)) df[, (e + 1):ncol(df), drop = FALSE] else df[, 0, drop = FALSE]
  keep_det <- if (detritus_handling == "keep" && length(detritus_cols)) df[, detritus_cols, drop = FALSE] else NULL

  dplyr::bind_cols(keep_left, grouped, keep_det, keep_right)
}

load_project_NBSS <- function(project_dir, project_name, env_vars) {

  size_spectra_dir <- file.path(project_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded")
  df_list_NBSS <- list.files(size_spectra_dir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)

  message(project_name, ": found ", length(df_list_NBSS), " NBSS files")

  if (length(df_list_NBSS) == 0L) {
    warning("No NBSS files found for ", project_name, " in ", size_spectra_dir)
    return(tibble())
  }

  dat <- purrr::map_dfr(
    df_list_NBSS,
    ~ readr::read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        Cruise, Transect_ID, depth_mid, Latitude, Longitude, group_id,
        total_biovolume, biovol_um3_mid, biovol_um3_width,
        time_Pacific, time_of_day, rel_dist, distance_sample_km,
        slope, intercept,
        dplyr::any_of(c("tow_id", "tow_ID", "Tow_ID", "R2", "residual_logNB", env_vars))
      )
    )
  ) %>%
    dplyr::arrange(time_Pacific) %>%
    dplyr::rename(group_id_old = group_id) %>%
    dplyr::mutate(
      project = project_name,
      group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
      dplyr::across(dplyr::any_of(c("total_biovolume", "biovol_um3_mid", "biovol_um3_width", "slope", "intercept", "R2", "residual_logNB")), as.numeric)
    )

  dat <- standardize_tow_id(dat)

  dat <- dat %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(
      project = dplyr::first(project),
      Cruise = dplyr::first(Cruise),
      Transect_ID = dplyr::first(Transect_ID),
      dplyr::across(dplyr::any_of("tow_id"), dplyr::first),
      time_Pacific = dplyr::first(time_Pacific),
      time_of_day = dplyr::first(time_of_day),
      rel_dist = dplyr::first(rel_dist),
      distance_sample_km = dplyr::first(distance_sample_km),
      depth_mid = dplyr::first(depth_mid),
      Latitude = dplyr::first(Latitude),
      Longitude = dplyr::first(Longitude),
      slope = dplyr::first(slope),
      intercept = dplyr::first(intercept),
      total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
      R2 = dplyr::first(R2),
      residual_logNB = dplyr::first(residual_logNB),
      dplyr::across(dplyr::any_of(env_vars), dplyr::first),
      .groups = "drop"
    ) %>%
    dplyr::mutate(total_biovolume_sum = as.numeric(total_biovolume_sum)) %>%
    dplyr::filter(!is.na(slope), slope < 0) %>%
    dplyr::arrange(Latitude, Longitude, time_Pacific) %>%
    dplyr::distinct(group_id, .keep_all = TRUE)

  dat
}

load_project_abundance <- function(project_dir, project_name) {

  abundance_dir <- file.path(project_dir, "data", "binned_abundance")
  df_list_abundance <- list.files(abundance_dir, pattern = "chunk_1\\.csv$", full.names = TRUE, recursive = TRUE)
  df_list_abundance <- df_list_abundance[!grepl("test", df_list_abundance, ignore.case = TRUE)]

  message(project_name, ": found ", length(df_list_abundance), " abundance files")

  if (length(df_list_abundance) == 0L) {
    warning("No abundance files found for ", project_name, " in ", abundance_dir)
    return(tibble())
  }

  dat_abundance <- purrr::map_dfr(
    df_list_abundance,
    ~ readr::read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        Cruise_name, Transect, cast_id, depth_run_id, depth_mid,
        mean_time_Pacific, mean_Latitude, mean_Longitude,
        dplyr::any_of(c("tow_id", "tow_ID", "Tow_ID")),
        dplyr::matches("^ind_m3")
      )
    )
  )

  dat_abundance <- standardize_tow_id(dat_abundance)

  dat_abundance <- dat_abundance %>%
    dplyr::rename(
      group_id = depth_run_id,
      Cruise = Cruise_name,
      Transect_ID = Transect,
      time_Pacific = mean_time_Pacific,
      Latitude = mean_Latitude,
      Longitude = mean_Longitude,
      Depth = depth_mid
    ) %>%
    dplyr::mutate(
      project = project_name,
      Cruise = dplyr::case_when(
        project_name == "SPECTRA" ~ sub("^SPECTRA_", "", Cruise),
        project_name == "MEZCAL" ~ sub("^MEZCAL_", "", Cruise),
        project_name == "MBON_ISIIS-2" ~ sub("^MBON_ISIIS-2_", "", Cruise),
        project_name == "MBON_ISIIS-3" ~ sub("^MBON_ISIIS-3_", "", Cruise),
        TRUE ~ Cruise
      )
    )

  dat_abundance <- dat_abundance %>%
    depth_bin_calc(depth_bin_width = 5) %>%
    dplyr::relocate(depth_mid, .before = 5)

  abundance_group_cols <- c("Cruise", "Transect_ID", "cast_id", if ("tow_id" %in% names(dat_abundance)) "tow_id", "depth_mid")

  dat_abundance_grouped <- dat_abundance %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(abundance_group_cols))) %>%
    dplyr::summarise(
      project = dplyr::first(project),
      time_Pacific = {
        valid_times <- sort(time_Pacific[!is.na(time_Pacific)])
        if (length(valid_times) == 0L) as.POSIXct(NA) else valid_times[ceiling(length(valid_times) / 2)]
      },
      Latitude = median(Latitude, na.rm = TRUE),
      Longitude = median(Longitude, na.rm = TRUE),
      dplyr::across(dplyr::starts_with("ind_m3_"), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Cruise, Transect_ID, cast_id, depth_mid)

  dat_abundance_grouped
}

flag_time_outliers <- function(dat, gap_multiplier = 20, minimum_gap_minutes = 30, maximum_outlier_block_size = 2) {

  dat <- dat %>%
    dplyr::select(-dplyr::any_of(c("gap_minutes", "gap_threshold_minutes", "time_block", "n_times_in_block", "time_outlier", "typical_interval_minutes", "starts_new_block")))

  required_cols <- c("Cruise", "Transect_ID", "time_Pacific")
  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0L) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  group_cols <- c("Cruise", "Transect_ID", if ("tow_id" %in% names(dat)) "tow_id")

  time_blocks <- dat %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(group_cols)), time_Pacific) %>%
    dplyr::filter(!is.na(time_Pacific)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::arrange(time_Pacific, .by_group = TRUE) %>%
    dplyr::mutate(
      gap_minutes = as.numeric(difftime(time_Pacific, dplyr::lag(time_Pacific), units = "mins")),
      typical_interval_minutes = stats::median(gap_minutes[is.finite(gap_minutes) & gap_minutes > 0], na.rm = TRUE),
      gap_threshold_minutes = pmax(minimum_gap_minutes, gap_multiplier * typical_interval_minutes, na.rm = TRUE),
      starts_new_block = dplyr::if_else(is.na(gap_minutes), TRUE, gap_minutes > gap_threshold_minutes),
      time_block = cumsum(starts_new_block)
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), time_block) %>%
    dplyr::mutate(n_times_in_block = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time_outlier = n_times_in_block <= maximum_outlier_block_size) %>%
    dplyr::select(dplyr::all_of(group_cols), time_Pacific, gap_minutes, gap_threshold_minutes, time_block, n_times_in_block, time_outlier)

  dat %>%
    dplyr::left_join(time_blocks, by = c(group_cols, "time_Pacific"))
}

add_nearest_abundance <- function(dat_NBSS, dat_abundance, lat_tol = 0.005, lon_tol = 0.005, time_tol_min = 30) {

  abundance_cols <- names(dat_abundance)[grepl("^ind_m3", names(dat_abundance))]
  use_tow_id <- "tow_id" %in% names(dat_NBSS) && "tow_id" %in% names(dat_abundance)
  abundance_select <- c("Cruise", "Transect_ID", "depth_mid", if (use_tow_id) "tow_id")

  dat_abundance2 <- dat_abundance %>%
    dplyr::select(
      dplyr::all_of(abundance_select),
      time_Pacific_abundance = time_Pacific,
      Latitude_abundance = Latitude,
      Longitude_abundance = Longitude,
      dplyr::all_of(abundance_cols)
    )

  purrr::map_dfr(seq_len(nrow(dat_NBSS)), function(i) {

    xi <- dat_NBSS[i, ]
    cand <- dat_abundance2 %>%
      dplyr::filter(Cruise == xi$Cruise, Transect_ID == xi$Transect_ID, depth_mid == xi$depth_mid)

    if (use_tow_id) cand <- cand %>% dplyr::filter(tow_id == xi$tow_id)

    cand <- cand %>%
      dplyr::mutate(
        lat_diff = abs(Latitude_abundance - xi$Latitude),
        lon_diff = abs(Longitude_abundance - xi$Longitude),
        time_diff_min = abs(as.numeric(difftime(time_Pacific_abundance, xi$time_Pacific, units = "mins"))),
        match_score = lat_diff / max(lat_tol, 1e-9) + lon_diff / max(lon_tol, 1e-9) + time_diff_min / max(time_tol_min, 1e-9)
      ) %>%
      dplyr::filter(lat_diff <= lat_tol, lon_diff <= lon_tol, time_diff_min <= time_tol_min) %>%
      dplyr::arrange(match_score, time_diff_min, lat_diff, lon_diff)

    if (nrow(cand) == 0L) {
      for (cc in abundance_cols) xi[[cc]] <- NA_real_
      xi$time_Pacific_abundance <- as.POSIXct(NA, tz = "America/Los_Angeles")
      xi$Latitude_abundance <- NA_real_
      xi$Longitude_abundance <- NA_real_
      xi$time_diff_min <- NA_real_
      xi$lat_diff <- NA_real_
      xi$lon_diff <- NA_real_
      return(xi)
    }

    best <- cand[1, ]

    xi %>%
      dplyr::bind_cols(
        best %>%
          dplyr::select(time_Pacific_abundance, Latitude_abundance, Longitude_abundance, dplyr::all_of(abundance_cols), time_diff_min, lat_diff, lon_diff)
      )
  })
}

save_pairing_diagnostics <- function(dat_merged_list, output_root, width = 10, height = 6, dpi = 200) {

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  summary_list <- list()

  for (project_name in names(dat_merged_list)) {

    dat <- dat_merged_list[[project_name]] %>%
      dplyr::mutate(dplyr::across(dplyr::any_of("tow_id"), as.character))

    if (!"time_Pacific_abundance" %in% names(dat)) stop("time_Pacific_abundance is missing from merged data for ", project_name)

    project_dir <- file.path(output_root, project_name)
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

    use_tow_id <- "tow_id" %in% names(dat)
    group_cols <- c("Cruise", "Transect_ID", if (use_tow_id) "tow_id")

    dat <- dat %>%
      dplyr::mutate(match_status = dplyr::if_else(is.na(time_Pacific_abundance), "No abundance match", "Matched"))

    group_keys <- dat %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)))

    project_summary <- dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
      dplyr::summarise(
        n_NBSS = dplyr::n(),
        n_matched = sum(match_status == "Matched"),
        n_unmatched = sum(match_status == "No abundance match"),
        percent_unmatched = 100 * n_unmatched / n_NBSS,
        .groups = "drop"
      )

    summary_list[[project_name]] <- project_summary %>%
      dplyr::mutate(project = project_name, dplyr::across(dplyr::any_of("tow_id"), as.character), .before = 1)

    for (i in seq_len(nrow(group_keys))) {

      key <- group_keys[i, ]
      dat_plot <- dat %>%
        dplyr::semi_join(key, by = group_cols) %>%
        dplyr::filter(!is.na(time_Pacific), !is.na(depth_mid), !is.na(Longitude)) %>%
        dplyr::arrange(time_Pacific, match_status == "No abundance match")

      if (nrow(dat_plot) == 0L) next

      cruise <- as.character(key$Cruise[[1]])
      transect <- as.character(key$Transect_ID[[1]])
      tow <- if (use_tow_id) as.character(key$tow_id[[1]]) else NULL
      n_total <- nrow(dat_plot)
      n_unmatched <- sum(dat_plot$match_status == "No abundance match")
      pct_unmatched <- 100 * n_unmatched / n_total
      title_text <- paste0(project_name, " | ", cruise, " | ", transect, if (use_tow_id) paste0(" | tow ", tow) else "")
      subtitle_text <- paste0("NBSS samples: ", n_total, " | Unmatched: ", n_unmatched, " (", round(pct_unmatched, 1), "%)")

      p_time <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = time_Pacific, y = depth_mid, color = match_status)) +
        ggplot2::geom_point(size = 0.6, alpha = 0.35) +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_color_manual(values = c("Matched" = "grey40", "No abundance match" = "red")) +
        ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d\n%H:%M") +
        ggplot2::labs(title = title_text, subtitle = subtitle_text, x = "Time Pacific", y = "Depth midpoint (m)", color = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")

      p_lon <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = Longitude, y = depth_mid, color = match_status)) +
        ggplot2::geom_point(size = 0.6, alpha = 0.35) +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_color_manual(values = c("Matched" = "grey40", "No abundance match" = "red")) +
        ggplot2::labs(title = title_text, subtitle = subtitle_text, x = "Longitude", y = "Depth midpoint (m)", color = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "top")

      filename_base <- paste0(project_name, "_", cruise, "_", transect, if (use_tow_id) paste0("_tow_", tow) else "")
      filename_base <- gsub("[^A-Za-z0-9_.-]", "_", filename_base)

      ggplot2::ggsave(file.path(project_dir, paste0(filename_base, "_TIME.png")), p_time, width = width, height = height, dpi = dpi)
      ggplot2::ggsave(file.path(project_dir, paste0(filename_base, "_LONGITUDE.png")), p_lon, width = width, height = height, dpi = dpi)
    }

    readr::write_csv(project_summary, file.path(project_dir, paste0(project_name, "_pairing_summary.csv")))
  }

  pairing_summary <- dplyr::bind_rows(summary_list)
  readr::write_csv(pairing_summary, file.path(output_root, "all_projects_pairing_summary.csv"))
  pairing_summary
}

# Marginal 95% confidence intervals for the mean PC1 and PC2 score of each sampling unit.
summarize_pca_centroid_ci <- function(df, group_cols) {

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n = dplyr::n(),
      PC1_mean = mean(PC1, na.rm = TRUE),
      PC2_mean = mean(PC2, na.rm = TRUE),
      PC1_sd = stats::sd(PC1, na.rm = TRUE),
      PC2_sd = stats::sd(PC2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      tcrit = dplyr::if_else(n > 1, stats::qt(0.975, df = n - 1), NA_real_),
      PC1_se = PC1_sd / sqrt(n),
      PC2_se = PC2_sd / sqrt(n),
      PC1_low = PC1_mean - tcrit * PC1_se,
      PC1_high = PC1_mean + tcrit * PC1_se,
      PC2_low = PC2_mean - tcrit * PC2_se,
      PC2_high = PC2_mean + tcrit * PC2_se
    )
}

plot_pca_centroid_ci <- function(summary_df, title_text) {

  ggplot2::ggplot(summary_df, ggplot2::aes(x = PC1_mean, y = PC2_mean, color = ISIIS_version)) +
    ggplot2::geom_segment(ggplot2::aes(x = PC1_low, xend = PC1_high, y = PC2_mean, yend = PC2_mean), linewidth = 0.25, alpha = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = PC1_mean, xend = PC1_mean, y = PC2_low, yend = PC2_high), linewidth = 0.25, alpha = 0.5) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8) +
    ggplot2::facet_wrap(~ project, scales = "free") +
    ggplot2::labs(title = title_text, subtitle = "Points are sampling-unit centroids; bars are marginal 95% CIs for mean PC1 and PC2", x = "Mean PC1", y = "Mean PC2", color = "ISIIS version") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}

# Full-data Euclidean PERMDISP2 centroid test without an n x n distance matrix.
# For Euclidean Hellinger coordinates, the distance from each observation to its
# group centroid can be calculated directly in the original coordinate space.
# The permutation step reproduces vegan::permutest.betadisper's residual-
# permutation F test once those centroid distances are available.
permdisp_euclidean_centroid <- function(X, group, permutations = 9999, seed = 88, cores = 1) {

  X <- as.matrix(X)
  group <- droplevels(as.factor(group))

  if (nlevels(group) < 2L) stop("PERMDISP requires at least two groups.")
  if (anyNA(X)) stop("X contains missing values.")
  if (anyNA(group)) stop("group contains missing values.")

  group_index <- as.integer(group)
  group_sizes <- tabulate(group_index, nbins = nlevels(group))
  centroid_sums <- rowsum(X, group_index, reorder = FALSE)
  centroids <- centroid_sums / group_sizes
  observation_centroids <- centroids[group_index, , drop = FALSE]
  distances <- sqrt(rowSums((X - observation_centroids)^2))

  mod <- stats::lm(distances ~ group)
  mod_aov <- stats::anova(mod)
  mod_qr <- mod$qr
  n <- length(distances)
  p <- mod_qr$rank
  rdf <- n - p
  resids <- qr.resid(mod_qr, distances)
  F0 <- unname(summary(mod)$fstatistic[[1]])
  EPS <- sqrt(.Machine$double.eps)

  one_perm <- function() {
    take <- sample.int(n)
    p_resid <- resids[take]
    fitted_perm <- qr.fitted(mod_qr, p_resid)
    mss <- sum((fitted_perm - mean(fitted_perm))^2)
    residual_perm <- qr.resid(mod_qr, p_resid)
    rss <- sum(residual_perm^2)
    resvar <- rss / rdf
    (mss / (p - 1)) / resvar
  }

  set.seed(seed)

  if (cores > 1L && .Platform$OS.type == "unix") {
    n_cores <- min(cores, permutations)
    chunks <- split(seq_len(permutations), rep(seq_len(n_cores), length.out = permutations))
    Fperm_list <- parallel::mclapply(
      seq_along(chunks),
      function(j) {
        set.seed(seed + j * 100003L)
        vapply(chunks[[j]], function(k) one_perm(), numeric(1))
      },
      mc.cores = n_cores,
      mc.preschedule = TRUE
    )
    Fperm <- unlist(Fperm_list, use.names = FALSE)
  } else {
    Fperm <- vapply(seq_len(permutations), function(k) one_perm(), numeric(1))
  }

  pval <- (sum(Fperm >= F0 - EPS) + 1) / (length(Fperm) + 1)

  group_summary <- tibble::tibble(group = group, distance_to_centroid = distances) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n = dplyr::n(), mean_distance = mean(distance_to_centroid), sd_distance = stats::sd(distance_to_centroid), median_distance = stats::median(distance_to_centroid), .groups = "drop")

  permutation_table <- tibble::tibble(
    term = c("Groups", "Residuals"),
    Df = c(mod_aov$Df[1], mod_aov$Df[2]),
    SumSq = c(mod_aov$`Sum Sq`[1], mod_aov$`Sum Sq`[2]),
    MeanSq = c(mod_aov$`Mean Sq`[1], mod_aov$`Mean Sq`[2]),
    F = c(F0, NA_real_),
    N_Perm = c(permutations, NA_integer_),
    p_value = c(pval, NA_real_)
  )

  list(
    distances = distances,
    group = group,
    group_summary = group_summary,
    permutation_table = permutation_table,
    observed_F = F0,
    p_value = pval,
    n_permutations = permutations
  )
}

fit_rf_ale <- function(dat, response, predictors, seed = 88, folds = 5, tune_length = 5) {

  model_data <- dat %>%
    dplyr::select(dplyr::all_of(c(predictors, response))) %>%
    tidyr::drop_na() %>%
    droplevels()

  message("Response: ", response)
  message("Rows used: ", nrow(model_data))
  message("Rows removed because of NA: ", nrow(dat) - nrow(model_data))

  factor_level_check <- model_data %>%
    dplyr::select(where(is.factor)) %>%
    purrr::map_dfr(
      ~ tibble::tibble(n_levels_defined = nlevels(.x), n_levels_present = dplyr::n_distinct(.x), levels_present = paste(unique(as.character(.x)), collapse = ", ")),
      .id = "variable"
    )

  print(factor_level_check)

  zero_variance_vars <- model_data %>%
    dplyr::select(-dplyr::all_of(response)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ dplyr::n_distinct(.x))) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "n_unique") %>%
    dplyr::filter(n_unique < 2) %>%
    dplyr::pull(variable)

  if (length(zero_variance_vars) > 0L) {
    message("Removing zero-variance predictors: ", paste(zero_variance_vars, collapse = ", "))
    model_data <- model_data %>% dplyr::select(-dplyr::all_of(zero_variance_vars))
  }

  set.seed(seed)
  ctrl <- caret::trainControl(method = "cv", number = folds)
  formula_rf <- stats::as.formula(paste(response, "~ ."))

  rf_model <- caret::train(
    formula_rf,
    data = model_data,
    method = "ranger",
    trControl = ctrl,
    tuneLength = tune_length,
    importance = "permutation",
    metric = "RMSE"
  )

  var_imp <- caret::varImp(rf_model)
  X <- model_data %>% dplyr::select(-dplyr::all_of(response))

  pred <- iml::Predictor$new(
    model = rf_model,
    data = X,
    y = model_data[[response]],
    predict.function = function(m, newdata) predict(m, newdata)
  )

  ale_data <- function(feature) {
    fe <- iml::FeatureEffect$new(pred, feature = feature, method = "ale")
    out <- fe$results
    out$feature <- feature
    out
  }

  list(model_data = model_data, model = rf_model, variable_importance = var_imp, predictor = pred, ale = ale_data, removed_zero_variance_predictors = zero_variance_vars)
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

message("Loading NBSS data...")
dat_NBSS_list <- purrr::imap(project_dirs, ~ load_project_NBSS(project_dir = .x, project_name = .y, env_vars = env_vars))

message("Loading abundance data...")
dat_abundance_list <- purrr::imap(project_dirs, ~ load_project_abundance(project_dir = .x, project_name = .y))

# ==============================================================================
# REMOVE ISOLATED ABUNDANCE-TIME OUTLIERS
# ==============================================================================

message("Flagging abundance time outliers...")
dat_abundance_list <- purrr::map(dat_abundance_list, ~ flag_time_outliers(.x, gap_multiplier = 20, minimum_gap_minutes = 30, maximum_outlier_block_size = 2))

outlier_summary <- purrr::imap_dfr(dat_abundance_list, function(dat, project_name) {
  tibble::tibble(project = project_name, n_rows = nrow(dat), n_time_outliers = sum(dat$time_outlier, na.rm = TRUE))
})
readr::write_csv(outlier_summary, file.path(table_dir, "abundance_time_outlier_summary.csv"))

message("Removing flagged abundance time outliers...")
dat_abundance_list <- purrr::map(
  dat_abundance_list,
  ~ .x %>%
    dplyr::filter(!time_outlier) %>%
    dplyr::select(-gap_minutes, -gap_threshold_minutes, -time_block, -n_times_in_block, -time_outlier)
)

# ==============================================================================
# PAIR NBSS AND ABUNDANCE
# ==============================================================================

message("Pairing NBSS and abundance samples...")
dat_merged_list <- purrr::imap(
  dat_NBSS_list,
  ~ add_nearest_abundance(
    dat_NBSS = .x,
    dat_abundance = dat_abundance_list[[.y]],
    lat_tol = LAT_TOL,
    lon_tol = LON_TOL,
    time_tol_min = TIME_TOL_MIN
  )
)

pairing_loss <- purrr::imap_dfr(dat_merged_list, function(dat, project_name) {
  tibble::tibble(
    project = project_name,
    n_NBSS_original = nrow(dat),
    n_without_abundance = sum(is.na(dat$time_Pacific_abundance)),
    n_with_abundance = sum(!is.na(dat$time_Pacific_abundance)),
    percent_removed = 100 * n_without_abundance / n_NBSS_original
  )
})
readr::write_csv(pairing_loss, file.path(table_dir, "NBSS_abundance_pairing_loss.csv"))

if (RUN_PAIRING_DIAGNOSTICS) {
  message("Saving all-data pairing diagnostics...")
  save_pairing_diagnostics(dat_merged_list, file.path(pairing_dir, "all_data"))
}

dat_merged_list_unfiltered <- dat_merged_list

dat_merged_list <- purrr::map(dat_merged_list, ~ .x %>% dplyr::filter(!is.na(time_Pacific_abundance)))

if (RUN_PAIRING_DIAGNOSTICS) {
  message("Saving paired-only diagnostics...")
  save_pairing_diagnostics(dat_merged_list, file.path(pairing_dir, "paired_only"))
}

# ==============================================================================
# DISCOVER TAXONOMIC COLUMN LIMITS AND COLLAPSE TO FUNCTIONAL GROUPS
# ==============================================================================

message("Discovering project-specific abundance-column limits...")
column_info_by_project <- purrr::imap_dfr(dat_merged_list, function(dat, project_name) {
  tibble::tibble(column_name = names(dat), project = project_name, column_order = seq_along(names(dat)))
}) %>%
  tidyr::pivot_wider(names_from = project, values_from = column_order, names_glue = "{project}_order") %>%
  dplyr::mutate(
    MEZCAL = dplyr::if_else(!is.na(`MEZCAL_order`), column_name, NA_character_),
    `MBON_ISIIS-2` = dplyr::if_else(!is.na(`MBON_ISIIS-2_order`), column_name, NA_character_),
    `MBON_ISIIS-3` = dplyr::if_else(!is.na(`MBON_ISIIS-3_order`), column_name, NA_character_),
    SPECTRA = dplyr::if_else(!is.na(`SPECTRA_order`), column_name, NA_character_)
  ) %>%
  dplyr::select(column_name, MEZCAL, `MEZCAL_order`, `MBON_ISIIS-2`, `MBON_ISIIS-2_order`, `MBON_ISIIS-3`, `MBON_ISIIS-3_order`, SPECTRA, `SPECTRA_order`)

readr::write_csv(column_info_by_project, file.path(table_dir, "column_info_by_project.csv"))

project_names <- names(dat_merged_list)

taxa_anchor_list <- purrr::map(project_names, function(project_name) {

  name_col <- project_name
  order_col <- paste0(project_name, "_order")

  taxa_cols <- column_info_by_project %>%
    dplyr::filter(!is.na(.data[[name_col]])) %>%
    dplyr::filter(grepl("^ind_m3", .data[[name_col]])) %>%
    dplyr::filter(!is.na(.data[[order_col]]))

  if (nrow(taxa_cols) == 0L) stop("No ind_m3 columns found for project ", project_name)

  start_row <- taxa_cols %>% dplyr::slice_min(order_by = .data[[order_col]], n = 1, with_ties = FALSE)
  end_row <- taxa_cols %>% dplyr::slice_max(order_by = .data[[order_col]], n = 1, with_ties = FALSE)

  c(start = start_row[[name_col]][[1]], end = end_row[[name_col]][[1]])
}) %>%
  rlang::set_names(project_names)

taxa_anchor_table <- purrr::imap_dfr(taxa_anchor_list, function(anchors, project_name) {
  tibble::tibble(project = project_name, start = anchors[["start"]], end = anchors[["end"]])
})
readr::write_csv(taxa_anchor_table, file.path(table_dir, "taxa_anchor_table.csv"))

message("Collapsing taxonomic abundance columns...")
dat_merged_list_collapsed <- purrr::imap(dat_merged_list, function(dat, project_name) {
  anchors <- taxa_anchor_list[[project_name]]
  collapse_taxa_cols(df = dat, start_name = anchors[["start"]], end_name = anchors[["end"]], detritus_handling = "remove")
})

collapse_check <- purrr::imap_dfr(dat_merged_list, function(dat, project_name) {
  tibble::tibble(
    project = project_name,
    rows_before = nrow(dat),
    rows_after = nrow(dat_merged_list_collapsed[[project_name]]),
    cols_before = ncol(dat),
    cols_after = ncol(dat_merged_list_collapsed[[project_name]]),
    group_id_unique = nrow(dat_merged_list_collapsed[[project_name]]) == dplyr::n_distinct(dat_merged_list_collapsed[[project_name]]$group_id)
  )
})
readr::write_csv(collapse_check, file.path(table_dir, "collapse_check.csv"))

# ==============================================================================
# SHELF LOCATION
# ==============================================================================

message("Assigning shelf location...")
dat_merged_list_collapsed <- purrr::map(
  dat_merged_list_collapsed,
  ~ assign_shelf_location(.x, shelf_break_path = shelf_break_path, lon_col_name = "Longitude")
)

shelf_location_check <- purrr::imap_dfr(dat_merged_list_collapsed, function(dat, project_name) {
  dat %>%
    dplyr::count(shelf_loc_200) %>%
    dplyr::mutate(project = project_name, .before = 1)
})
readr::write_csv(shelf_location_check, file.path(table_dir, "shelf_location_check.csv"))

# ==============================================================================
# POOLED MODEL DATA
# ==============================================================================

message("Building pooled model data...")
dat_merged_list_collapsed <- purrr::map(dat_merged_list_collapsed, ~ .x %>% dplyr::mutate(dplyr::across(dplyr::any_of("tow_id"), as.character)))

dat_model_all <- dplyr::bind_rows(dat_merged_list_collapsed) %>%
  dplyr::mutate(
    project = as.factor(project),
    Cruise = as.factor(Cruise),
    Transect_ID = as.factor(Transect_ID),
    shelf_loc_200 = as.factor(shelf_loc_200),
    dplyr::across(dplyr::any_of("tow_id"), as.factor),
    dplyr::across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA_real_)),
    ISIIS_version = dplyr::case_when(
      project %in% c("MEZCAL", "MBON_ISIIS-2") ~ "ISIIS-2",
      project %in% c("SPECTRA", "MBON_ISIIS-3") ~ "ISIIS-3",
      TRUE ~ NA_character_
    ),
    ISIIS_version = as.factor(ISIIS_version)
  )

# Functional-group abundance variables.
plankton_vars <- c("protist", "diatom", "copepod", "crustacean", "predatory_gelatinous", "grazer_gelatinous", "fish", "other")

# Remove phytoplankton because the current grouping showed that this column is
# structurally zero after diatoms are separated.
dat_model_all <- dat_model_all %>%
  dplyr::select(-dplyr::any_of(c("phytoplankton", "phytoplankton_scaled")))

# Instrument-specific z-scoring retained for the RF analysis only.
dat_model_all <- dat_model_all %>%
  dplyr::group_by(ISIIS_version) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(plankton_vars), ~ as.numeric(scale(.x)), .names = "{.col}_scaled")) %>%
  dplyr::ungroup()

sample_counts <- dat_model_all %>%
  dplyr::count(ISIIS_version, project, Cruise, name = "n_samples") %>%
  dplyr::arrange(ISIIS_version, project, Cruise)
readr::write_csv(sample_counts, file.path(table_dir, "sample_counts.csv"))

saveRDS(dat_model_all, file.path(analysis_dir, "dat_model_all.rds"), compress = TRUE)

# ==============================================================================
# HELLINGER COMMUNITY PCA -- RAW FUNCTIONAL-GROUP ABUNDANCES
# ==============================================================================

message("Running Hellinger transformation and PCA...")
community_vars <- plankton_vars

missing_community_vars <- setdiff(community_vars, names(dat_model_all))
if (length(missing_community_vars) > 0L) stop("Missing community variables: ", paste(missing_community_vars, collapse = ", "))

community_complete <- dat_model_all %>%
  dplyr::select(group_id, project, ISIIS_version, Cruise, Transect_ID, tow_id, time_Pacific, depth_mid, Latitude, Longitude, slope, intercept, dplyr::all_of(community_vars)) %>%
  tidyr::drop_na(dplyr::all_of(community_vars)) %>%
  dplyr::filter(rowSums(dplyr::across(dplyr::all_of(community_vars))) > 0)

community_matrix <- community_complete %>%
  dplyr::select(dplyr::all_of(community_vars))

if (any(as.matrix(community_matrix) < 0, na.rm = TRUE)) stop("Community matrix contains negative abundance values; Hellinger transformation requires non-negative abundances.")

community_hellinger <- vegan::decostand(community_matrix, method = "hellinger")
community_hellinger_matrix <- as.matrix(community_hellinger)
community_pca <- stats::prcomp(community_hellinger_matrix, center = TRUE, scale. = FALSE)

pca_variance_table <- tibble::tibble(
  PC = paste0("PC", seq_along(community_pca$sdev)),
  variance_explained = 100 * community_pca$sdev^2 / sum(community_pca$sdev^2)
) %>%
  dplyr::mutate(cumulative_variance = cumsum(variance_explained))

pca_variance <- pca_variance_table$variance_explained
pca_loadings <- as.data.frame(community_pca$rotation) %>%
  tibble::rownames_to_column("functional_group")

pca_scores <- as.data.frame(community_pca$x) %>%
  dplyr::bind_cols(community_complete %>% dplyr::select(group_id, project, ISIIS_version, Cruise, Transect_ID, tow_id, time_Pacific, depth_mid, Latitude, Longitude, slope, intercept))

readr::write_csv(pca_variance_table, file.path(table_dir, "PCA_variance_explained.csv"))
readr::write_csv(pca_loadings, file.path(table_dir, "PCA_loadings.csv"))

# ----------------------------- PCA FIGURES ------------------------------------

message("Saving PCA figures...")

p_pca_instrument <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2, color = ISIIS_version)) +
  ggplot2::geom_point(alpha = 0.06, size = 0.25) +
  ggplot2::stat_ellipse(level = 0.95, linewidth = 0.8) +
  ggplot2::labs(
    x = paste0("PC1 (", round(pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance[2], 1), "%)"),
    color = "ISIIS version",
    title = "Hellinger-transformed functional community composition by ISIIS version",
    subtitle = "Small transparent points show all samples; ellipses contain approximately 95% of each fitted bivariate distribution"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top")

p_pca_project <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2, color = project)) +
  ggplot2::geom_point(alpha = 0.05, size = 0.25) +
  ggplot2::stat_ellipse(level = 0.95, linewidth = 0.7) +
  ggplot2::labs(
    x = paste0("PC1 (", round(pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance[2], 1), "%)"),
    color = "Project",
    title = "Hellinger-transformed functional community composition by project",
    subtitle = "Project is partially confounded with instrument version, sampling year, and spatial coverage"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top")

p_scree <- pca_variance_table %>%
  dplyr::mutate(PC = factor(PC, levels = PC)) %>%
  ggplot2::ggplot(ggplot2::aes(x = PC, y = variance_explained)) +
  ggplot2::geom_col() +
  ggplot2::geom_line(ggplot2::aes(y = cumulative_variance, group = 1), linewidth = 0.6) +
  ggplot2::geom_point(ggplot2::aes(y = cumulative_variance), size = 1.5) +
  ggplot2::labs(title = "PCA variance explained", y = "Variance explained / cumulative variance (%)", x = NULL) +
  ggplot2::theme_bw()

p_loadings <- pca_loadings %>%
  dplyr::select(functional_group, PC1, PC2) %>%
  tidyr::pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "loading") %>%
  ggplot2::ggplot(ggplot2::aes(x = reorder(functional_group, loading), y = loading)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~ PC, scales = "free_y") +
  ggplot2::labs(title = "Functional-group PCA loadings", x = "Functional group", y = "Loading") +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(figure_dir, "PCA_by_ISIIS_version.png"), p_pca_instrument, width = 9, height = 7, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "PCA_by_project.png"), p_pca_project, width = 9, height = 7, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "PCA_variance_scree.png"), p_scree, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "PCA_loadings_PC1_PC2.png"), p_loadings, width = 9, height = 7, dpi = 300)

# -------------------- PCA CENTROIDS + 95% MARGINAL CIs -----------------------

message("Calculating PCA sampling-unit centroid confidence intervals...")

pca_ci_cruise <- summarize_pca_centroid_ci(pca_scores, c("ISIIS_version", "project", "Cruise"))
pca_ci_transect <- summarize_pca_centroid_ci(pca_scores, c("ISIIS_version", "project", "Cruise", "Transect_ID"))
pca_ci_tow <- pca_scores %>%
  dplyr::filter(!is.na(tow_id)) %>%
  summarize_pca_centroid_ci(c("ISIIS_version", "project", "Cruise", "Transect_ID", "tow_id"))

readr::write_csv(pca_ci_cruise, file.path(table_dir, "PCA_centroids_Cruise_95CI.csv"))
readr::write_csv(pca_ci_transect, file.path(table_dir, "PCA_centroids_Transect_95CI.csv"))
readr::write_csv(pca_ci_tow, file.path(table_dir, "PCA_centroids_Tow_95CI.csv"))

p_ci_cruise <- plot_pca_centroid_ci(pca_ci_cruise, "PCA community centroids by Cruise")
p_ci_transect <- plot_pca_centroid_ci(pca_ci_transect, "PCA community centroids by Cruise × Transect")
p_ci_tow <- plot_pca_centroid_ci(pca_ci_tow, "PCA community centroids by Cruise × Transect × Tow")

ggplot2::ggsave(file.path(figure_dir, "PCA_centroids_Cruise_95CI.png"), p_ci_cruise, width = 10, height = 7, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "PCA_centroids_Transect_95CI.png"), p_ci_transect, width = 10, height = 7, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "PCA_centroids_Tow_95CI.png"), p_ci_tow, width = 10, height = 7, dpi = 300)

# ==============================================================================
# FULL-DATA PERMANOVA: MEMORY-EFFICIENT RDA EQUIVALENT
# ==============================================================================

message("Running full-data PERMANOVA-equivalent permutation test via RDA...")
message("This uses every community sample; no subsampling is performed.")

# vegan documentation states that with Euclidean distances adonis2 tests are
# identical to anova.cca tests of rda. Because the response is already a
# Hellinger-transformed Euclidean matrix, this avoids materializing the huge
# pairwise distance matrix while retaining the same full-data test.

set.seed(RANDOM_SEED)
rda_ISIIS <- vegan::rda(community_hellinger_matrix ~ ISIIS_version, data = community_complete)
permanova_ISIIS_test <- stats::anova(rda_ISIIS, permutations = N_PERM, parallel = PERM_CORES)
permanova_ISIIS_r2 <- vegan::RsquareAdj(rda_ISIIS)

permanova_ISIIS_table <- as.data.frame(permanova_ISIIS_test) %>%
  tibble::rownames_to_column("term")

permanova_ISIIS_effect <- tibble::tibble(
  predictor = "ISIIS_version",
  R2 = unname(permanova_ISIIS_r2$r.squared),
  adjusted_R2 = unname(permanova_ISIIS_r2$adj.r.squared),
  n_samples = nrow(community_complete),
  permutations = N_PERM
)

readr::write_csv(permanova_ISIIS_table, file.path(table_dir, "PERMANOVA_ISIIS_version_full_data.csv"))
readr::write_csv(permanova_ISIIS_effect, file.path(table_dir, "PERMANOVA_ISIIS_version_R2_full_data.csv"))
saveRDS(rda_ISIIS, file.path(analysis_dir, "rda_ISIIS_full_data.rds"), compress = TRUE)

# Additional project-level diagnostic. This does NOT resolve confounding, but it
# quantifies how much total community variation is associated with the four
# project labels themselves.
message("Running full-data project-level PERMANOVA diagnostic via RDA...")
set.seed(RANDOM_SEED)
rda_project <- vegan::rda(community_hellinger_matrix ~ project, data = community_complete)
permanova_project_test <- stats::anova(rda_project, permutations = N_PERM, parallel = PERM_CORES)
permanova_project_r2 <- vegan::RsquareAdj(rda_project)

permanova_project_table <- as.data.frame(permanova_project_test) %>%
  tibble::rownames_to_column("term")

permanova_project_effect <- tibble::tibble(
  predictor = "project",
  R2 = unname(permanova_project_r2$r.squared),
  adjusted_R2 = unname(permanova_project_r2$adj.r.squared),
  n_samples = nrow(community_complete),
  permutations = N_PERM
)

readr::write_csv(permanova_project_table, file.path(table_dir, "PERMANOVA_project_full_data.csv"))
readr::write_csv(permanova_project_effect, file.path(table_dir, "PERMANOVA_project_R2_full_data.csv"))

# ==============================================================================
# FULL-DATA PERMDISP: EUCLIDEAN CENTROID VERSION WITHOUT PAIRWISE DIST MATRIX
# ==============================================================================

message("Running full-data Euclidean centroid PERMDISP permutation test...")
message("This uses every community sample and avoids the O(n^2) pairwise-distance object.")

set.seed(RANDOM_SEED)
permdisp_ISIIS <- permdisp_euclidean_centroid(
  X = community_hellinger_matrix,
  group = community_complete$ISIIS_version,
  permutations = N_PERM,
  seed = RANDOM_SEED,
  cores = PERM_CORES
)

readr::write_csv(permdisp_ISIIS$permutation_table, file.path(table_dir, "PERMDISP_ISIIS_version_full_data.csv"))
readr::write_csv(permdisp_ISIIS$group_summary, file.path(table_dir, "PERMDISP_ISIIS_version_group_dispersion.csv"))

permdisp_plot_data <- tibble::tibble(
  ISIIS_version = permdisp_ISIIS$group,
  distance_to_centroid = permdisp_ISIIS$distances
)

p_permdisp <- ggplot2::ggplot(permdisp_plot_data, ggplot2::aes(x = ISIIS_version, y = distance_to_centroid)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.25) +
  ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.45) +
  ggplot2::labs(
    title = "Multivariate dispersion by ISIIS version",
    subtitle = paste0("Full-data Euclidean centroid PERMDISP; ", format(nrow(community_complete), big.mark = ","), " samples"),
    x = "ISIIS version",
    y = "Euclidean distance to group centroid in Hellinger space"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(figure_dir, "PERMDISP_ISIIS_version.png"), p_permdisp, width = 7, height = 6, dpi = 300)

# ==============================================================================
# RANDOM FOREST SLOPE MODEL
# ==============================================================================

# The abundance predictors here are the instrument-specific z-scored versions.
# Hellinger transformation is used for the community-composition PCA only.
rf_environment_vars <- c("shelf_loc_200", "depth_mid", "Temperature", "Salinity", "Oxygen_ml_l", "chl_a_ug_l")
rf_plankton_vars <- paste0(plankton_vars, "_scaled")
predictor_vars_rf <- c(rf_environment_vars, rf_plankton_vars)

missing_rf_vars <- setdiff(predictor_vars_rf, names(dat_model_all))
if (length(missing_rf_vars) > 0L) stop("Missing RF predictor variables: ", paste(missing_rf_vars, collapse = ", "))

if (RUN_RF_SLOPE) {
  message("Fitting random-forest model for NBSS slope...")
  rf_slope_results <- fit_rf_ale(dat = dat_model_all, response = "slope", predictors = predictor_vars_rf, seed = RANDOM_SEED)
  saveRDS(rf_slope_results$model, file.path(analysis_dir, "rf_slope_caret_model.rds"), compress = TRUE)
  rf_slope_importance <- rf_slope_results$variable_importance$importance %>%
    as.data.frame() %>%
    tibble::rownames_to_column("predictor") %>%
    dplyr::arrange(dplyr::desc(Overall))
  readr::write_csv(rf_slope_importance, file.path(table_dir, "RF_slope_variable_importance.csv"))

  p_rf_imp <- ggplot2::ggplot(rf_slope_importance, ggplot2::aes(x = reorder(predictor, Overall), y = Overall)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Random-forest variable importance for NBSS slope", x = NULL, y = "Permutation importance") +
    ggplot2::theme_bw()

  ggplot2::ggsave(file.path(figure_dir, "RF_slope_variable_importance.png"), p_rf_imp, width = 8, height = 7, dpi = 300)
}

# ==============================================================================
# SAVE CORE ANALYSIS OBJECTS
# ==============================================================================

report_results <- list(
  sample_counts = sample_counts,
  pairing_loss = pairing_loss,
  shelf_location_check = shelf_location_check,
  pca_variance_table = pca_variance_table,
  pca_loadings = pca_loadings,
  pca_ci_cruise = pca_ci_cruise,
  pca_ci_transect = pca_ci_transect,
  pca_ci_tow = pca_ci_tow,
  permanova_ISIIS_table = permanova_ISIIS_table,
  permanova_ISIIS_effect = permanova_ISIIS_effect,
  permanova_project_table = permanova_project_table,
  permanova_project_effect = permanova_project_effect,
  permdisp_table = permdisp_ISIIS$permutation_table,
  permdisp_group_summary = permdisp_ISIIS$group_summary,
  n_permutations = N_PERM,
  n_community_samples = nrow(community_complete)
)

if (RUN_RF_SLOPE) report_results$rf_slope_importance <- rf_slope_importance

saveRDS(report_results, file.path(report_dir, "NBSS_community_preanalysis_results.rds"), compress = TRUE)
saveRDS(community_pca, file.path(analysis_dir, "community_PCA.rds"), compress = TRUE)

# ==============================================================================
# GENERATE R MARKDOWN REPORT
# ==============================================================================

message("Writing R Markdown report source...")

rmd_path <- file.path(report_dir, "NBSS_community_preanalysis_report.Rmd")
html_path <- file.path(report_dir, "NBSS_community_preanalysis_report.html")

rmd_lines <- c(
  "---",
  "title: \"NBSS Functional Community Composition: Instrument Pre-analysis\"",
  "output:",
  "  html_document:",
  "    self_contained: true",
  "    toc: true",
  "    toc_float: true",
  "    df_print: paged",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "library(tidyverse)",
  "res <- readRDS('NBSS_community_preanalysis_results.rds')",
  "fig_dir <- '../figures'",
  "```",
  "",
  "# Scope and interpretation",
  "",
  "This report evaluates whether Hellinger-transformed functional community composition differs systematically between ISIIS-2 and ISIIS-3 data. All available paired NBSS-abundance samples are used; no PERMANOVA or PERMDISP subsampling is performed.",
  "",
  "**Important limitation:** ISIIS version is confounded with project, sampling year, and spatial coverage. MEZCAL and MBON_ISIIS-2 comprise ISIIS-2, whereas SPECTRA and MBON_ISIIS-3 comprise ISIIS-3, and these projects do not share equivalent spatiotemporal coverage. Therefore, separation by ISIIS version cannot be interpreted as a purely causal instrument effect.",
  "",
  "# Sample coverage",
  "",
  "```{r}",
  "knitr::kable(res$sample_counts, digits = 3)",
  "```",
  "",
  "# NBSS-abundance pairing",
  "",
  "```{r}",
  "knitr::kable(res$pairing_loss, digits = 3)",
  "```",
  "",
  "# Shelf-location coverage",
  "",
  "```{r}",
  "knitr::kable(res$shelf_location_check)",
  "```",
  "",
  "# Hellinger PCA",
  "",
  "Functional-group abundances were converted to Hellinger-transformed relative community composition before PCA. PCA was centered but not z-scaled after Hellinger transformation.",
  "",
  "## Variance explained",
  "",
  "```{r}",
  "knitr::kable(res$pca_variance_table, digits = 3)",
  "```",
  "",
  "```{r, out.width='90%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_variance_scree.png'))",
  "```",
  "",
  "## PC1 and PC2 loadings",
  "",
  "```{r}",
  "knitr::kable(res$pca_loadings %>% dplyr::select(functional_group, PC1, PC2) %>% dplyr::arrange(dplyr::desc(abs(PC1))), digits = 4)",
  "```",
  "",
  "```{r, out.width='90%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_loadings_PC1_PC2.png'))",
  "```",
  "",
  "## PCA scores by instrument version",
  "",
  "The ellipses in the scatterplot are 95% fitted bivariate data ellipses, not confidence intervals for group centroids.",
  "",
  "```{r, out.width='90%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_by_ISIIS_version.png'))",
  "```",
  "",
  "## PCA scores by project",
  "",
  "```{r, out.width='90%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_by_project.png'))",
  "```",
  "",
  "# Sampling-unit centroids and 95% confidence intervals",
  "",
  "These plots show the mean PC1-PC2 position for Cruises, Cruise × Transect combinations, and Cruise × Transect × Tow combinations. Horizontal and vertical bars are marginal 95% confidence intervals for the mean PC1 and PC2 scores, respectively.",
  "",
  "## Cruise",
  "```{r, out.width='95%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_centroids_Cruise_95CI.png'))",
  "```",
  "",
  "## Cruise × Transect",
  "```{r, out.width='95%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_centroids_Transect_95CI.png'))",
  "```",
  "",
  "## Cruise × Transect × Tow",
  "```{r, out.width='95%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PCA_centroids_Tow_95CI.png'))",
  "```",
  "",
  "# Full-data PERMANOVA diagnostic",
  "",
  "Because Euclidean PERMANOVA on the Hellinger matrix is identical to the permutation test from RDA, the analysis uses `vegan::rda` + `anova.cca`. This avoids constructing the enormous all-pairs distance matrix while using every observation.",
  "",
  "## ISIIS version",
  "```{r}",
  "knitr::kable(res$permanova_ISIIS_table, digits = 6)",
  "knitr::kable(res$permanova_ISIIS_effect, digits = 6)",
  "```",
  "",
  "## Project diagnostic",
  "```{r}",
  "knitr::kable(res$permanova_project_table, digits = 6)",
  "knitr::kable(res$permanova_project_effect, digits = 6)",
  "```",
  "",
  "The R2 values are the primary effect-size diagnostics. With tens of thousands of samples, permutation p-values can be extremely small even for modest compositional differences.",
  "",
  "# Full-data PERMDISP diagnostic",
  "",
  "For Euclidean Hellinger coordinates, distance to an ISIIS-version centroid can be calculated directly in the original transformed coordinate space. The test uses the full dataset and a residual-permutation F test corresponding to the centroid form of PERMDISP, avoiding an n-by-n PCoA distance matrix.",
  "",
  "```{r}",
  "knitr::kable(res$permdisp_table, digits = 6)",
  "knitr::kable(res$permdisp_group_summary, digits = 6)",
  "```",
  "",
  "```{r, out.width='80%'}",
  "knitr::include_graphics(file.path(fig_dir, 'PERMDISP_ISIIS_version.png'))",
  "```",
  "",
  "# Random-forest slope model",
  "",
  "For the RF model only, functional-group abundances are z-scored separately within ISIIS-2 and ISIIS-3. This standardization is distinct from the Hellinger transformation used for the community PCA.",
  "",
  if (RUN_RF_SLOPE) "```{r}\nknitr::kable(res$rf_slope_importance, digits = 5)\n```" else "Random-forest fitting was disabled for this run.",
  "",
  if (RUN_RF_SLOPE) "```{r, out.width='85%'}\nknitr::include_graphics(file.path(fig_dir, 'RF_slope_variable_importance.png'))\n```" else "",
  "",
  "# Output provenance",
  "",
  paste0("Permutation count: **", N_PERM, "**."),
  "",
  paste0("Community samples analyzed: **", nrow(community_complete), "**."),
  "",
  paste0("Generated on: **", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "**.")
)

writeLines(rmd_lines, rmd_path)

if (RENDER_RMARKDOWN) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    warning("Package 'rmarkdown' is not installed. Rmd source was created but HTML was not rendered: ", rmd_path)
  } else if (!rmarkdown::pandoc_available()) {
    warning("Pandoc is not available. Rmd source was created but HTML was not rendered: ", rmd_path)
  } else {
    message("Rendering self-contained HTML report...")
    rmarkdown::render(input = rmd_path, output_file = basename(html_path), output_dir = report_dir, quiet = FALSE)
    message("HTML report: ", html_path)
  }
}

# ==============================================================================
# FINAL SESSION INFORMATION
# ==============================================================================

capture.output(sessionInfo(), file = file.path(analysis_dir, "sessionInfo.txt"))

message("DONE")
message("Results directory: ", analysis_dir)
message("R Markdown source: ", rmd_path)
if (file.exists(html_path)) message("Rendered HTML report: ", html_path)
