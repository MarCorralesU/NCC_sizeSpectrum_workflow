#!/usr/bin/env Rscript

setwd('/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(lubridate)
  library(optparse)
  library(rlang)
  library(tidyr)
  library(grid)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(spdep)
})

base_dir <- getwd()

option_list <- list(
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables",metavar = "character"),
  make_option(c("--thinned"), action = "store_true", default = FALSE, help = "Run sensitivity analysis using thinned on-shelf seston and plankton data")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

env_vars <- strsplit(opt$envVars, ",")[[1]] #env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]

sizeSpectraDirPlankton <- file.path(base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded")
sizeSpectraDirAll      <- file.path(base_dir, "data", "ISIIS-LISST_mergedNBSS_5m_allParticles")
sizeSpectraDirDetritus <- file.path(base_dir, "data", "size_spectra_detritus_5-10m-bins_thresholded")

output_dir_figs_NB <- file.path(base_dir, "figures", "size_spectra_mixed_models_cast_ID")
if (!dir.exists(output_dir_figs_NB)) {
  dir.create(output_dir_figs_NB, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_NB)
}

# specific dir for sensitivity analysis
analysis_output_dir <- if (opt$thinned) {
  file.path(output_dir_figs_NB, "revisions")
} else {
  output_dir_figs_NB
} #analysis_output_dir <-output_dir_figs_NB

dir.create(analysis_output_dir, recursive = TRUE, showWarnings = FALSE)

df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_full     <- list.files(sizeSpectraDirAll,      pattern = "LISST-ISIIS_NBSS\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_detritus <- list.files(sizeSpectraDirDetritus, pattern = "by_taxa\\.csv$", full.names = TRUE, recursive = TRUE)

load_ISIIS_data <- function(df_list,
                            env_vars = c("Temperature","Salinity","Oxygen_ml_l","pH","PAR","chl_a_ug_l"),
                            x1_full = 12.78849,
                            x2_full = 4.547442e13) {
  
  dat_ISIIS <- purrr::map2_dfr(
    df_list,
    seq_along(df_list),
    ~ readr::read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        "group_id",
        dplyr::any_of(c("total_biovolume", "biovolume_um3")),
        "biovol_um3_mid", "biovol_um3_width",
        "time_Pacific", "rel_dist", "depth_mid", "Latitude", "Longitude",
        "slope", "intercept", "Cruise", "Transect_ID", "cast_id",
        "logNB", "ECD_um_mid",
        dplyr::any_of(c("grouped_taxon", "season", "distance_sample_km")),
        dplyr::any_of(env_vars)
      )
    ) %>%
      dplyr::mutate(
        source_file = .x,
        source_id = .y
      )
  )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(time_Pacific = lubridate::ymd_hms(time_Pacific,quiet = TRUE)) %>%
    dplyr::arrange(Cruise,Transect_ID,cast_id,time_Pacific)
  
  cast_types <- dat_ISIIS %>%
    dplyr::select(Cruise, Transect_ID,cast_id,time_Pacific,depth_mid) %>%
    dplyr::filter(!is.na(Cruise),!is.na(Transect_ID),!is.na(cast_id),!is.na(time_Pacific),!is.na(depth_mid)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Cruise,Transect_ID,cast_id,time_Pacific) %>%
    dplyr::group_by(Cruise,Transect_ID,cast_id) %>%
    dplyr::summarise(
      n_depths = dplyr::n_distinct(depth_mid),
      first_depth = dplyr::first(depth_mid),
      last_depth = dplyr::last(depth_mid),
      depth_change = last_depth - first_depth,
      cast_type = dplyr::case_when(
        n_depths == 1 ~ "one_depth",
        depth_change > 0 ~ "down_cast",
        depth_change < 0 ~ "up_cast",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::left_join(cast_types %>%
    dplyr::select(Cruise, Transect_ID, cast_id, cast_type),
      by = c("Cruise","Transect_ID","cast_id"))
  
  if (!"total_biovolume" %in% names(dat_ISIIS)) dat_ISIIS$total_biovolume <- NA_real_
  if (!"biovolume_um3"  %in% names(dat_ISIIS)) dat_ISIIS$biovolume_um3  <- NA_real_
  if (!"grouped_taxon"  %in% names(dat_ISIIS)) dat_ISIIS$grouped_taxon  <- NA_character_
  if (!"season"         %in% names(dat_ISIIS)) dat_ISIIS$season         <- NA_character_
  if (!"distance_sample_km" %in% names(dat_ISIIS)) dat_ISIIS$distance_sample_km <- NA_real_
  
  dat_ISIIS <- dat_ISIIS %>%
    mutate(
      total_biovolume = coalesce(total_biovolume, biovolume_um3)
    ) %>%
    dplyr::select(-any_of("biovolume_um3")) %>%
    mutate(
      time_Pacific_orig = time_Pacific,
      mean_time_Pacific = time_Pacific - hours(1)
    ) %>%
    arrange(time_Pacific) %>%
    rename(group_id_old = group_id) %>%
    mutate(
      group_id = dplyr::cur_group_id(), .by = c(source_id, group_id_old),
      total_biovolume  = as.numeric(total_biovolume),
      biovolume_um3_m3 = as.numeric((10^logNB)*biovol_um3_width),
      biovol_um3_mid   = as.numeric(biovol_um3_mid),
      biovol_um3_width = as.numeric(biovol_um3_width),
      slope            = as.numeric(slope),
      intercept        = as.numeric(intercept),
      lower_size_edge  = biovol_um3_mid - biovol_um3_width / 2,
      upper_size_edge  = biovol_um3_mid + biovol_um3_width / 2,
      b_lin            = 10^intercept
    )
  
  group_summ <- dat_ISIIS %>%
    group_by(group_id) %>%
    summarise(
      total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
      total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE),
      min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
      max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
      slope_g = first(slope),
      intercept_g = first(intercept),
      b_lin_g = first(b_lin),
      .groups = "drop"
    ) %>%
    mutate(
      integrated_biovol_all = if_else(
        is.finite(slope_g) & is.finite(b_lin_g),
        ifelse(
          abs(slope_g + 1) < 1e-12,
          b_lin_g * log(x2_full / x1_full),
          (b_lin_g / (slope_g + 1)) * (x2_full^(slope_g + 1) - x1_full^(slope_g + 1))
        ),
        NA_real_
      ),
      integrated_biovol_sample = if_else(
        is.finite(slope_g) & is.finite(b_lin_g) &
          is.finite(min_size) & is.finite(max_size) &
          min_size > 0 & max_size > 0 & max_size >= min_size,
        ifelse(
          abs(slope_g + 1) < 1e-12,
          b_lin_g * log(max_size / min_size),
          (b_lin_g / (slope_g + 1)) * (max_size^(slope_g + 1) - min_size^(slope_g + 1))
        ),
        NA_real_
      )
    ) %>%
    dplyr::select(
      group_id,
      total_biovolume_sum,
      total_biovolume_um3_m3,
      integrated_biovol_all,
      integrated_biovol_sample
    )
  
  dat_ISIIS <- dat_ISIIS %>%
    left_join(group_summ, by = "group_id") %>%
    mutate(
      slope_anomaly = slope + 1,
      NB = 10^(logNB),
      total_biovolume = total_biovolume_sum
    )
  
  return(dat_ISIIS)
}

dat_ISIIS_full_raw     <- load_ISIIS_data(df_list_full, env_vars = env_vars)
dat_ISIIS_plankton_raw <- load_ISIIS_data(df_list_plankton, env_vars = env_vars)
dat_ISIIS_detritus_raw <- load_ISIIS_data(df_list_detritus, env_vars = env_vars)

remove_duplicate_NBSS_samples <- function(dat, verbose = TRUE) {
  
  # ------------------------------------------------------------
  # 1. Create sample-level dataset
  # ------------------------------------------------------------
  
  sample_data <- dat %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(
      Cruise = dplyr::first(Cruise),
      Transect_ID = dplyr::first(Transect_ID),
      cast_id = dplyr::first(cast_id),
      depth_mid = dplyr::first(depth_mid),
      Longitude = dplyr::first(Longitude),
      time_Pacific = dplyr::first(time_Pacific),
      slope = dplyr::first(slope),
      intercept = dplyr::first(intercept),
      total_biovolume_um3_m3 = dplyr::first(total_biovolume_um3_m3),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------
  # 2. Identify duplicated NBSS samples
  # ------------------------------------------------------------
  
  duplicate_samples <- sample_data %>%
    dplyr::filter(
      !is.na(depth_mid),
      !is.na(Longitude),
      !is.na(slope),
      !is.na(intercept),
      !is.na(total_biovolume_um3_m3)
    ) %>%
    dplyr::group_by(
      Cruise, Transect_ID, depth_mid, Longitude,
      slope, intercept, total_biovolume_um3_m3
    ) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::arrange(time_Pacific, group_id, .by_group = TRUE) %>%
    dplyr::mutate(
      duplicate_rank = dplyr::row_number(),
      action = dplyr::if_else(
        duplicate_rank == 1,
        "retain",
        "remove"
      )
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------
  # 3. Identify group_ids to remove
  # ------------------------------------------------------------
  
  groups_to_remove <- duplicate_samples %>%
    dplyr::filter(action == "remove") %>%
    dplyr::pull(group_id)
  
  # ------------------------------------------------------------
  # 4. Remove all size-bin rows for duplicated samples
  # ------------------------------------------------------------
  
  dat_clean <- dat %>%
    dplyr::filter(!group_id %in% groups_to_remove)
  
  # ------------------------------------------------------------
  # 5. Print summary
  # ------------------------------------------------------------
  
  if (verbose) {
    cat("\n========================================\n")
    cat("NBSS DUPLICATE SAMPLE CLEANING\n")
    cat("========================================\n")
    cat("Original samples:", dplyr::n_distinct(dat$group_id), "\n")
    cat("Duplicate sets:", sum(duplicate_samples$duplicate_rank == 1), "\n")
    cat("Samples removed:", length(groups_to_remove), "\n")
    cat("Remaining samples:", dplyr::n_distinct(dat_clean$group_id), "\n")
    cat("========================================\n\n")
  }
  
  return(list(
    data = dat_clean,
    duplicate_report = duplicate_samples,
    removed_group_ids = groups_to_remove
  ))
}

Seston_cleaning <- remove_duplicate_NBSS_samples(dat_ISIIS_full_raw)
dat_ISIIS_full_raw <- Seston_cleaning$data

plankton_cleaning <- remove_duplicate_NBSS_samples(dat_ISIIS_plankton_raw)
dat_ISIIS_plankton_raw <- plankton_cleaning$data

detritus_cleaning <- remove_duplicate_NBSS_samples(dat_ISIIS_detritus_raw)
dat_ISIIS_detritus_raw <- detritus_cleaning$data



clean_ISIIS_dat <- function(dat_ISIIS) {
  dat_ISIIS <- dat_ISIIS %>%
    mutate(
      shelf_loc_200 = case_when(
        rel_dist >= 0 ~ "on_shelf",
        rel_dist <  0 ~ "off_shelf",
        TRUE          ~ NA_character_
      ),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf")),
      Cruise        = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      Transect_ID   = factor(Transect_ID, levels = c("GH", "OB", "CR", "CM", "NH", "HH", "RR")),
      season        = factor(season, levels = c("Winter", "Summer"))
    ) %>%
    filter(!(Cruise == "w22" & Transect_ID == "CM" & Longitude < -124.75))
  
  req_cols <- c("season", "distance_sample_km")
  miss <- req_cols[!req_cols %in% names(dat_ISIIS)]
  if (length(miss) > 0) {
    stop("Missing required columns for mixed models: ", paste(miss, collapse = ", "))
  }
  
  if (all(is.na(dat_ISIIS$season))) {
    stop("Column 'season' exists but is all NA.")
  }
  if (all(is.na(dat_ISIIS$distance_sample_km))) {
    stop("Column 'distance_sample_km' exists but is all NA.")
  }
  
  dat_ISIIS
}

dat_ISIIS_full     <- clean_ISIIS_dat(dat_ISIIS_full_raw)
dat_ISIIS_plankton <- clean_ISIIS_dat(dat_ISIIS_plankton_raw)
dat_ISIIS_detritus <- clean_ISIIS_dat(dat_ISIIS_detritus_raw)


# create diagnostic tables to show unbalanced design between inshore and offshore

cast_diagnostics_seston <- dat_ISIIS_full %>%
  dplyr::filter(shelf_loc_200 %in% c("on_shelf", "off_shelf")) %>%
  dplyr::distinct(Cruise, Transect_ID, group_id, shelf_loc_200, cast_type) %>%
  dplyr::count(Cruise, Transect_ID, shelf_loc_200, cast_type, name = "n_samples") %>%
  tidyr::complete(Cruise, Transect_ID, shelf_loc_200, cast_type = c( "down_cast", "up_cast", "one_depth"), fill = list(n_samples = 0)) %>%
  tidyr::pivot_wider(names_from = cast_type, values_from = n_samples, values_fill = 0) %>%
  dplyr::mutate(total_samples = down_cast + up_cast + one_depth) %>%
  dplyr::select(Cruise, Transect_ID, shelf_loc_200, total_samples, down_cast, up_cast, one_depth) %>%
  dplyr::arrange(Cruise, Transect_ID, shelf_loc_200)

cast_diagnostics_plankton <- dat_ISIIS_plankton %>%
  dplyr::filter(shelf_loc_200 %in% c("on_shelf", "off_shelf")) %>%
  dplyr::distinct(Cruise, Transect_ID, group_id, shelf_loc_200, cast_type) %>%
  dplyr::count(Cruise, Transect_ID, shelf_loc_200, cast_type, name = "n_samples") %>%
  tidyr::complete(Cruise, Transect_ID, shelf_loc_200, cast_type = c( "down_cast", "up_cast", "one_depth"), fill = list(n_samples = 0)) %>%
  tidyr::pivot_wider(names_from = cast_type, values_from = n_samples, values_fill = 0) %>%
  dplyr::mutate(total_samples = down_cast + up_cast + one_depth) %>%
  dplyr::select(Cruise, Transect_ID, shelf_loc_200, total_samples, down_cast, up_cast, one_depth) %>%
  dplyr::arrange(Cruise, Transect_ID, shelf_loc_200)

cast_diagnostics_pct_seston <- cast_diagnostics_seston %>%
  dplyr::mutate(
    pct_down_cast = 100 * down_cast / total_samples,
    pct_up_cast = 100 * up_cast / total_samples,
    pct_one_depth = 100 * one_depth / total_samples
  )

cast_diagnostics_pct_plankton <- cast_diagnostics_plankton %>%
  dplyr::mutate(
    pct_down_cast = 100 * down_cast / total_samples,
    pct_up_cast = 100 * up_cast / total_samples,
    pct_one_depth = 100 * one_depth / total_samples
  )

dir.create(file.path(output_dir_figs_NB, "revisions"), recursive = TRUE, showWarnings = FALSE)

sample_size_dir <- file.path(output_dir_figs_NB, "revisions", "sample_sizes")
dir.create(sample_size_dir, recursive = TRUE, showWarnings = FALSE)





# new process for graphing sample sizes 
# ============================================================
# Diagnostic figures: number of samples on/off shelf
# ============================================================

transect_levels <- c("GH","CR","CM","NH","HH","RR")

cruise_levels <- c("w22", "s22", "w23", "s23")

cruise_label_map <- c( w22 = "Winter 2022", s22 = "Summer 2022", w23 = "Winter 2023", s23 = "Summer 2023")

transect_label_map <- c(GH = "Grays Harbor", CR = "Columbia River", CM = "Cape Meares",  NH = "Newport Line", HH = "Heceta Head", RR = "Rogue River")

season_colors <- c("Winter" = "#0072B2", "Summer" = "#E69F00")

make_cast_diagnostic_figure <- function(dat, figure_title) {
  
  plot_dat <- dat %>%
    dplyr::filter(Cruise %in% cruise_levels, Transect_ID %in% transect_levels, shelf_loc_200 %in% c("off_shelf", "on_shelf")) %>%
    dplyr::mutate(
      season = dplyr::case_when(
        Cruise %in% c("w22", "w23") ~ "Winter",
        Cruise %in% c("s22", "s23") ~ "Summer",
        TRUE ~ NA_character_
      ),
      season = factor(season, levels = c("Winter", "Summer")),
      Cruise = factor(Cruise, levels = cruise_levels, labels = cruise_label_map[cruise_levels]),
      Transect_ID = factor(Transect_ID, levels = transect_levels, labels = transect_label_map[transect_levels]),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"), labels = c("Off shelf","On shelf"))
    )
  
  ggplot2::ggplot(plot_dat, ggplot2::aes(x = shelf_loc_200, y = total_samples, fill = season)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = total_samples), vjust = -0.4, size = 4) +
    ggplot2::facet_grid(rows = ggplot2::vars(Transect_ID), cols = ggplot2::vars(Cruise), drop = FALSE, switch = "y") +
    ggplot2::scale_fill_manual(values = season_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = NULL, y = "Number of samples", title = figure_title, fill = "Season") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 13, face = "bold"),
      strip.text.y.left = ggplot2::element_text(size = 12, face = "bold", angle = 90),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 13, face = "bold"),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      legend.position = "top",
      panel.spacing = grid::unit(0.6, "lines")
    )
}

p_cast_diagnostics_seston <- make_cast_diagnostic_figure(
  cast_diagnostics_seston,
  "Seston"
)

p_cast_diagnostics_plankton <- make_cast_diagnostic_figure(
  cast_diagnostics_plankton,
  "Plankton"
)

ggplot2::ggsave(file.path(sample_size_dir, "SPECTRA_sample_sizes_seston.png"), p_cast_diagnostics_seston, width = 14, height = 16, dpi = 300, bg = "white"
)

ggplot2::ggsave(file.path(sample_size_dir, "SPECTRA_sample_sizes_plankton.png"), p_cast_diagnostics_plankton, width = 14, height = 16, dpi = 300, bg = "white")

#counts of season X shelf location
cast_counts_seston <- cast_diagnostics_seston %>%
  dplyr::group_by(Cruise,shelf_loc_200) %>%
  dplyr::summarise(total_samples = sum(total_samples),.groups = "drop") %>%
  dplyr::arrange(factor(Cruise, levels = cruise_levels),factor(shelf_loc_200,levels = c("off_shelf","on_shelf")))

cast_counts_plankton <- cast_diagnostics_plankton %>%
  dplyr::group_by(Cruise,shelf_loc_200) %>%
  dplyr::summarise(total_samples = sum(total_samples),.groups = "drop") %>%
  dplyr::arrange(factor(Cruise, levels = cruise_levels),factor(shelf_loc_200,levels = c("off_shelf","on_shelf")))


readr::write_csv(cast_diagnostics_seston, file.path(sample_size_dir, "SPECTRA_sample_sizes_byCruise_Transect_seston.csv"))

readr::write_csv(cast_diagnostics_plankton, file.path(sample_size_dir, "SPECTRA_sample_sizes_byCruise_Transect_plankton.csv"))

readr::write_csv(cast_counts_seston, file.path(sample_size_dir, "SPECTRA_sample_sizes_byCruise_seston.csv"))

readr::write_csv(cast_counts_plankton, file.path(sample_size_dir, "SPECTRA_sample_sizes_byCruise_plankton.csv"))







# Systematically thin on-shelf casts by retaining every third cast

identify_coarse_casts <- function(x) {
  
  x <- x %>%
    dplyr::arrange(time_Pacific) %>%
    dplyr::mutate(cast_position = dplyr::row_number(),
                  keep_cast = (cast_position - 1) %% 3 == 0)
  
  x
}
thin_onshelf_casts <- function(dat) {
  
  cast_sequence <- dat %>%
    dplyr::filter(shelf_loc_200 == "on_shelf", cast_type %in% c("down_cast", "up_cast")) %>%
    dplyr::group_by(Cruise, Transect_ID, cast_id, cast_type) %>%
    dplyr::summarise(time_Pacific = min(time_Pacific, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Cruise, Transect_ID, time_Pacific) %>%
    dplyr::group_by(Cruise, Transect_ID) %>%
    dplyr::group_modify(~ identify_coarse_casts(.x)) %>%
    dplyr::ungroup()
  
  casts_to_keep <- cast_sequence %>%
    dplyr::filter(keep_cast) %>%
    dplyr::select(Cruise, Transect_ID, cast_id) %>%
    dplyr::distinct()
  
  dat_thinned <- dat %>%
    dplyr::filter(cast_type != "one_depth") %>%
    dplyr::filter(shelf_loc_200 != "on_shelf" |
                    paste(Cruise, Transect_ID, cast_id) %in%
                    paste(casts_to_keep$Cruise, casts_to_keep$Transect_ID, casts_to_keep$cast_id))
  
  return(list(
    data = dat_thinned,
    cast_sequence = cast_sequence
  ))
}

#apply function
seston_thinning <- thin_onshelf_casts(dat_ISIIS_full)
plankton_thinning <- thin_onshelf_casts(dat_ISIIS_plankton)

#extract datasets
dat_ISIIS_full_thinned <- seston_thinning$data
dat_ISIIS_plankton_thinned <- plankton_thinning$data

#extract diagnostics
cast_sequence_seston <- seston_thinning$cast_sequence
cast_sequence_plankton <- plankton_thinning$cast_sequence

readr::write_csv(cast_sequence_seston, file.path(output_dir_figs_NB, "revisions", "SPECTRA_cast_subsettingReport_seston.csv"))

readr::write_csv(cast_sequence_plankton, file.path(output_dir_figs_NB, "revisions", "SPECTRA_cast_subsettingReport_plankton.csv"))


#ANALYSIS STARTS HERE: first, we decide whether the subsetted or the full dataset
dat_ISIIS_full_analysis <- if (opt$thinned) dat_ISIIS_full_thinned else dat_ISIIS_full
dat_ISIIS_plankton_analysis <- if (opt$thinned) dat_ISIIS_plankton_thinned else dat_ISIIS_plankton


#function added /10/05/2026, determine if there is spatial autocorrelation

diagnose_mixed_model_residuals <- function(model, df_mod, model_prefix, output_dir,
                                           horizontal_breaks = seq(0, 100, by = 5),
                                           vertical_breaks = seq(0, 100, by = 5),
                                           max_neighbors = 20,
                                           min_pairs = 10) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------
  # 1. Extract observations actually used by the fitted model
  # ------------------------------------------------------------
  
  model_rows <- suppressWarnings(as.integer(rownames(stats::model.frame(model))))
  
  if (anyNA(model_rows) || length(model_rows) != stats::nobs(model)) {
    stop("Could not align model observations with df_mod.")
  }
  
  diag_dat <- df_mod[model_rows, , drop = FALSE] %>%
    dplyr::mutate(
      residual = as.numeric(stats::residuals(model)),
      fitted_value = as.numeric(stats::fitted(model)),
      standardized_residual = residual / stats::sigma(model)
    ) %>%
    dplyr::filter(
      is.finite(standardized_residual),
      is.finite(fitted_value)
    )
  
  if (!"cast_id" %in% names(diag_dat)) {
    stop("cast_id is required for vertical residual diagnostics.")
  }
  
  # ------------------------------------------------------------
  # 2. Calculate neighboring residual pairs without constructing
  #    every possible pair of observations
  # ------------------------------------------------------------
  
  make_lag_pairs <- function(dat, coordinate, breaks, direction) {
    
    dat <- dat[is.finite(dat[[coordinate]]), , drop = FALSE]
    dat <- dat[order(dat[[coordinate]]), , drop = FALSE]
    
    n <- nrow(dat)
    
    if (n < 2) {
      return(tibble::tibble())
    }
    
    pair_list <- vector("list", min(max_neighbors, n - 1))
    
    for (k in seq_len(min(max_neighbors, n - 1))) {
      
      i <- seq_len(n - k)
      j <- i + k
      
      pair_list[[k]] <- tibble::tibble(
        lag = abs(dat[[coordinate]][j] - dat[[coordinate]][i]),
        residual_1 = dat$standardized_residual[i],
        residual_2 = dat$standardized_residual[j]
      )
    }
    
    dplyr::bind_rows(pair_list) %>%
      dplyr::filter(
        is.finite(lag),
        lag > 0,
        lag <= max(breaks)
      ) %>%
      dplyr::mutate(
        direction = direction,
        lag_bin = cut(
          lag,
          breaks = breaks,
          include.lowest = TRUE,
          right = FALSE),
        lag_mid = ((breaks[-1] + breaks[-length(breaks)]) / 2)[as.integer(lag_bin)]
      )%>%
      dplyr::filter(!is.na(lag_bin))
  }
  
  # ------------------------------------------------------------
  # 3. Horizontal autocorrelation
  #    Compare observations at the same depth within transects
  # ------------------------------------------------------------
  
  horizontal_pairs <- diag_dat %>%
    dplyr::filter(
      is.finite(distance_sample_km),
      is.finite(depth_mid)
    ) %>%
    dplyr::group_by(Cruise, Transect_ID, depth_mid) %>%
    dplyr::group_modify(
      ~ make_lag_pairs(
        .x,
        coordinate = "distance_sample_km",
        breaks = horizontal_breaks,
        direction = "Horizontal"
      )
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------
  # 4. Vertical autocorrelation
  #    Compare depths within individual casts
  # ------------------------------------------------------------
  
  vertical_pairs <- diag_dat %>%
    dplyr::filter(
      !is.na(cast_id),
      is.finite(depth_mid)
    ) %>%
    dplyr::group_by(Cruise, Transect_ID, cast_id) %>%
    dplyr::group_modify(
      ~ make_lag_pairs(
        .x,
        coordinate = "depth_mid",
        breaks = vertical_breaks,
        direction = "Vertical"
      )
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------
  # 5. Summarize residual correlations by separation distance
  # ------------------------------------------------------------
  
  correlation_data <- dplyr::bind_rows(
    horizontal_pairs,
    vertical_pairs
  ) %>%
    dplyr::group_by(direction, lag_bin, lag_mid) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      correlation = if (
        dplyr::n() >= min_pairs &&
        stats::sd(residual_1) > 0 &&
        stats::sd(residual_2) > 0
      ) {
        stats::cor(residual_1, residual_2)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  
  readr::write_csv(
    correlation_data,
    file.path(output_dir, paste0(model_prefix, "_residual_correlations.csv"))
  )
  
  # Summarize short-distance residual autocorrelation
  
  short_lag_summary <- correlation_data %>%
    dplyr::filter(
      (direction == "Horizontal" & lag_mid == 2.5) |
        (direction == "Vertical" & lag_mid == 2.5)
    ) %>%
    dplyr::select(
      direction,
      lag_mid,
      n_pairs,
      correlation
    )
  
  cat("\nSHORT-DISTANCE RESIDUAL AUTOCORRELATION\n")
  cat("Model:", model_prefix, "\n")
  print(short_lag_summary)
  cat("========================================\n\n")
  
  # ------------------------------------------------------------
  # 6. Panel A: Horizontal residual correlation
  # ------------------------------------------------------------
  
  p_horizontal <- correlation_data %>%
    dplyr::filter(direction == "Horizontal") %>%
    ggplot2::ggplot(
      ggplot2::aes(x = lag_mid, y = correlation)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_point(color = "#0072B2", size = 2) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "A. Horizontal residual correlation",
      x = "Along-transect separation (km)",
      y = "Pearson correlation"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 7. Panel B: Vertical residual correlation
  # ------------------------------------------------------------
  
  p_vertical <- correlation_data %>%
    dplyr::filter(direction == "Vertical") %>%
    ggplot2::ggplot(
      ggplot2::aes(x = lag_mid, y = correlation)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_line(color = "#D55E00", linewidth = 0.8) +
    ggplot2::geom_point(color = "#D55E00", size = 2) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "B. Vertical residual correlation",
      x = "Depth separation (m)",
      y = "Pearson correlation"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 8. Panel C: Standardized residuals versus fitted values
  # ------------------------------------------------------------
  
  p_fitted <- ggplot2::ggplot(
    diag_dat,
    ggplot2::aes(x = fitted_value, y = standardized_residual)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_point(
      alpha = 0.2, size = 0.7, color = "#0072B2"
    ) +
    ggplot2::geom_smooth(
      method = "loess", se = FALSE, color = "#D55E00",
      linewidth = 0.9
    ) +
    ggplot2::labs(
      title = "C. Residuals vs. fitted values",
      x = "Fitted values",
      y = "Standardized residuals"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 9. Panel D: Normal Q-Q plot
  # ------------------------------------------------------------
  
  p_qq <- ggplot2::ggplot(
    diag_dat,
    ggplot2::aes(sample = standardized_residual)
  ) +
    ggplot2::stat_qq(
      alpha = 0.35, size = 0.8, color = "#0072B2"
    ) +
    ggplot2::stat_qq_line(
      color = "#D55E00", linewidth = 0.9
    ) +
    ggplot2::labs(
      title = "D. Normal Q-Q plot",
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 10. Assemble and save the four-panel diagnostic figure
  # ------------------------------------------------------------
  
  p_diagnostics <- (p_horizontal | p_vertical) /
    (p_fitted | p_qq)
  
  p_diagnostics <- p_diagnostics +
    patchwork::plot_annotation(
      title = paste("Mixed-effects model diagnostics:", model_prefix)
    )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(model_prefix, "_residual_diagnostics_4panel.png")
    ),
    plot = p_diagnostics,
    width = 13,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  # ------------------------------------------------------------
  # 11. Export model residuals and additional diagnostics
  # ------------------------------------------------------------
  
  readr::write_csv(
    diag_dat,
    file.path(output_dir, paste0(model_prefix, "_model_residuals.csv"))
  )
  
  convergence_messages <- model@optinfo$conv$lme4$messages
  
  model_checks <- tibble::tibble(
    n_observations = stats::nobs(model),
    residual_sd = stats::sigma(model),
    singular_fit = lme4::isSingular(model),
    convergence_messages = if (length(convergence_messages) == 0) {
      "None"
    } else {
      paste(convergence_messages, collapse = "; ")
    }
  )
  
  readr::write_csv(
    model_checks,
    file.path(output_dir, paste0(model_prefix, "_model_checks.csv"))
  )
  
  message("Saved residual diagnostics for: ", model_prefix)
  
  invisible(list(
    residual_data = diag_dat,
    correlations = correlation_data,
    model_checks = model_checks,
    plot = p_diagnostics
  ))
}


#then, the anaylsis follows

PSS_violin_plus_mixedModel <- function(dat,
                                       string_ID,
                                       PSS_var = "slope",
                                       depth_cut = 100,
                                       log10_transform = FALSE,
                                       output_dir = output_dir_figs_NB,
                                       y_limits = NULL,
                                       random_effect = "(1 | Cruise/Transect_ID/cast_id)",
                                       include_covariates = TRUE,
                                       season_colors = c("Winter" = "#0072B2",
                                                         "Summer" = "#E69F00"),
                                       gap = 0.35,
                                       ci_cap_width = 0.03,
                                       width = 12,
                                       height = 6,
                                       dpi = 300,
                                       save_plot = TRUE,
                                       save_model = TRUE,
                                       model_dir = output_dir_figs_NB,
                                       model_prefix = NULL) {
  
  needed <- c(
    PSS_var, "season", "shelf_loc_200", "depth_mid",
    "Cruise", "Transect_ID", "cast_id", "group_id",
    "distance_sample_km"
  )
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))
  
  yvar <- sym(PSS_var)
  
  df <- dat %>%
    filter(.data$depth_mid < depth_cut) %>%
    mutate(
      season = factor(.data$season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      Cruise = factor(.data$Cruise),
      Transect_ID = factor(.data$Transect_ID),
      .y_raw = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    filter(!is.na(.data$season), !is.na(.data$shelf_loc_200), is.finite(.data$.y_raw))
  
  if (nrow(df) == 0) stop("After filtering depth/NA/non-finite, no rows remain.")
  
  if (log10_transform) {
    df <- df %>%
      mutate(.y = if_else(.data$.y_raw > 0, log10(.data$.y_raw), NA_real_)) %>%
      filter(is.finite(.data$.y))
  } else {
    df <- df %>% mutate(.y = .data$.y_raw)
  }
  
  if (PSS_var %in% c("total_biovolume", "total_biovolume_um3_m3", "integrated_biovol_all", "integrated_biovol_sample")) {
    y_label_expr <- expression(Biovolume~(mu*m^3~m^-3))
  } else if (PSS_var == "intercept") {
    y_label_expr <- expression(Intercept~(mu*m^3~m^-3~mu*m^-3))
  } else if (PSS_var == "slope") {
    y_label_expr <- expression(Slope~(m^-3~mu*m^-3))
  } else {
    y_label_expr <- PSS_var
  }
  
  y_scale_obj <- if (log10_transform || PSS_var == "intercept") {
    ggplot2::scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", round(x, 2)))
    )
  } else {
    NULL
  }
  
  cast_consistency <- df %>%
    dplyr::group_by(Cruise, Transect_ID, group_id) %>%
    dplyr::summarise(
      n_cast_ids = dplyr::n_distinct(cast_id, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cast_ids != 1)
  
  cat(
    "\nNumber of group_ids without exactly one valid cast_id:",
    nrow(cast_consistency),
    "\n"
  )
  
  df_mod <- df %>%
    dplyr::group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    dplyr::summarise(
      y = dplyr::first(.y),
      depth_mid = dplyr::first(depth_mid),
      distance_sample_km = dplyr::first(distance_sample_km),
      cast_id = dplyr::first(cast_id),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(y))
  
  df_mod <- df_mod %>%
    dplyr::filter(!is.na(cast_id)) %>%
    dplyr::mutate(cast_id = factor(cast_id))
  
  cat("\nCAST-LEVEL MODEL STRUCTURE\n")
  cat("Observations:", nrow(df_mod), "\n")
  cat(
    "Unique casts:",
    dplyr::n_distinct(
      interaction(df_mod$Cruise, df_mod$Transect_ID, df_mod$cast_id, drop = TRUE)
    ),
    "\n"
  )

  
  cat("\n========================================\n")
  cat("MODEL DATA DEPTH RANGE\n")
  cat("Variable:", PSS_var, "| Dataset:", string_ID, "\n")
  cat("Requested depth cutoff:", depth_cut, "m\n")
  cat("df_mod depth range (m):\n")
  print(range(df_mod$depth_mid, na.rm = TRUE))
  cat("========================================\n\n")
  
  
  fixed_part <- if (include_covariates) {
    "y ~ season * shelf_loc_200 + scale(depth_mid) + scale(distance_sample_km)"
  } else {
    "y ~ season * shelf_loc_200"
  }
  fml <- stats::as.formula(paste(fixed_part, "+", random_effect))
  
  m1 <- lmerTest::lmer(fml, data = df_mod, REML = TRUE)
  
  at_list <- if (include_covariates) {
    list(`scale(depth_mid)` = 0, `scale(distance_sample_km)` = 0)
  } else {
    list()
  }
  
  emm <- emmeans::emmeans(m1, ~ season * shelf_loc_200, at = at_list)
  pred_df <- as.data.frame(emm)
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$asymp.LCL, upper.CL = .data$asymp.UCL)
  } else if (all(c("LCL", "UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$LCL, upper.CL = .data$UCL)
  }
  
  dodge_w <- 0.75
  
  p_violin <- ggplot(df, aes(x = .data$shelf_loc_200, y = .data$.y)) +
    geom_violin(
      aes(fill = .data$season),
      position = position_dodge(width = dodge_w),
      trim = FALSE,
      alpha = 0.65,
      linewidth = 0.3
    ) +
    geom_boxplot(
      aes(group = interaction(.data$shelf_loc_200, .data$season)),
      position = position_dodge(width = dodge_w),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    scale_x_discrete(
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      )
    ) +
    scale_fill_manual(values = season_colors, name = "Season") +
    labs(
      title = paste0("Raw data (depth < ", depth_cut, " m) — ", string_ID),
      x = NULL,
      y = y_label_expr
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  pred_df2 <- pred_df %>%
    mutate(
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      season = factor(.data$season, levels = c("Winter", "Summer")),
      x_base_comp = 1 + (as.numeric(.data$shelf_loc_200) - 1) * gap,
      x = .data$x_base_comp + ifelse(.data$season == "Winter", -0.12, 0.12)
    )
  
  seg_df <- pred_df2 %>%
    dplyr::select(.data$shelf_loc_200, .data$season, .data$x, .data$emmean) %>%
    pivot_wider(names_from = .data$season, values_from = c(.data$x, .data$emmean)) %>%
    transmute(
      shelf_loc_200,
      x = .data$x_Winter,
      xend = .data$x_Summer,
      y = .data$emmean_Winter,
      yend = .data$emmean_Summer
    )
  
  p_model <- ggplot() +
    geom_segment(
      data = seg_df,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      color = "black", linetype = "dashed", linewidth = 0.9
    ) +
    geom_errorbar(
      data = pred_df2,
      aes(x = .data$x, y = .data$emmean, ymin = .data$lower.CL, ymax = .data$upper.CL, color = .data$season),
      width = ci_cap_width, linewidth = 0.6
    ) +
    geom_point(
      data = pred_df2,
      aes(x = .data$x, y = .data$emmean, color = .data$season),
      size = 3
    ) +
    scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("off shelf", "on shelf")
    )+
    scale_color_manual(values = season_colors, name = "Season") +
    labs(
      title = "Mixed model (emmeans ± CI)",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  if (is.null(y_limits)) {
    y_limits <- range(df$.y, pred_df2$lower.CL, pred_df2$upper.CL, na.rm = TRUE)
  }
  
  p_violin <- p_violin + coord_cartesian(ylim = y_limits)
  p_model  <- p_model  + coord_cartesian(ylim = y_limits)
  
  if (!is.null(y_scale_obj)) {
    p_violin <- p_violin + y_scale_obj
    p_model  <- p_model  + y_scale_obj
  }
  
  p_combined <- p_violin + p_model + patchwork::plot_layout(widths = c(1.3, 1))
  print(p_combined)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_", PSS_var,
             if (log10_transform) "_log10" else "",
             "_violin_plusMixedModel_depthlt", depth_cut, "_", string_ID, ".png")
    )
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p_combined, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved plot to: ", out_path)
  }
  
  if (save_model) {
    if (is.null(model_prefix)) {
      model_prefix <- paste0(
        "SPECTRA_", PSS_var,
        if (log10_transform) "_log10" else "",
        "_mixedModel_depthlt", depth_cut, "_", string_ID
      )
    }
    
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    
    
    model_diagnostics <- diagnose_mixed_model_residuals(
      model = m1,
      df_mod = df_mod,
      model_prefix = model_prefix,
      output_dir = model_dir
    )
    
    
    cat("\nRESIDUAL DIAGNOSTICS DEPTH RANGE\n")
    cat("Model:", model_prefix, "\n")
    
    print(
      dplyr::summarise(
        model_diagnostics$residual_data,
        min_depth = min(depth_mid, na.rm = TRUE),
        max_depth = max(depth_mid, na.rm = TRUE)
      )
    )
    
    cat("========================================\n\n")
    
    
    cat("\nVERTICAL RESIDUAL CORRELATION LAG RANGE\n")
    cat("Model:", model_prefix, "\n")
    
    vertical_lags <- model_diagnostics$correlations %>%
      dplyr::filter(direction == "Vertical", is.finite(lag_mid))
    
    if (nrow(vertical_lags) > 0) {
      print(
        vertical_lags %>%
          dplyr::summarise(
            min_lag = min(lag_mid),
            max_lag = max(lag_mid)
          )
      )
    } else {
      cat("No populated vertical lag bins were available.\n")
    }
    
    cat("========================================\n\n")
    
    rds_path <- file.path(model_dir, paste0(model_prefix, ".rds"))
    save_obj <- list(
      model = m1,
      formula = fml,
      df_mod = df_mod,
      emmeans_raw = pred_df,
      emmeans_plot = pred_df2,
      seg_df = seg_df
    )
    saveRDS(save_obj, rds_path)
    
    sum_path   <- file.path(model_dir, paste0(model_prefix, "_summary.txt"))
    anova_path <- file.path(model_dir, paste0(model_prefix, "_anova.txt"))
    emm_path   <- file.path(model_dir, paste0(model_prefix, "_emmeans.csv"))
    
    capture.output(summary(m1), file = sum_path)
    capture.output(anova(m1), file = anova_path)
    utils::write.csv(pred_df, emm_path, row.names = FALSE)
  }
  
  invisible(list(
    plot = p_combined,
    model = m1,
    emmeans = pred_df2,
    df_mod = df_mod
  ))
}

run_mm_set <- function(dat, string_ID, output_dir) {
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "slope",
    depth_cut = 100,
    output_dir = output_dir,
    model_dir = output_dir,
    y_limits = c(-2.5, 0),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "intercept",
    depth_cut = 100,
    output_dir = output_dir,
    model_dir = output_dir,
    y_limits = c(2, 20),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "total_biovolume_um3_m3",
    depth_cut = 100,
    log10_transform = TRUE,
    output_dir = output_dir,
    model_dir = output_dir,
    y_limits = c(7, 19),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "integrated_biovol_all",
    depth_cut = 100,
    log10_transform = TRUE,
    output_dir = output_dir,
    model_dir = output_dir,
    y_limits = c(7, 19),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "integrated_biovol_sample",
    depth_cut = 100,
    log10_transform = TRUE,
    output_dir = output_dir,
    model_dir = output_dir,
    y_limits = c(7, 19),
    save_plot = TRUE,
    save_model = TRUE
  )
}


prep_mm_for_compare <- function(dat,
                                string_ID,
                                PSS_var = "slope",
                                depth_cut = 100,
                                log10_transform = FALSE,
                                random_effect = "(1 | Cruise/Transect_ID/cast_id)",
                                include_covariates = TRUE,
                                gap = 0.35) {
  
  needed <- c(PSS_var, "season", "shelf_loc_200", "depth_mid",
    "Cruise", "Transect_ID", "cast_id", "group_id",
    "distance_sample_km")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))
  
  yvar <- sym(PSS_var)
  
  df <- dat %>%
    filter(.data$depth_mid < depth_cut) %>%
    mutate(
      season = factor(.data$season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      Cruise = factor(.data$Cruise),
      Transect_ID = factor(.data$Transect_ID),
      .y_raw = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    filter(!is.na(.data$season), !is.na(.data$shelf_loc_200), is.finite(.data$.y_raw))
  
  if (nrow(df) == 0) stop("After filtering depth/NA/non-finite, no rows remain.")
  
  if (log10_transform) {
    df <- df %>%
      mutate(.y = if_else(.data$.y_raw > 0, log10(.data$.y_raw), NA_real_)) %>%
      filter(is.finite(.data$.y))
  } else {
    df <- df %>% mutate(.y = .data$.y_raw)
  }
  
  
  df_mod <- df %>%
    dplyr::group_by(
      Cruise, Transect_ID, group_id,
      season, shelf_loc_200
    ) %>%
    dplyr::summarise(
      y = dplyr::first(.y),
      depth_mid = dplyr::first(depth_mid),
      distance_sample_km = dplyr::first(distance_sample_km),
      cast_id = dplyr::first(cast_id),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      is.finite(y),
      !is.na(cast_id)
    ) %>%
    dplyr::mutate(cast_id = factor(cast_id))
  
  
  fixed_part <- if (include_covariates) {
    "y ~ season * shelf_loc_200 + scale(depth_mid) + scale(distance_sample_km)"
  } else {
    "y ~ season * shelf_loc_200"
  }
  fml <- stats::as.formula(paste(fixed_part, "+", random_effect))
  
  m1 <- lmerTest::lmer(fml, data = df_mod, REML = TRUE)
  
  at_list <- if (include_covariates) {
    list(`scale(depth_mid)` = 0, `scale(distance_sample_km)` = 0)
  } else {
    list()
  }
  
  emm <- emmeans::emmeans(m1, ~ season * shelf_loc_200, at = at_list)
  pred_df <- as.data.frame(emm)
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$asymp.LCL, upper.CL = .data$asymp.UCL)
  } else if (all(c("LCL", "UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$LCL, upper.CL = .data$UCL)
  }
  
  pred_df2 <- pred_df %>%
    mutate(
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      season = factor(.data$season, levels = c("Winter", "Summer")),
      x_base_comp = 1 + (as.numeric(.data$shelf_loc_200) - 1) * gap
    )
  
  list(
    df = df,
    df_mod = df_mod,
    model = m1,
    pred_df = pred_df,
    pred_df2 = pred_df2,
    string_ID = string_ID,
    PSS_var = PSS_var,
    depth_cut = depth_cut,
    log10_transform = log10_transform
  )
}

plot_compare_sources_row <- function(dat_all,
                                     dat_plankton,
                                     PSS_var = "slope",
                                     depth_cut = 100,
                                     log10_transform = FALSE,
                                     y_limits = NULL,
                                     season_colors = c("Winter" = "#0072B2",
                                                       "Summer" = "#E69F00"),
                                     source_colors = c("allParticles" = "black",
                                                       "plankton" = "gray"),
                                     source_shapes = c("allParticles" = 16,
                                                       "plankton" = 17),
                                     violin_dodge = 0.75,
                                     emm_dodge = 0.22,
                                     ci_cap_width = 0.03,
                                     gap = 0.35) {
  
  obj_all <- prep_mm_for_compare(
    dat = dat_all,
    string_ID = "ISIIS_full",
    PSS_var = PSS_var,
    depth_cut = depth_cut,
    log10_transform = log10_transform,
    gap = gap
  )
  
  obj_plk <- prep_mm_for_compare(
    dat = dat_plankton,
    string_ID = "ISIIS_plankton",
    PSS_var = PSS_var,
    depth_cut = depth_cut,
    log10_transform = log10_transform,
    gap = gap
  )
  
  df_both <- bind_rows(
    obj_all$df %>% mutate(data_source = "allParticles"),
    obj_plk$df %>% mutate(data_source = "plankton")
  ) %>%
    mutate(
      data_source = factor(data_source, levels = c("allParticles", "plankton")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  pred_both <- bind_rows(
    obj_all$pred_df2 %>% mutate(data_source = "allParticles"),
    obj_plk$pred_df2 %>% mutate(data_source = "plankton")
  ) %>%
    mutate(
      data_source = factor(data_source, levels = c("allParticles", "plankton")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      
      # two source-centered subgroups within each shelf location
      source_offset = case_when(
        data_source == "allParticles" ~ -0.09,
        data_source == "plankton"     ~  0.09,
        TRUE ~ 0
      ),
      
      # winter first, summer second within each source subgroup
      season_offset = case_when(
        season == "Winter" ~ -0.03,
        season == "Summer" ~  0.03,
        TRUE ~ 0
      ),
      
      x = x_base_comp + source_offset + season_offset
    )
  
  seg_both <- pred_both %>%
    dplyr::select(.data$shelf_loc_200, .data$data_source, .data$season, .data$x, .data$emmean) %>%
    pivot_wider(
      names_from = .data$season,
      values_from = c(.data$x, .data$emmean)
    ) %>%
    transmute(
      shelf_loc_200,
      data_source,
      x = .data$x_Winter,
      xend = .data$x_Summer,
      y = .data$emmean_Winter,
      yend = .data$emmean_Summer
    )
  
  if (PSS_var %in% c("total_biovolume", "total_biovolume_um3_m3", "integrated_biovol_all", "integrated_biovol_sample")) {
    y_label_expr <- expression(Biovolume~(mu*m^3~m^-3))
  } else if (PSS_var == "intercept") {
    y_label_expr <- expression(Intercept~(mu*m^3~m^-3~mu*m^-3))
  } else if (PSS_var == "slope") {
    y_label_expr <- expression(Slope~(m^-3~mu*m^-3))
  } else {
    y_label_expr <- PSS_var
  }
  
  y_scale_obj <- if (log10_transform || PSS_var == "intercept") {
    ggplot2::scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", round(x, 2)))
    )
  } else {
    NULL
  }
  
  if (is.null(y_limits)) {
    y_limits <- range(
      df_both$.y,
      pred_both$lower.CL,
      pred_both$upper.CL,
      na.rm = TRUE
    )
  }
  
  p_violin <- ggplot(
    df_both,
    aes(
      x = .data$shelf_loc_200,
      y = .data$.y,
      group = interaction(.data$shelf_loc_200, .data$season, .data$data_source)
    )
  ) +
    geom_violin(
      aes(fill = .data$season, color = .data$data_source),
      position = position_dodge(width = violin_dodge),
      trim = FALSE,
      alpha = 0.65,
      linewidth = 0.3
    ) +
    geom_boxplot(
      aes(
        fill = .data$season,
        color = .data$data_source,
        group = interaction(.data$shelf_loc_200, .data$season, .data$data_source)
      ),
      position = position_dodge(width = violin_dodge),
      width = 0.14,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    scale_x_discrete(
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      )
    ) +
    scale_fill_manual(values = season_colors, name = "Season") +
    scale_color_manual(
      values = source_colors,
      name = "Source",
      labels = c(
        allParticles = "Seston",
        plankton = "Plankton"
      )
    ) +
    labs(
      x = NULL,
      y = y_label_expr
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  p_emm <- ggplot(
    pred_both,
    aes(
      x = .data$x,
      y = .data$emmean,
      color = .data$season,
      shape = .data$data_source
    )
  ) +
    geom_segment(
      data = seg_both,
      aes(
        x = .data$x,
        xend = .data$xend,
        y = .data$y,
        yend = .data$yend,
        group = interaction(.data$shelf_loc_200, .data$data_source)
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.8
    ) +
    geom_errorbar(
      aes(ymin = .data$lower.CL, ymax = .data$upper.CL),
      width = ci_cap_width,
      linewidth = 0.6
    ) +
    geom_point(size = 3) +
    scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("off shelf", "on shelf")
    ) +
    scale_color_manual(values = season_colors, name = "Season") +
    scale_shape_manual(
      values = source_shapes,
      name = "Source",
      labels = c(
        allParticles = "Seston",
        plankton = "Plankton"
      )
    )+
    labs(
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  if (!is.null(y_scale_obj)) {
    p_violin <- p_violin + y_scale_obj
    p_emm    <- p_emm + y_scale_obj
  }
  
  p_violin + p_emm + patchwork::plot_layout(widths = c(1.3, 1))
}

make_allParticles_vs_plankton_3x2 <- function(dat_all,
                                              dat_plankton,
                                              output_dir = output_dir_figs_NB,
                                              width = 13,
                                              height = 16,
                                              dpi = 300) {
  
  row_slope <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "slope",
    depth_cut = 100,
    log10_transform = FALSE,
    y_limits = c(-2.5, 0)
  )
  
  row_intercept <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "intercept",
    depth_cut = 100,
    log10_transform = FALSE,
    y_limits = c(2, 20)
  )
  
  row_integrated <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "integrated_biovol_all",
    depth_cut = 100,
    log10_transform = TRUE,
    y_limits = c(7, 19)
  )
  
  p_final <- (row_slope / row_intercept / row_integrated) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(p_final)
  
  out_path <- file.path(
    output_dir,
    "SPECTRA_allParticles_vs_plankton_violin_plusMixedModel_3x2.png"
  )
  
  ggsave(
    out_path,
    p_final,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  message("Saved combined comparison figure to: ", out_path)
  
  invisible(p_final)
}

make_allParticles_vs_plankton_totalBiovol_3x2 <- function(dat_all,
                                                          dat_plankton,
                                                          output_dir = output_dir_figs_NB,
                                                          width = 13,
                                                          height = 16,
                                                          dpi = 300) {
  
  row_slope <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "slope",
    depth_cut = 100,
    log10_transform = FALSE,
    y_limits = c(-2.5, 0)
  )
  
  row_intercept <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "intercept",
    depth_cut = 100,
    log10_transform = FALSE,
    y_limits = c(2, 20)
  )
  
  row_total_biovol <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "total_biovolume_um3_m3",
    depth_cut = 100,
    log10_transform = TRUE,
    y_limits = c(7, 19)
  )
  
  p_final <- (row_slope / row_intercept / row_total_biovol) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(p_final)
  
  out_path <- file.path(
    output_dir,
    "SPECTRA_allParticles_vs_plankton_violin_plusMixedModel_totalBiovolume_3x2.png"
  )
  
  ggsave(
    out_path,
    p_final,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  message("Saved total-biovolume comparison figure to: ", out_path)
  
  invisible(p_final)
}

run_mm_set(dat_ISIIS_full_analysis, "ISIIS_full", analysis_output_dir)
run_mm_set(dat_ISIIS_plankton_analysis, "ISIIS_plankton", analysis_output_dir)

if (!opt$thinned) {
  run_mm_set(dat_ISIIS_detritus, "ISIIS_detritus", analysis_output_dir)
}

make_allParticles_vs_plankton_3x2(
  dat_all = dat_ISIIS_full_analysis,
  dat_plankton = dat_ISIIS_plankton_analysis,
  output_dir = analysis_output_dir
)

make_allParticles_vs_plankton_totalBiovol_3x2(
  dat_all = dat_ISIIS_full_analysis,
  dat_plankton = dat_ISIIS_plankton_analysis,
  output_dir = analysis_output_dir
)

# ============================================================
# 2x4 figure:
#   Y = log10(integrated_biovol_all)
#   X = slope
#   rows = Seston, Plankton
#   cols = w22, s22, w23, s23
#   separate linear fits for shelf_loc_200
# Also saves regression summary table.
# ============================================================
# ============================================================
# Exploratory environmental and biovolume mixed regressions
# One model per assemblage x cruise x shelf location
# ============================================================
fit_exploratory_mixed_models <- function(reg_dat, response, predictor, analysis_name) {
  groups <- reg_dat %>%
    dplyr::group_by(source, Cruise, shelf_loc_200) %>%
    dplyr::group_split(.keep = TRUE)

  purrr::map_dfr(groups, function(d) {
    d <- d %>%
      dplyr::filter(is.finite(.data[[response]]), is.finite(.data[[predictor]]),
                    !is.na(Transect_ID), !is.na(cast_id)) %>%
      dplyr::mutate(Transect_ID = droplevels(factor(Transect_ID)),
                    cast_id = factor(cast_id))
    source_i <- as.character(d$source[1])
    cruise_i <- as.character(d$Cruise[1])
    shelf_i <- as.character(d$shelf_loc_200[1])
    n_i <- nrow(d)
    n_transects <- dplyr::n_distinct(d$Transect_ID)
    n_casts <- dplyr::n_distinct(interaction(d$Transect_ID, d$cast_id, drop = TRUE))
    base <- tibble::tibble(
      analysis = analysis_name, response = response, predictor = predictor,
      source = source_i, Cruise = cruise_i, shelf_loc_200 = shelf_i,
      n = n_i, n_transects = n_transects, n_casts = n_casts)
    blank <- function(status) {
      dplyr::bind_cols(base, tibble::tibble(
        intercept = NA_real_, intercept_se = NA_real_, intercept_ci_low = NA_real_,
        intercept_ci_high = NA_real_, intercept_p_value = NA_real_,
        slope_coef = NA_real_, slope_se = NA_real_, ci_low = NA_real_,
        ci_high = NA_real_, p_value = NA_real_,
        var_intercept = NA_real_, var_slope = NA_real_, cov_intercept_slope = NA_real_,
        marginal_r2 = NA_real_, conditional_r2 = NA_real_, singular_fit = NA,
        convergence_messages = NA_character_, model_status = status,
        equation_label = status))
    }
    if (n_i < 6 || n_transects < 2 || n_casts < 3 ||
        dplyr::n_distinct(d[[predictor]]) < 3 ||
        n_i <= n_casts + 2) {
      return(blank("Insufficient replication"))
    }
    fml <- stats::as.formula(paste(response, "~", predictor,
                                    "+ (1 | Transect_ID/cast_id)"))
    fit <- tryCatch(lmerTest::lmer(fml, data = d, REML = TRUE),
                    error = function(e) e)
    if (inherits(fit, "error")) return(blank(paste("Fit failed:", conditionMessage(fit))))
    coef_table <- coef(summary(fit))
    if (!all(c("(Intercept)", predictor) %in% rownames(coef_table))) {
      return(blank("Fixed effect not estimable"))
    }
    V <- as.matrix(stats::vcov(fit))[c("(Intercept)", predictor),
                                     c("(Intercept)", predictor), drop = FALSE]
    ci <- tryCatch(suppressMessages(stats::confint(fit, parm = c("(Intercept)", predictor),
                                    method = "Wald")), error = function(e) NULL)
    b0 <- unname(lme4::fixef(fit)["(Intercept)"])
    b1 <- unname(lme4::fixef(fit)[predictor])
    se0 <- coef_table["(Intercept)", "Std. Error"]
    se1 <- coef_table[predictor, "Std. Error"]
    # Satterthwaite p-values supplied by lmerTest
    p0 <- coef_table["(Intercept)", "Pr(>|t|)"]
    p1 <- coef_table[predictor, "Pr(>|t|)"]
    if (is.null(ci)) {
      ci0 <- b0 + c(-1, 1) * 1.96 * se0
      ci1 <- b1 + c(-1, 1) * 1.96 * se1
    } else {
      ci0 <- ci["(Intercept)", ]
      ci1 <- ci[predictor, ]
    }
    # Nakagawa R2: optional dependency; does not prevent fitting
    r2m <- r2c <- NA_real_
    if (requireNamespace("performance", quietly = TRUE)) {
      r2 <- tryCatch(suppressWarnings(performance::r2_nakagawa(fit)),
                     error = function(e) NULL)
      if (!is.null(r2)) {
        r2m <- as.numeric(r2$R2_marginal)
        r2c <- as.numeric(r2$R2_conditional)
      }
    }
    conv <- fit@optinfo$conv$lme4$messages
    conv_text <- if (length(conv)) paste(conv, collapse = "; ") else "None"
    singular <- lme4::isSingular(fit)
    eq <- paste0(ifelse(shelf_i == "off_shelf", "Off shelf", "On shelf"),
                 ": beta = ", formatC(b1, digits = 3, format = "fg"),
                 " (95% CI ", formatC(ci1[1], digits = 3, format = "fg"),
                 ", ", formatC(ci1[2], digits = 3, format = "fg"), ")",
                 "\np = ", format.pval(p1, digits = 2, eps = 1e-4))
    dplyr::bind_cols(base, tibble::tibble(
      intercept = b0, intercept_se = se0,
      intercept_ci_low = unname(ci0[1]), intercept_ci_high = unname(ci0[2]),
      intercept_p_value = p0,
      slope_coef = b1, slope_se = se1,
      ci_low = unname(ci1[1]), ci_high = unname(ci1[2]), p_value = p1,
      var_intercept = unname(V[1, 1]), var_slope = unname(V[2, 2]),
      cov_intercept_slope = unname(V[1, 2]),
      marginal_r2 = r2m, conditional_r2 = r2c,
      singular_fit = singular, convergence_messages = conv_text,
      model_status = if (singular) "Singular fit" else if (length(conv)) "Convergence warning" else "OK",
      equation_label = eq))
  })
}

# Population-level mixed-model lines and fixed-effect 95% confidence bands.
# These are confidence intervals for the mean, not prediction intervals.
mixed_prediction_grid <- function(df_panel, fit_panel, xvar, n_points = 100) {
  purrr::map_dfr(seq_len(nrow(fit_panel)), function(i) {
    f <- fit_panel[i, ]
    if (!is.finite(f$slope_coef) || !is.finite(f$var_intercept) ||
        !is.finite(f$var_slope) || !is.finite(f$cov_intercept_slope)) {
      return(tibble::tibble())
    }
    xs <- df_panel %>%
      dplyr::filter(as.character(shelf_loc_200) == as.character(f$shelf_loc_200)) %>%
      dplyr::pull(.data[[xvar]])
    xs <- xs[is.finite(xs)]
    if (length(xs) < 2 || diff(range(xs)) == 0) return(tibble::tibble())
    x <- seq(min(xs), max(xs), length.out = n_points)
    y <- f$intercept + f$slope_coef * x
    se <- sqrt(pmax(0, f$var_intercept + 2*x*f$cov_intercept_slope + x*x*f$var_slope))
    tibble::tibble(x = x, y = y, ymin = y - 1.96*se, ymax = y + 1.96*se,
                   shelf_loc_200 = f$shelf_loc_200)
  })
}

fit_biovol_vs_slope_models <- function(reg_dat) {
  fit_exploratory_mixed_models(reg_dat, response = "log10_biovolume",
                               predictor = "slope", analysis_name = "Biovolume vs slope")
}

make_annotation_positions <- function(df_panel,
                                      fit_panel,
                                      x_limits = c(-2, NA),
                                      y_limits = c(NA, 20)) {
  
  x_data_min <- min(df_panel$slope, na.rm = TRUE)
  x_data_max <- max(df_panel$slope, na.rm = TRUE)
  y_data_min <- min(df_panel$log10_biovolume, na.rm = TRUE)
  y_data_max <- max(df_panel$log10_biovolume, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_max - 0.05 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_max - 0.22 * (y_max - y_min),
        TRUE ~ y_max
      )
    )
}

plot_biovol_vs_slope_panel <- function(df_panel,
                                       fit_tbl,
                                       source_i,
                                       cruise_i,
                                       biovol_var = "total_biovolume_um3_m3",
                                       x_limits = c(-2, NA),
                                       y_limits = c(NA, 20),
                                       point_size = 0.8,
                                       point_alpha = 0.18,
                                       show_x_title = FALSE,
                                       show_y_title = FALSE,
                                       show_legend = FALSE,
                                       strip_title = NULL,
                                       panel_letter = NULL) {
  
  y_axis_title <- dplyr::case_when(
    biovol_var == "integrated_biovol_all" ~ "Integrated biovolume",
    biovol_var == "total_biovolume_um3_m3" ~ "Total biovolume",
    TRUE ~ "Biovolume"
  )
  
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  pred_panel <- mixed_prediction_grid(df_panel, fit_panel, "slope")

  p <- ggplot(
    df_panel,
    aes(x = slope, y = log10_biovolume, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_ribbon(
      data = pred_panel,
      aes(x = x, ymin = ymin, ymax = ymax, fill = shelf_loc_200),
      inherit.aes = FALSE, alpha = 0.13, color = NA
    ) +
    geom_line(
      data = pred_panel,
      aes(x = x, y = y, color = shelf_loc_200),
      inherit.aes = FALSE, linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", x))
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(Slope~(m^-3~mu*m^-3)) else NULL,
      y = if (show_y_title) {
        bquote(.(y_axis_title)~(mu*m^3~m^-3))
      } else {
        NULL
      }
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

prep_biovol_vs_slope_data <- function(dat, source_name, depth_cut = 100, biovol_var = "total_biovolume_um3_m3") {
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
           season = factor(season, levels = c("Winter", "Summer")),
           shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))) %>%
    filter(!is.na(Cruise), !is.na(Transect_ID), !is.na(group_id), !is.na(shelf_loc_200)) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(slope = first(slope),
              biovolume = first(.data[[biovol_var]]),
              cast_id = first(cast_id),
              .groups = "drop") %>%
    mutate(source = source_name,
           log10_biovolume = if_else(is.finite(biovolume) & biovolume > 0, log10(biovolume), NA_real_)) %>%
    filter(is.finite(slope), is.finite(log10_biovolume))
}

plot_biovol_vs_slope_2x4 <- function(dat_all,
                                     dat_plankton,
                                     depth_cut = 100,
                                     biovol_var = "total_biovolume_um3_m3",
                                     output_name = "total_biovolume_um3_m3",
                                     output_dir = analysis_output_dir,
                                     width = 16,
                                     height = 8,
                                     dpi = 300,
                                     save_plot = TRUE,
                                     save_table = TRUE) {
  
  reg_all <- prep_biovol_vs_slope_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut,
    biovol_var = biovol_var
  )
  
  reg_plk <- prep_biovol_vs_slope_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut,
    biovol_var = biovol_var
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for biovolume-vs-slope plotting after filtering.")
  }
  
  fit_tbl <- fit_biovol_vs_slope_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_biovol_vs_slope_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      biovol_var = biovol_var,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(-2, NA)
  ylim_default <- c(7.5, 20)
  
  xlim_s23_plk <- c(NA, NA)
  ylim_s23_plk <- c(7.5, 20)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_s23_plk,
    y_limits = ylim_s23_plk,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_", output_name, "_vs_slope_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved biovolume-vs-slope figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      output_dir,
      "tables",
      paste0("SPECTRA_", output_name, "_vs_slope_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}

# Exploratory slope-biovolume analysis: observed total biovolume only.
# Integrated biovolume is derived from the fitted spectrum and is excluded here.
biovol_slope_out <- plot_biovol_vs_slope_2x4(
  dat_all = dat_ISIIS_full_analysis,
  dat_plankton = dat_ISIIS_plankton_analysis,
  depth_cut = 100,
  biovol_var = "total_biovolume_um3_m3",
  output_name = "total_biovolume_um3_m3",
  output_dir = analysis_output_dir,
  save_plot = TRUE,
  save_table = TRUE
)

#Table S6:
suppressPackageStartupMessages({
  library(gt)
})

# ------------------------------------------------
# Collapse to group-level observations
# ------------------------------------------------

prep_group_level <- function(dat, source_name, depth_cut = 100) {
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      intercept_untrans = 10^intercept) %>%
    filter(
      !is.na(season),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      intercept = first(intercept_untrans),
      integrated_biovolume = first(integrated_biovol_all),
      .groups = "drop"
    ) %>%
    mutate(source = source_name)
  
}

group_all <- prep_group_level(dat_ISIIS_full_analysis, "Seston")
group_plk <- prep_group_level(dat_ISIIS_plankton_analysis, "Plankton")

group_all_data <- bind_rows(group_all, group_plk)

# ------------------------------------------------
# Summary statistics
# ------------------------------------------------

summary_stats <- group_all_data %>%
  group_by(source, season, shelf_loc_200) %>%
  summarise(
    
    slope_med = median(slope, na.rm = TRUE),
    slope_sd  = sd(slope, na.rm = TRUE),
    
    int_med = median(intercept, na.rm = TRUE),
    int_q1  = quantile(intercept, 0.25, na.rm = TRUE, names = FALSE),
    int_q3  = quantile(intercept, 0.75, na.rm = TRUE, names = FALSE),
    
    bio_med = median(integrated_biovolume, na.rm = TRUE),
    bio_q1  = quantile(integrated_biovolume, 0.25, na.rm = TRUE, names = FALSE),
    bio_q3  = quantile(integrated_biovolume, 0.75, na.rm = TRUE, names = FALSE),
    
    .groups = "drop"
  )

# ------------------------------------------------
# Formatting helpers
# ------------------------------------------------

fmt_slope <- function(med, sd) {
  if (!is.finite(med)) return(NA_character_)
  if (!is.finite(sd)) {
    return(formatC(med, format = "f", digits = 2))
  }
  paste0(
    formatC(med, format = "f", digits = 2),
    " ± ",
    formatC(sd, format = "f", digits = 2)
  )
}

fmt_range <- function(med, q1, q3) {
  if (!is.finite(med)) return(NA_character_)
  if (!is.finite(q1) || !is.finite(q3)) {
    return(formatC(med, format = "e", digits = 2))
  }
  paste0(
    formatC(med, format = "e", digits = 2),
    " [",
    formatC(q1, format = "e", digits = 2),
    "–",
    formatC(q3, format = "e", digits = 2),
    "]"
  )
}


# ------------------------------------------------
# Convert to long metric rows
# ------------------------------------------------

table_long <- summary_stats %>%
  mutate(
    slope = mapply(fmt_slope, slope_med, slope_sd),
    intercept = mapply(fmt_range, int_med, int_q1, int_q3),
    integrated_biovolume = mapply(fmt_range, bio_med, bio_q1, bio_q3)
  ) %>%
  dplyr::select(source, season, shelf_loc_200, slope, intercept, integrated_biovolume) %>%
  pivot_longer(
    cols = c(slope, intercept, integrated_biovolume),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    source = factor(source, levels = c("Seston", "Plankton")),
    metric = factor(metric, levels = c("slope", "intercept", "integrated_biovolume")),
    row_name = paste(
      recode(metric,
             slope = "Slope",
             intercept = "Intercept",
             integrated_biovolume = "Integrated biovolume"
      ),
      source,
      sep = " – "
    ),
    column_name = paste(
      as.character(season),
      ifelse(shelf_loc_200 == "off_shelf", "off shelf", "on shelf"),
      sep = " – "
    )
  )

# ------------------------------------------------
# Final wide table, ordered by metric first
# ------------------------------------------------

table_wide <- table_long %>%
  mutate(
    metric_order = factor(metric, levels = c("slope", "intercept", "integrated_biovolume")),
    source_order = factor(source, levels = c("Seston", "Plankton"))
  ) %>%
  arrange(metric_order, source_order) %>%
  dplyr::select(row_name, column_name, value, metric_order, source_order) %>%
  pivot_wider(
    names_from = column_name,
    values_from = value
  ) %>%
  arrange(metric_order, source_order) %>%
  dplyr::select(
    row_name,
    `Winter – off shelf`,
    `Winter – on shelf`,
    `Summer – off shelf`,
    `Summer – on shelf`
  )

# ------------------------------------------------
# Save CSV
# ------------------------------------------------
dir.create(file.path(analysis_output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  table_wide,
  file.path(analysis_output_dir, "tables", "season_shelf_PSS_coefficients_table.csv"),
  na = ""
)

# ------------------------------------------------
# GT formatted table
# ------------------------------------------------

gt_table <- table_wide %>%
  gt(rowname_col = "row_name") %>%
  cols_align("center", everything()) %>%
  tab_options(
    table.font.size = px(13),
    column_labels.font.weight = "bold"
  )

gtsave(
  gt_table,
  file.path(analysis_output_dir, "tables", "season_shelf_PSS_coefficients_table.html")
)

gt_table


# ------------------------------------------------
# Chla and slope
# ------------------------------------------------
prep_slope_vs_chla_data <- function(dat,
                                    source_name,
                                    depth_cut = 100) {
  
  needed <- c("slope", "chl_a_ug_l", "Cruise", "Transect_ID", "group_id", "cast_id", "season", "shelf_loc_200")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns for slope vs chl_a_ug_l plot: ",
         paste(miss, collapse = ", "))
  }
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      slope = as.numeric(slope),
      chl_a_ug_l = as.numeric(chl_a_ug_l)
    ) %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      chl_a_ug_l = first(chl_a_ug_l),
      cast_id = first(cast_id),
      .groups = "drop"
    ) %>%
    mutate(
      source = source_name,
      log10_chl_a_ug_l = if_else(
        is.finite(chl_a_ug_l) & chl_a_ug_l > 0,
        log10(chl_a_ug_l),
        NA_real_
      )
    ) %>%
    filter(
      is.finite(slope),
      is.finite(log10_chl_a_ug_l)
    )
}

fit_slope_vs_chla_models <- function(reg_dat) {
  fit_exploratory_mixed_models(reg_dat, response = "slope",
                               predictor = "log10_chl_a_ug_l", analysis_name = "Slope vs chlorophyll-a")
}

make_annotation_positions_slope_chla <- function(df_panel,
                                                 fit_panel,
                                                 x_limits = c(NA, NA),
                                                 y_limits = c(-2.5, 0)) {
  
  x_data_min <- min(df_panel$log10_chl_a_ug_l, na.rm = TRUE)
  x_data_max <- max(df_panel$log10_chl_a_ug_l, na.rm = TRUE)
  y_data_min <- min(df_panel$slope, na.rm = TRUE)
  y_data_max <- max(df_panel$slope, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_min + 0.28 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_min + 0.14 * (y_max - y_min),
        TRUE ~ y_min
      )
    )
}

plot_slope_vs_chla_panel <- function(df_panel,
                                     fit_tbl,
                                     source_i,
                                     cruise_i,
                                     x_limits = c(NA, NA),
                                     y_limits = c(-2.5, 0),
                                     point_size = 0.8,
                                     point_alpha = 0.18,
                                     show_x_title = FALSE,
                                     show_y_title = FALSE,
                                     show_legend = FALSE,
                                     strip_title = NULL,
                                     panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions_slope_chla(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  pred_panel <- mixed_prediction_grid(df_panel, fit_panel, "log10_chl_a_ug_l")

  p <- ggplot(
    df_panel,
    aes(x = log10_chl_a_ug_l, y = slope, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_ribbon(
      data = pred_panel,
      aes(x = x, ymin = ymin, ymax = ymax, fill = shelf_loc_200),
      inherit.aes = FALSE, alpha = 0.13, color = NA
    ) +
    geom_line(
      data = pred_panel,
      aes(x = x, y = y, color = shelf_loc_200),
      inherit.aes = FALSE, linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(log[10]~chl[a]~(mu*g~L^-1)) else NULL,
      y = if (show_y_title) expression(Slope~(m^-3~mu*m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_slope_vs_chla_2x4 <- function(dat_all,
                                   dat_plankton,
                                   depth_cut = 100,
                                   output_dir = output_dir_figs_NB,
                                   width = 16,
                                   height = 8,
                                   dpi = 300,
                                   save_plot = TRUE,
                                   save_table = TRUE) {
  
  reg_all <- prep_slope_vs_chla_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_slope_vs_chla_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for slope-vs-chla plotting after filtering.")
  }
  
  fit_tbl <- fit_slope_vs_chla_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_slope_vs_chla_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(NA, NA)
  ylim_default <- c(-2.5, 0)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_slope_vs_log10_chla_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved slope-vs-chla figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      output_dir,
      "tables",
      paste0("SPECTRA_slope_vs_log10_chla_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}

slope_chla_out <- plot_slope_vs_chla_2x4(
  dat_all = dat_ISIIS_full_analysis,
  dat_plankton = dat_ISIIS_plankton_analysis,
  depth_cut = 100,
  output_dir = analysis_output_dir,
  save_plot = TRUE,
  save_table = TRUE
)

# ------------------------------------------------
# temperature and slope
# ------------------------------------------------

prep_slope_vs_temperature_data <- function(dat,
                                           source_name,
                                           depth_cut = 100) {
  
  needed <- c("slope", "Temperature", "Cruise", "Transect_ID", "group_id", "cast_id", "season", "shelf_loc_200")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns for slope vs Temperature plot: ",
         paste(miss, collapse = ", "))
  }
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      slope = as.numeric(slope),
      Temperature = as.numeric(Temperature)
    ) %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      Temperature = first(Temperature),
      cast_id = first(cast_id),
      .groups = "drop"
    ) %>%
    mutate(source = source_name) %>%
    filter(
      is.finite(slope),
      is.finite(Temperature)
    )
}

fit_slope_vs_temperature_models <- function(reg_dat) {
  fit_exploratory_mixed_models(reg_dat, response = "slope",
                               predictor = "Temperature", analysis_name = "Slope vs temperature")
}

make_annotation_positions_slope_temperature <- function(df_panel,
                                                        fit_panel,
                                                        x_limits = c(NA, NA),
                                                        y_limits = c(-2.5, 0)) {
  
  x_data_min <- min(df_panel$Temperature, na.rm = TRUE)
  x_data_max <- max(df_panel$Temperature, na.rm = TRUE)
  y_data_min <- min(df_panel$slope, na.rm = TRUE)
  y_data_max <- max(df_panel$slope, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_min + 0.28 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_min + 0.14 * (y_max - y_min),
        TRUE ~ y_min
      )
    )
}

plot_slope_vs_temperature_panel <- function(df_panel,
                                            fit_tbl,
                                            source_i,
                                            cruise_i,
                                            x_limits = c(NA, NA),
                                            y_limits = c(-2.5, 0),
                                            point_size = 0.8,
                                            point_alpha = 0.18,
                                            show_x_title = FALSE,
                                            show_y_title = FALSE,
                                            show_legend = FALSE,
                                            strip_title = NULL,
                                            panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions_slope_temperature(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  pred_panel <- mixed_prediction_grid(df_panel, fit_panel, "Temperature")

  p <- ggplot(
    df_panel,
    aes(x = Temperature, y = slope, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_ribbon(
      data = pred_panel,
      aes(x = x, ymin = ymin, ymax = ymax, fill = shelf_loc_200),
      inherit.aes = FALSE, alpha = 0.13, color = NA
    ) +
    geom_line(
      data = pred_panel,
      aes(x = x, y = y, color = shelf_loc_200),
      inherit.aes = FALSE, linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(Temperature~(degree*C)) else NULL,
      y = if (show_y_title) expression(Slope~(m^-3~mu*m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_slope_vs_temperature_2x4 <- function(dat_all,
                                          dat_plankton,
                                          depth_cut = 100,
                                          output_dir = output_dir_figs_NB,
                                          width = 16,
                                          height = 8,
                                          dpi = 300,
                                          save_plot = TRUE,
                                          save_table = TRUE) {
  
  reg_all <- prep_slope_vs_temperature_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_slope_vs_temperature_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for slope-vs-temperature plotting after filtering.")
  }
  
  fit_tbl <- fit_slope_vs_temperature_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_slope_vs_temperature_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(NA, NA)
  ylim_default <- c(-2.5, 0)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_slope_vs_temperature_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved slope-vs-temperature figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      output_dir,
      "tables",
      paste0("SPECTRA_slope_vs_temperature_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}

slope_temperature_out <- plot_slope_vs_temperature_2x4(
  dat_all = dat_ISIIS_full_analysis,
  dat_plankton = dat_ISIIS_plankton_analysis,
  depth_cut = 100,
  output_dir = analysis_output_dir,
  save_plot = TRUE,
  save_table = TRUE
)

# ------------------------------------------------
# oxygen and slope
# ------------------------------------------------

prep_slope_vs_oxygen_data <- function(dat,
                                           source_name,
                                           depth_cut = 100) {
  
  needed <- c("slope", "Oxygen_ml_l", "Cruise", "Transect_ID", "group_id", "cast_id", "season", "shelf_loc_200")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns for slope vs oxygen plot: ",
         paste(miss, collapse = ", "))
  }
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      slope = as.numeric(slope),
      Oxygen_ml_l = as.numeric(Oxygen_ml_l)
    ) %>%
    dplyr::filter(.data$Cruise != "s23") %>%   # skip incorrect Summer 2023 oxygen %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      Oxygen_ml_l = first(Oxygen_ml_l),
      cast_id = first(cast_id),
      .groups = "drop"
    ) %>%
    mutate(source = source_name) %>%
    filter(
      is.finite(slope),
      is.finite(Oxygen_ml_l)
    )
}

fit_slope_vs_oxygen_models <- function(reg_dat) {
  fit_exploratory_mixed_models(reg_dat, response = "slope",
                               predictor = "Oxygen_ml_l", analysis_name = "Slope vs oxygen")
}

make_annotation_positions_slope_oxygen <- function(df_panel,
                                                        fit_panel,
                                                        x_limits = c(NA, NA),
                                                        y_limits = c(-2.5, 0)) {
  
  x_data_min <- min(df_panel$Oxygen_ml_l, na.rm = TRUE)
  x_data_max <- max(df_panel$Oxygen_ml_l, na.rm = TRUE)
  y_data_min <- min(df_panel$slope, na.rm = TRUE)
  y_data_max <- max(df_panel$slope, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_min + 0.28 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_min + 0.14 * (y_max - y_min),
        TRUE ~ y_min
      )
    )
}

plot_slope_vs_oxygen_panel <- function(df_panel,
                                            fit_tbl,
                                            source_i,
                                            cruise_i,
                                            x_limits = c(NA, NA),
                                            y_limits = c(-2.5, 0),
                                            point_size = 0.8,
                                            point_alpha = 0.18,
                                            show_x_title = FALSE,
                                            show_y_title = FALSE,
                                            show_legend = FALSE,
                                            strip_title = NULL,
                                            panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions_slope_oxygen(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  pred_panel <- mixed_prediction_grid(df_panel, fit_panel, "Oxygen_ml_l")

  p <- ggplot(
    df_panel,
    aes(x = Oxygen_ml_l, y = slope, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_ribbon(
      data = pred_panel,
      aes(x = x, ymin = ymin, ymax = ymax, fill = shelf_loc_200),
      inherit.aes = FALSE, alpha = 0.13, color = NA
    ) +
    geom_line(
      data = pred_panel,
      aes(x = x, y = y, color = shelf_loc_200),
      inherit.aes = FALSE, linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(Oxygen~(mL~L^-1)) else NULL,
      y = if (show_y_title) expression(Slope~(m^-3~mu*m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_slope_vs_oxygen_2x4 <- function(dat_all,
                                          dat_plankton,
                                          depth_cut = 100,
                                          output_dir = output_dir_figs_NB,
                                          width = 16,
                                          height = 8,
                                          dpi = 300,
                                          save_plot = TRUE,
                                          save_table = TRUE) {
  
  reg_all <- prep_slope_vs_oxygen_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_slope_vs_oxygen_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for slope-vs-oxygen plotting after filtering.")
  }
  
  fit_tbl <- fit_slope_vs_oxygen_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_slope_vs_oxygen_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(NA, NA)
  ylim_default <- c(-2.5, 0)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_slope_vs_oxygen_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved slope-vs-oxygen figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      output_dir,
      "tables",
      paste0("SPECTRA_slope_vs_oxygen_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}


slope_oxygen_out <- plot_slope_vs_oxygen_2x4(
  dat_all = dat_ISIIS_full_analysis,
  dat_plankton = dat_ISIIS_plankton_analysis,
  depth_cut = 100,
  output_dir = analysis_output_dir,
  save_plot = TRUE,
  save_table = TRUE
)



# ============================================================
# Combined fixed-effects table: one row per coefficient per model
# ============================================================
dir.create(file.path(analysis_output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

exploratory_model_results <- dplyr::bind_rows(
  biovol_slope_out$regression_table,
  slope_chla_out$regression_table,
  slope_temperature_out$regression_table,
  slope_oxygen_out$regression_table
)

exploratory_fixed_effects <- dplyr::bind_rows(
  exploratory_model_results %>%
    dplyr::transmute(analysis, response, predictor, source, Cruise, shelf_loc_200,
                     n, n_transects, n_casts, term = "(Intercept)",
                     estimate = intercept, std_error = intercept_se,
                     conf_low = intercept_ci_low, conf_high = intercept_ci_high,
                     p_value = intercept_p_value, marginal_r2, conditional_r2,
                     singular_fit, convergence_messages, model_status),
  exploratory_model_results %>%
    dplyr::transmute(analysis, response, predictor, source, Cruise, shelf_loc_200,
                     n, n_transects, n_casts, term = predictor,
                     estimate = slope_coef, std_error = slope_se,
                     conf_low = ci_low, conf_high = ci_high,
                     p_value, marginal_r2, conditional_r2,
                     singular_fit, convergence_messages, model_status)
) %>%
  dplyr::arrange(analysis, source, Cruise, shelf_loc_200, term)

readr::write_csv(
  exploratory_fixed_effects,
  file.path(analysis_output_dir, "tables", "SPECTRA_exploratory_mixed_models_fixed_effects.csv")
)
readr::write_csv(
  exploratory_model_results,
  file.path(analysis_output_dir, "tables", "SPECTRA_exploratory_mixed_models_by_panel.csv")
)
