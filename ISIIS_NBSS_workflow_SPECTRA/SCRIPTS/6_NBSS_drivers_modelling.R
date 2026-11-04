setwd('/Users/user/SCRIPTS_backups/')
library(tidyverse)
library(gstat)
library(sp)
library(akima)
library(metR)
library(ggthemes)
library(fields)
library(sf)
library(colorRamps)
library(caret)
library(ranger)
library(iml)
library(MBA)
library(scales)

root_dir <- getwd()
source(file.path(root_dir, 'ISIIS_NBSS_workflow_SPECTRA', 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))
source(file.path(root_dir, 'ISIIS_NBSS_workflow_SPECTRA', 'SCRIPTS', 'functions', 'model_functions.R'))

#projects dirs
project_dirs <- c(
  MEZCAL        = file.path(root_dir, "ISIIS_NBSS_workflow_MEZCAL"),
  `MBON_ISIIS-2` = file.path(root_dir, "ISIIS_NBSS_workflow_MBON_ISIIS-2"),
  `MBON_ISIIS-3` = file.path(root_dir, "ISIIS_NBSS_workflow_MBON_ISIIS-3"),
  SPECTRA       = file.path(root_dir, "ISIIS_NBSS_workflow_SPECTRA")
)

# create paths and file lists

env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
covariates_list<- c('depth_mid', 'time_of_day')
random_effects_list<- c('Cruise', 'Transect_ID')


tow_id_candidates <- c(
  "tow_id",
  "tow_ID",
  "Tow_ID",
  "TOW_ID"
)

#functions
# Standardize the tow-ID column name when one is present
standardize_tow_id <- function(df) {
  
  tow_col <- intersect(tow_id_candidates, names(df))
  
  if (length(tow_col) > 1L) {
    stop(
      "More than one possible tow-ID column was found: ",
      paste(tow_col, collapse = ", ")
    )
  }
  
  if (length(tow_col) == 1L && tow_col != "tow_id") {
    df <- df %>%
      dplyr::rename(tow_id = dplyr::all_of(tow_col))
  }
  
  df
}


load_project_NBSS <- function(project_dir, project_name, env_vars) {
  
  size_spectra_dir <- file.path(project_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded")
  
  df_list_NBSS <- list.files(size_spectra_dir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
  
  message(project_name, ": found ", length(df_list_NBSS), " NBSS files")
  
  if (length(df_list_NBSS) == 0L) {
    warning("No NBSS files found for ", project_name, " in ", size_spectra_dir)
    return(tibble())
  }
  
  dat <- purrr::map_dfr(df_list_NBSS, ~ readr::read_csv( .x, show_col_types = FALSE,
      col_select = c(
        Cruise,
        Transect_ID,
        depth_mid,
        Latitude,
        Longitude,
        group_id,
        total_biovolume,
        biovol_um3_mid,
        biovol_um3_width,
        time_Pacific,
        time_of_day,
        rel_dist,
        distance_sample_km,
        slope,
        intercept,
        dplyr::any_of(c("tow_id","tow_ID","Tow_ID","R2", "residual_logNB",env_vars))))) %>%
    dplyr::arrange(time_Pacific) %>%
    dplyr::rename(group_id_old = group_id) %>%
    dplyr::mutate(project = project_name,
                  group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
      dplyr::across(dplyr::any_of(c("total_biovolume", "biovol_um3_mid", "biovol_um3_width", "slope", "intercept","R2", "residual_logNB")), as.numeric)
    ) 
  
  dat <- standardize_tow_id(dat)

  dat <- dat %>% dplyr::group_by(group_id) %>%
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
      total_biovolume_sum = sum(total_biovolume,na.rm = TRUE),
      R2 = dplyr::first(R2),
      residual_logNB = dplyr::first(residual_logNB),
      dplyr::across(dplyr::any_of(env_vars), dplyr::first),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      total_biovolume_sum = as.numeric(total_biovolume_sum)
    ) %>%
    dplyr::filter(
      !is.na(slope),
      slope < 0
    ) %>%
    dplyr::arrange(
      Latitude,
      Longitude,
      time_Pacific
    ) %>%
    dplyr::distinct(
      group_id,
      .keep_all = TRUE
    )

  dat
}

load_project_abundance <- function(project_dir, project_name) {
  
  abundance_dir <- file.path(
    project_dir,
    "data",
    "binned_abundance"
  )
  
  df_list_abundance <- list.files(
    abundance_dir,
    pattern = "chunk_1\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  df_list_abundance <- df_list_abundance[
    !grepl("test", df_list_abundance, ignore.case = TRUE)
  ]
  
  
  
  message(
    project_name,
    ": found ",
    length(df_list_abundance),
    " abundance files"
  )
  
  if (length(df_list_abundance) == 0L) {
    warning(
      "No abundance files found for ",
      project_name,
      " in ",
      abundance_dir
    )
    
    return(tibble())
  }
  
  dat_abundance <- purrr::map_dfr(
    df_list_abundance,
    ~ readr::read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        Cruise_name,
        Transect,
        cast_id,
        depth_run_id,
        depth_mid,
        mean_time_Pacific,
        mean_Latitude,
        mean_Longitude,
        dplyr::any_of(c(
          "tow_id",
          "tow_ID",
          "Tow_ID"
        )),
        dplyr::matches("^ind_m3")
      )
    )
  )
  
  
 dat_abundance <- standardize_tow_id(dat_abundance)
  
 dat_abundance<- dat_abundance %>% dplyr::rename(
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
  
  # Average values within 5-m depth bins
  dat_abundance <- dat_abundance %>%
    depth_bin_calc(depth_bin_width = 5) %>%
    dplyr::relocate(depth_mid, .before = 5)
  
  abundance_group_cols <- c(
    "Cruise",
    "Transect_ID",
    "cast_id",
    if ("tow_id" %in% names(dat_abundance)) "tow_id",
    "depth_mid"
  )
  
  dat_abundance_grouped <- dat_abundance %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(abundance_group_cols))
    ) %>%
    dplyr::summarise(
      project = dplyr::first(project),
      
      time_Pacific = {
        valid_times <- sort(time_Pacific[!is.na(time_Pacific)])
        
        if (length(valid_times) == 0L) {
          as.POSIXct(NA)
        } else {
          valid_times[ceiling(length(valid_times) / 2)]
        }
      },
      
      Latitude = median(Latitude, na.rm = TRUE),
      Longitude = median(Longitude, na.rm = TRUE),
      
      dplyr::across(
        dplyr::starts_with("ind_m3_"),
        ~ mean(.x, na.rm = TRUE)
      ),
      
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      Cruise,
      Transect_ID,
      cast_id,
      depth_mid
    )
  
  return(dat_abundance_grouped)
}

# create dat lists
dat_NBSS_list <- purrr::imap(
  project_dirs,
  ~ load_project_NBSS(
    project_dir = .x,
    project_name = .y,
    env_vars = env_vars
  )
)

dat_abundance_list <- purrr::imap(
  project_dirs,
  ~ load_project_abundance(
    project_dir = .x,
    project_name = .y
  )
)

#Table to check discrepancies in time stamps:

summarize_time_ranges <- function(df, project_name, dataset_name) {
  
  df <- standardize_tow_id(df)
  
  required_cols <- c(
    "Cruise",
    "Transect_ID",
    "time_Pacific"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      project_name,
      " ",
      dataset_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Use tow_id only when it exists
  grouping_cols <- c(
    "Cruise",
    "Transect_ID",
    if ("tow_id" %in% names(df)) "tow_id"
  )
  
  df %>%
    dplyr::filter(!is.na(time_Pacific)) %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(grouping_cols))
    ) %>%
    dplyr::summarise(
      min_time = min(time_Pacific, na.rm = TRUE),
      max_time = max(time_Pacific, na.rm = TRUE),
      n_rows = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      project = project_name,
      dataset = dataset_name,
      .before = 1
    )
}

compare_project_time_ranges <- function(
    dat_NBSS,
    dat_abundance,
    project_name
) {
  
  dat_NBSS <- standardize_tow_id(dat_NBSS)
  dat_abundance <- standardize_tow_id(dat_abundance)
  
  tow_required <- project_name %in% c(
    "MEZCAL",
    "MBON_ISIIS-2"
  )
  
  if (tow_required && !"tow_id" %in% names(dat_NBSS)) {
    stop(
      "tow_id is required but missing from NBSS data for ",
      project_name
    )
  }
  
  if (tow_required && !"tow_id" %in% names(dat_abundance)) {
    stop(
      "tow_id is required but missing from abundance data for ",
      project_name
    )
  }
  
  grouping_cols <- c(
    "Cruise",
    "Transect_ID",
    if (tow_required) "tow_id"
  )
  
  nbss_summary <- dat_NBSS %>%
    dplyr::filter(!is.na(time_Pacific)) %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(grouping_cols))
    ) %>%
    dplyr::summarise(
      min_time_NBSS = min(time_Pacific, na.rm = TRUE),
      max_time_NBSS = max(time_Pacific, na.rm = TRUE),
      n_rows_NBSS = dplyr::n(),
      .groups = "drop"
    )
  
  abundance_summary <- dat_abundance %>%
    dplyr::filter(!is.na(time_Pacific)) %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(grouping_cols))
    ) %>%
    dplyr::summarise(
      min_time_abundance = min(time_Pacific, na.rm = TRUE),
      max_time_abundance = max(time_Pacific, na.rm = TRUE),
      n_rows_abundance = dplyr::n(),
      .groups = "drop"
    )
  
  dplyr::full_join(
    nbss_summary,
    abundance_summary,
    by = grouping_cols
  ) %>%
    dplyr::mutate(
      project = project_name,
      .before = 1,
      
      min_time_difference_min = as.numeric(
        difftime(
          min_time_NBSS,
          min_time_abundance,
          units = "mins"
        )
      ),
      
      max_time_difference_min = as.numeric(
        difftime(
          max_time_NBSS,
          max_time_abundance,
          units = "mins"
        )
      ),
      
      datasets_present = dplyr::case_when(
        !is.na(min_time_NBSS) &
          !is.na(min_time_abundance) ~ "both",
        
        !is.na(min_time_NBSS) ~ "NBSS only",
        
        !is.na(min_time_abundance) ~ "abundance only",
        
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::arrange(
      Cruise,
      Transect_ID,
      dplyr::across(dplyr::any_of("tow_id"))
    )
}

project_names <- intersect(
  names(dat_NBSS_list),
  names(dat_abundance_list)
)

time_comparison_list <- purrr::map(
  project_names,
  \(project_name) {
    compare_project_time_ranges(
      dat_NBSS = dat_NBSS_list[[project_name]],
      dat_abundance = dat_abundance_list[[project_name]],
      project_name = project_name
    )
  }
) %>%
  rlang::set_names(project_names)

#function to plot depth_mid as function of time_Pacific
plot_abundance_depth_time <- function(dat, cruise, transect, tow = NULL, point_alpha = 0.7,point_size = 1.5) {
  
  required_cols <- c( "Cruise", "Transect_ID", "time_Pacific", "depth_mid")
  
  missing_cols <- setdiff(required_cols, names(dat))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Dataset is missing required columns: ",paste(missing_cols, collapse = ", "))
  }
  
  # Filter by Cruise and Transect_ID
  dat_plot <- dat %>% dplyr::filter( Cruise == cruise, Transect_ID == transect)
  
  # Filter by tow_id when requested
  if (!is.null(tow)) {
    
    if (!"tow_id" %in% names(dat_plot)) {
      stop("A tow was supplied, but this dataset does not contain tow_id.")
    }
    
    dat_plot <- dat_plot %>% dplyr::filter(tow_id == tow)
  }
  
  if (nrow(dat_plot) == 0L) {
    stop("No rows found for Cruise = ", cruise,", Transect_ID = ", transect,
      if (!is.null(tow)) paste0(", tow_id = ", tow) else "")
  }
  
  # Ensure correct classes
  dat_plot <- dat_plot %>% dplyr::mutate(
      time_Pacific = as.POSIXct(time_Pacific),
      depth_mid = as.numeric(depth_mid)
    ) %>%
    dplyr::filter( !is.na(time_Pacific), !is.na(depth_mid))
  
  title_text <- paste0(cruise, " – ", transect, if (!is.null(tow)) paste0(" – tow ", tow) else "")
  
  ggplot2::ggplot(dat_plot, ggplot2::aes(x = time_Pacific, y = depth_mid)) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d\n%H:%M", date_breaks = "2 hours") +
    ggplot2::labs(title = title_text, x = "Time Pacific", y = "Depth midpoint (m)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) 
}

plot_abundance_depth_time(
  dat = dat_abundance_list[["MEZCAL"]],
  cruise = "s19",
  transect = "NH",
  tow = 425
)

plot_abundance_depth_time(
  dat = dat_abundance_list[["MEZCAL"]],
  cruise = "w19",
  transect = "TR",
  tow = 325
)

plot_abundance_depth_time(
  dat = dat_abundance_list[["MBON_ISIIS-2"]],
  cruise = "f20",
  transect = "LP",
  tow = 121
)

plot_abundance_depth_time(
     dat = dat_abundance_list[["SPECTRA"]],
     cruise = "s23",
     transect = "GH"
  )

plot_abundance_depth_time(
  dat = dat_abundance_list[["MBON_ISIIS-3"]],
  cruise = "f24",
  transect = "HH",
  tow = "D"
  )

# function to remove data outliers based on time:
flag_time_outliers <- function(dat, gap_multiplier = 20, minimum_gap_minutes = 30, maximum_outlier_block_size = 2) {
  
  required_cols <- c("Cruise", "Transect_ID", "time_Pacific")
  missing_cols <- setdiff(required_cols, names(dat))
  
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
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
    dplyr::select(
      dplyr::all_of(group_cols),
      time_Pacific,
      gap_minutes,
      gap_threshold_minutes,
      time_block,
      n_times_in_block,
      time_outlier
    )
  
  dat %>%
    dplyr::left_join(time_blocks, by = c(group_cols, "time_Pacific"))
}

#a
dat_abundance_list[["MEZCAL"]] <- flag_time_outliers(
  dat = dat_abundance_list[["MEZCAL"]],
  gap_multiplier = 20,
  minimum_gap_minutes = 30,
  maximum_outlier_block_size = 2
)

dat_abundance_list[["MBON_ISIIS-2"]] <- flag_time_outliers(
  dat = dat_abundance_list[["MBON_ISIIS-2"]],
  gap_multiplier = 20,
  minimum_gap_minutes = 30,
  maximum_outlier_block_size = 2
)

dat_abundance_list[["SPECTRA"]] <- flag_time_outliers(
  dat = dat_abundance_list[["SPECTRA"]],
  gap_multiplier = 20,
  minimum_gap_minutes = 30,
  maximum_outlier_block_size = 2
)

dat_abundance_list[["MBON_ISIIS-3"]] <- flag_time_outliers(
  dat = dat_abundance_list[["MBON_ISIIS-3"]],
  gap_multiplier = 20,
  minimum_gap_minutes = 30,
  maximum_outlier_block_size = 2
)

#remove outliers
dat_abundance_list[["MEZCAL"]] <- dat_abundance_list[["MEZCAL"]] %>%
  dplyr::filter(!time_outlier) %>%
  dplyr::select(-gap_minutes, -gap_threshold_minutes, -time_block, -n_times_in_block, -time_outlier)

dat_abundance_list[["MBON_ISIIS-2"]] <- dat_abundance_list[["MBON_ISIIS-2"]] %>%
  dplyr::filter(!time_outlier) %>%
  dplyr::select(-gap_minutes, -gap_threshold_minutes, -time_block, -n_times_in_block, -time_outlier)

dat_abundance_list[["SPECTRA"]] <- dat_abundance_list[["SPECTRA"]] %>%
  dplyr::filter(!time_outlier) %>%
  dplyr::select(-gap_minutes, -gap_threshold_minutes, -time_block, -n_times_in_block, -time_outlier)

dat_abundance_list[["MBON_ISIIS-3"]] <- dat_abundance_list[["MBON_ISIIS-3"]] %>%
  dplyr::filter(!time_outlier) %>%
  dplyr::select(-gap_minutes, -gap_threshold_minutes, -time_block, -n_times_in_block, -time_outlier)


# remove problematic transects




#load and arrange LISST data
dat_LISST <- map_dfr( df_list_LISST,~ read_csv(.x))
dat_LISST <- purrr::map_dfr(
  df_list_LISST,
  ~{
    fname <- basename(.x)
    parts <- strsplit(fname, "_")[[1]]
    
    Cruise <- parts[2]
    Transect_ID <- parts[3]
    
    readr::read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        group_id,
        mean_time_Pacific,
        depth_mid,
        Latitude,
        Longitude,
        slope
      )
    ) %>%
      dplyr::mutate(
        Cruise = Cruise,
        Transect_ID = Transect_ID
      )
  }
) %>%
  dplyr::rename(
    time_Pacific = mean_time_Pacific
  )

# merge datasets:
add_nearest_abundance <- function(dat_NBSS, dat_abundance, lat_tol = 0.005, lon_tol = 0.005, time_tol_min = 6) {
  
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
    
    if (use_tow_id) {
      cand <- cand %>%
        dplyr::filter(tow_id == xi$tow_id)
    }
    
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
      
      for (cc in abundance_cols) {
        xi[[cc]] <- NA_real_
      }
      
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
          dplyr::select(
            time_Pacific_abundance,
            Latitude_abundance,
            Longitude_abundance,
            dplyr::all_of(abundance_cols),
            time_diff_min,
            lat_diff,
            lon_diff
          )
      )
  })
}



#do the merging

dat_merged_list <- purrr::imap(
  dat_NBSS_list,
  ~ add_nearest_abundance(
    dat_NBSS = .x,
    dat_abundance = dat_abundance_list[[.y]],
    lat_tol = 0.005,
    lon_tol = 0.005,
    time_tol_min = 30
  )
)

# diagnostic to determine mismatches between NBSS and abundance data:

save_pairing_diagnostics <- function(dat_merged_list, output_root = "/Users/user/SCRIPTS_backups/NBSS_abundances_pairings", width = 10, height = 6, dpi = 300) {
  
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  
  summary_list <- list()
  
  for (project_name in names(dat_merged_list)) {
    
    dat <- dat_merged_list[[project_name]]
    
    if (!"time_Pacific_abundance" %in% names(dat)) {
      stop("time_Pacific_abundance is missing from merged data for ", project_name)
    }
    
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
      dplyr::mutate(project = project_name,  dplyr::across(dplyr::any_of("tow_id"), as.character), .before = 1)
    
    for (i in seq_len(nrow(group_keys))) {
      
      key <- group_keys[i, ]
      
      dat_plot <- dat %>%
        dplyr::semi_join(key, by = group_cols) %>%
        dplyr::filter(!is.na(time_Pacific), !is.na(depth_mid), !is.na(Longitude)) %>%
        dplyr::arrange(time_Pacific, match_status == "No abundance match")
      
      if (nrow(dat_plot) == 0L) {
        next
      }
      
      cruise <- as.character(key$Cruise[[1]])
      transect <- as.character(key$Transect_ID[[1]])
      tow <- if (use_tow_id) as.character(key$tow_id[[1]]) else NULL
      
      n_total <- nrow(dat_plot)
      n_unmatched <- sum(dat_plot$match_status == "No abundance match")
      pct_unmatched <- 100 * n_unmatched / n_total
      
      title_text <- paste0(project_name, " | ", cruise, " | ", transect, if (use_tow_id) paste0(" | tow ", tow) else "")
      subtitle_text <- paste0("NBSS samples: ", n_total, " | Unmatched: ", n_unmatched, " (", round(pct_unmatched, 1), "%)")
      
      p_time <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = time_Pacific, y = depth_mid, color = match_status)) +
        ggplot2::geom_point(size = 2, alpha = 0.8) +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_color_manual(values = c("Matched" = "grey40", "No abundance match" = "red")) +
        ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d\n%H:%M") +
        ggplot2::labs(title = title_text, subtitle = subtitle_text, x = "Time Pacific", y = "Depth midpoint (m)", color = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
      
      p_lon <- ggplot2::ggplot(dat_plot, ggplot2::aes(x = Longitude, y = depth_mid, color = match_status)) +
        ggplot2::geom_point(size = 2, alpha = 0.8) +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_color_manual(values = c("Matched" = "grey40", "No abundance match" = "red")) +
        ggplot2::labs(title = title_text, subtitle = subtitle_text, x = "Longitude", y = "Depth midpoint (m)", color = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "top")
      
      filename_time <- paste0(project_name, "_", cruise, "_", transect, if (use_tow_id) paste0("_tow_", tow) else "", "_NBSS_abundance_pairing_TIME.png")
      filename_lon <- paste0(project_name, "_", cruise, "_", transect, if (use_tow_id) paste0("_tow_", tow) else "", "_NBSS_abundance_pairing_LONGITUDE.png")
      
      filename_time <- gsub("[^A-Za-z0-9_.-]", "_", filename_time)
      filename_lon <- gsub("[^A-Za-z0-9_.-]", "_", filename_lon)
      
      ggplot2::ggsave(filename = file.path(project_dir, filename_time), plot = p_time, width = width, height = height, dpi = dpi)
      ggplot2::ggsave(filename = file.path(project_dir, filename_lon), plot = p_lon, width = width, height = height, dpi = dpi)
    }
    
    readr::write_csv(project_summary, file.path(project_dir, paste0(project_name, "_pairing_summary.csv")))
  }
  
  pairing_summary <- dplyr::bind_rows(summary_list)
  
  readr::write_csv(pairing_summary, file.path(output_root, "all_projects_pairing_summary.csv"))
  
  return(pairing_summary)
}

#run the function:
pairing_summary <- save_pairing_diagnostics(
  dat_merged_list = dat_merged_list,
  output_root = "/Users/user/SCRIPTS_backups/NBSS_abundances_pairings_allData"
)

# check for unique pairings, ALSO USEFUL TO KNOW HOW MANY SAMPLES WE GOT FROM EACH PROJECT
purrr::imap_dfr(dat_merged_list, function(dat, project_name) {
  tibble::tibble(
    project = project_name,
    n_rows = nrow(dat),
    n_group_ids = dplyr::n_distinct(dat$group_id),
    duplicated_group_ids = sum(duplicated(dat$group_id)),
    group_id_is_unique = nrow(dat) == dplyr::n_distinct(dat$group_id)
  )
})

#retain original merged dataset
dat_merged_list_unfiltered <- dat_merged_list

#quantify how much data we are removing:
purrr::imap_dfr(dat_merged_list_unfiltered, function(dat, project_name) {
  tibble::tibble(
    project = project_name,
    n_NBSS_original = nrow(dat),
    n_without_abundance = sum(is.na(dat$time_Pacific_abundance)),
    n_with_abundance = sum(!is.na(dat$time_Pacific_abundance)),
    percent_removed = 100 * n_without_abundance / n_NBSS_original
  )
})

# remove Nas
dat_merged_list <- purrr::map(dat_merged_list, ~ .x %>% dplyr::filter(!is.na(time_Pacific_abundance)))

pairing_summary <- save_pairing_diagnostics(
  dat_merged_list = dat_merged_list,
  output_root = "/Users/user/SCRIPTS_backups/NBSS_abundances_pairings_ONLYPaired"
)

#explore column names
column_info_by_project <- purrr::imap_dfr(dat_merged_list, function(dat, project_name) {
  
  tibble::tibble(
    column_name = names(dat),
    project = project_name,
    column_order = seq_along(names(dat))
  )
}) %>%
  tidyr::pivot_wider(
    names_from = project,
    values_from = column_order,
    names_glue = "{project}_order"
  ) %>%
  dplyr::mutate(
    MEZCAL = dplyr::if_else(!is.na(`MEZCAL_order`), column_name, NA_character_),
    `MBON_ISIIS-2` = dplyr::if_else(!is.na(`MBON_ISIIS-2_order`), column_name, NA_character_),
    `MBON_ISIIS-3` = dplyr::if_else(!is.na(`MBON_ISIIS-3_order`), column_name, NA_character_),
    SPECTRA = dplyr::if_else(!is.na(`SPECTRA_order`), column_name, NA_character_)
  ) %>%
  dplyr::select(column_name, MEZCAL, `MEZCAL_order`,`MBON_ISIIS-2`, `MBON_ISIIS-2_order`, `MBON_ISIIS-3`, `MBON_ISIIS-3_order`, SPECTRA, `SPECTRA_order`)

#Collapse  abundance data by functional group, THIS IS REALLY IMPORTANT FOR ANALYSIS!

# first, extract taxa column limits
project_names <- names(dat_merged_list)
taxa_anchor_list <- purrr::map(project_names, function(project_name) {
  
  name_col <- project_name
  order_col <- paste0(project_name, "_order")
  
  taxa_cols <- column_info_by_project %>%
    dplyr::filter(!is.na(.data[[name_col]])) %>%
    dplyr::filter(grepl("^ind_m3", .data[[name_col]])) %>%
    dplyr::filter(!is.na(.data[[order_col]]))
  
  if (nrow(taxa_cols) == 0L) {
    stop("No ind_m3 columns found for project ", project_name)
  }
  
  start_row <- taxa_cols %>%
    dplyr::slice_min(order_by = .data[[order_col]], n = 1, with_ties = FALSE)
  
  end_row <- taxa_cols %>%
    dplyr::slice_max(order_by = .data[[order_col]], n = 1, with_ties = FALSE)
  
  c(
    start = start_row[[name_col]][[1]],
    end = end_row[[name_col]][[1]]
  )
}) %>%
  rlang::set_names(project_names)

# now, implement the function
dat_merged_list_collapsed <- purrr::imap( dat_merged_list, function(dat, project_name) {
    
    anchors <- taxa_anchor_list[[project_name]]
    
    collapse_taxa_cols(df = dat, start_name = anchors[["start"]], end_name = anchors[["end"]], detritus_handling = "remove")
  }
)



# Random Forest Modelling 
# First, assign shelf_loc_200

dat_merged_list <- purrr::map(
  dat_merged_list,
  ~ .x %>%
    assign_shelf_location(
      "~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA/ancillary/shelf_break_200.csv",
      lon_col_name = Longitude
    )
)


dat_merged_list_collapsed <- purrr::map(
  dat_merged_list_collapsed,
  ~ .x %>%
    assign_shelf_location(
      "~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA/ancillary/shelf_break_200.csv",
      lon_col_name = Longitude
    )
)



#compiled dataset

dat_merged_list_collapsed <- purrr::map(
  dat_merged_list_collapsed, ~ .x %>% dplyr::mutate(dplyr::across(dplyr::any_of("tow_id"), as.character))) # standardize tow_id column

dat_model_all <- dplyr::bind_rows(dat_merged_list_collapsed) %>%
  dplyr::mutate(
    project = as.factor(project),
    Cruise = as.factor(Cruise),
    Transect_ID = as.factor(Transect_ID),
    shelf_loc_200 = as.factor(shelf_loc_200),
    dplyr::across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA_real_))
  )

predictor_vars <- c( "shelf_loc_200", "depth_mid", "Temperature", "Salinity", "Oxygen_ml_l", "chl_a_ug_l", "protist",  "diatom", "copepod", "crustacean", "predatory_gelatinous", "grazer_gelatinous", "fish", "other") # distance from shelf break (0 is shelf break)

setdiff(predictor_vars, names(dat_model_all))

#scale the abundance data due to strong differences in the total abundances across taxa (make a plot too to show this)
dat_model_all <- dat_model_all %>%
  dplyr::mutate(ISIIS_version = dplyr::case_when(project %in% c("MEZCAL", "MBON_ISIIS-2") ~ "ISIIS-2", project %in% c("SPECTRA", "MBON_ISIIS-3") ~ "ISIIS-3", TRUE ~ NA_character_))


plankton_vars <- c("protist",  "diatom", "copepod", "crustacean", "predatory_gelatinous", "grazer_gelatinous", "fish", "other")

dat_model_all <- dat_model_all %>%
  dplyr::group_by(ISIIS_version) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(plankton_vars), ~ as.numeric(scale(.x)), .names = "{.col}_scaled")) %>%
  dplyr::ungroup()

#necessary due to error in defining functional groups: all phytoplankton in ISIIS is diatoms
dat_model_all <- dat_model_all %>%
  dplyr::select(-phytoplankton, -dplyr::any_of("phytoplankton_scaled"))


##### Pre-analysis to determine if isntrument is affecting the community composition PCA with hellinger transformed abundances ####
# Functional-group abundance columns used to define community composition
community_vars <- c("protist",  "diatom", "copepod", "crustacean", "predatory_gelatinous", "grazer_gelatinous", "fish", "other")

# Check that all community variables exist
setdiff(community_vars, names(dat_model_all))

#create community data object
community_data <- dat_model_all %>%
  dplyr::select(group_id, project, ISIIS_version, dplyr::all_of(community_vars))

#Hellinger transform and PCA
community_data <- dat_model_all %>%
  dplyr::select(group_id, project, ISIIS_version, dplyr::all_of(community_vars))

community_complete <- community_data %>%
  tidyr::drop_na(dplyr::all_of(community_vars)) %>%
  dplyr::filter(rowSums(dplyr::across(dplyr::all_of(community_vars))) > 0)

community_matrix <- community_complete %>%
  dplyr::select(dplyr::all_of(community_vars))

community_hellinger <- vegan::decostand(community_matrix, method = "hellinger")

community_pca <- stats::prcomp(community_hellinger, center = TRUE, scale. = FALSE)

#get PCA variance
pca_variance <- community_pca$sdev^2 / sum(community_pca$sdev^2) * 100

pca_variance

#attach scores to metadata
pca_scores <- as.data.frame(community_pca$x) %>%
  dplyr::bind_cols(community_complete %>% dplyr::select(group_id, project, ISIIS_version))

pca_scores <- as.data.frame(community_pca$x) %>%
  dplyr::bind_cols( community_complete %>% dplyr::select(group_id, project, ISIIS_version) ) %>%
  dplyr::left_join( dat_model_all %>% dplyr::select(project, group_id, Cruise, Transect_ID, tow_id, time_Pacific, depth_mid, Latitude, Longitude, slope, intercept), by = c("project", "group_id"))


ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2, color = ISIIS_version)) +
  ggplot2::geom_point(alpha = 0.35, size = 1.5) +
  ggplot2::labs(x = paste0("PC1 (", round(pca_variance[1], 1), "%)"), y = paste0("PC2 (", round(pca_variance[2], 1), "%)"), color = "Instrument", title = "PCA of Hellinger-transformed functional community composition") +
  ggplot2::theme_bw()

ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2, color = project)) +
  ggplot2::geom_point(alpha = 0.35, size = 1.5) +
  ggplot2::labs(x = paste0("PC1 (", round(pca_variance[1], 1), "%)"), y = paste0("PC2 (", round(pca_variance[2], 1), "%)"), color = "Project", title = "PCA of Hellinger-transformed functional community composition") +
  ggplot2::theme_bw()


#get loadings 
pca_loadings <- as.data.frame(community_pca$rotation) %>%
  tibble::rownames_to_column("functional_group")

pca_loadings %>%
  dplyr::select(functional_group, PC1, PC2) %>%
  dplyr::arrange(dplyr::desc(abs(PC1)))

pca_loadings %>%
  dplyr::select(functional_group, PC1, PC2) %>%
  dplyr::arrange(dplyr::desc(abs(PC2)))

#PERMANOVA to determine the effect of ISIIS_version on community structure
#subsampling needed, too much data

set.seed(88)
community_complete_perm <- community_complete %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::group_by(ISIIS_version) %>%
  dplyr::slice_sample(n = 5000) %>%
  dplyr::ungroup()

community_hellinger_perm <- community_hellinger[community_complete_perm$row_id, , drop = FALSE]

permanova_ISIIS <- vegan::adonis2(
  community_hellinger_perm ~ ISIIS_version,
  data = community_complete_perm,
  method = "euclidean",
  permutations = 999
)

permanova_ISIIS


#PERMDISP
community_dist <- stats::dist(community_hellinger_perm, method = "euclidean")

dispersion_ISIIS <- vegan::betadisper(
  community_dist,
  group = community_complete_perm$ISIIS_version
)

anova(dispersion_ISIIS)

permutest(
  dispersion_ISIIS,
  permutations = 999
)

boxplot(
  dispersion_ISIIS,
  xlab = "ISIIS version",
  ylab = "Distance to group centroid"
)


#####.    random model function and ALE plots. ####
fit_rf_ale <- function(dat, response, predictors, seed = 88, folds = 5, tune_length = 5) {
  
  model_data <- dat %>%
    dplyr::select(dplyr::all_of(c(predictors, response))) %>%
    tidyr::drop_na()
  
  message("Response: ", response)
  message("Rows used: ", nrow(model_data))
  message("Rows removed because of NA: ", nrow(dat) - nrow(model_data))
  
  set.seed(seed)
  
  factor_level_check <- model_data %>%
    dplyr::select(where(is.factor)) %>%
    purrr::map_dfr(
      ~ tibble::tibble(
        n_levels_defined = nlevels(.x),
        n_levels_present = dplyr::n_distinct(.x),
        levels_present = paste(unique(as.character(.x)), collapse = ", ")
      ),
      .id = "variable"
    )
  
  print(factor_level_check)
  
  
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
  
  X <- model_data %>%
    dplyr::select(-dplyr::all_of(response))
  
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
  
  list(
    model_data = model_data,
    model = rf_model,
    variable_importance = var_imp,
    predictor = pred,
    ale = ale_data
  )
}


#run model
rf_slope_results <- fit_rf_ale(
  dat = dat_model_all,
  response = "slope",
  predictors = predictor_vars
)





model_data <- df_PSS_full %>%
  mutate(across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA_real_))) %>%
  mutate(Transect_ID = as.factor(Transect_ID),
         season = as.factor(season)) %>%
  select(all_of(predictors), slope) %>%
  na.omit()

# important additions


set.seed(88)
ctrl <- trainControl(method = "cv", number = 5)

rf_slope<- train(
  slope ~ .,
  data       = model_data,
  method     = "ranger",
  trControl  = ctrl,
  tuneLength = 5,
  importance = "permutation",
  metric     = "RMSE"
)
var_imp <- varImp(rf_slope)


# 5) Wrap model for iml
X <- model_data %>% select(-slope)
pred <- Predictor$new(
  model = rf_slope,
  data  = X,
  y     = model_data$slope,
  predict.function = function(m, newdata) predict(m, newdata)  # robust for caret
)

# 6) Helper to compute ALE and return its data
ale_data <- function(feature) {
  fe <- FeatureEffect$new(pred, feature = feature, method = "ale")
  out <- fe$results
  out$feature <- feature
  out
}

# 7) Examples: ALE for Temp, rel_dist, Depth
ale_Temp <- ale_data("Temp")
ale_protist <- ale_data("protist")
ale_Oxy <- ale_data("Oxygen")
ale_chl <- ale_data("chl.ug.l")
#

#without spatiotemporal variables:



# Plotting (prefer raw ALE line; add smooth if you like)
ggplot(ale_dist, aes(x = rel_dist, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "rel_dist", y = "ALE (effect on slope)")

ale_transect$Transect_ID <- fct_reorder(ale_transect$Transect_ID, ale_transect$.value)

ggplot(ale_transect, aes(x = Transect_ID, y = .value)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "Transect", y = "ALE (effect on slope)")

ale_season$season<- fct_reorder(ale_season$season, ale_season$.value)

ggplot(ale_season, aes(x = season, y = .value)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "season", y = "ALE (effect on slope)")

ggplot(ale_protist, aes(x = protist, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "protist", y = "ALE")

ggplot(ale_Oxy, aes(x = Oxygen, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Oxygen", y = "ALE")

ggplot(ale_chl, aes(x = chl.ug.l, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Chla", y = "ALE")

ggplot(ale_pred, aes(x = predatory_gelatinous, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "pred_jelly", y = "ALE")

ggplot(ale_crustacean, aes(x = crustacean, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "crustacean", y = "ALE")

ggplot(ale_other, aes(x = other, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "other", y = "ALE")

ggplot(ale_Temp, aes(x = Temp, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Temp", y = "ALE")

#Now intercept


model_data <- df_collapsed %>%
  select(all_of(predictors), intercept) %>%
  na.omit()

set.seed(88)
ctrl <- trainControl(method = "cv", number = 5)

rf_intercept<- train(
  intercept ~ .,
  data       = model_data,
  method     = "ranger",
  trControl  = ctrl,
  tuneLength = 5,
  importance = "permutation",
  metric     = "RMSE"
)
var_imp_intercept <- varImp(rf_intercept)


# 5) Wrap model for iml
X <- model_data %>% select(-intercept)
pred <- Predictor$new(
  model = rf_intercept,
  data  = X,
  y     = model_data$intercept,
  predict.function = function(m, newdata) predict(m, newdata)  # robust for caret
)

# 6) Helper to compute ALE and return its data
ale_data <- function(feature) {
  fe <- FeatureEffect$new(pred, feature = feature, method = "ale")
  out <- fe$results
  out$feature <- feature
  out
}

# 7) Examples: ALE for Temp, rel_dist, Depth
ale_Temp<- ale_data("Temp")
ale_O2 <- ale_data("Oxygen")
ale_Salinity <- ale_data("Salinity")
ale_chl <- ale_data("chl.ug.l")

# Plotting (prefer raw ALE line; add smooth if you like)
ggplot(ale_Temp, aes(x = Temp, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Temp", y = "ALE (effect on intercept)")

ggplot(ale_Salinity, aes(x = Salinity, y = .value)) +
  geom_line() +
  scale_x_reverse() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Salinity", y = "ALE")

ggplot(ale_O2, aes(x = Oxygen, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Oxygen", y = "ALE")

ggplot(ale_chl, aes(x = chl.ug.l, y = .value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Chla", y = "ALE")



# CODE AFTER THIS ROW IS LEFTOVER, REMOVE IF THERE IS NOTHING USEFUL

df_collapsed_MBON <- collapse_taxa_cols(df_combined_MBON,  start_name = "appendicularian", end_name   = "tunicate_doliolid")

df_collapsed_MBON <-df_collapsed_MBON  %>% assign_relative_distance(., file.path(base_dir, 'ancillary', 'transect_shelf_land_MBON.csv'))

#drop some irrelevant columns:
df_collapsed_MEZCAL <- df_collapsed_MEZCAL %>%
  dplyr::select(-any_of(c("down.up", "Wing.Fin.Angle")))

df_collapsed_MBON <- df_collapsed_MBON %>%
  dplyr::select(-any_of(c("down.up", "Wing.Fin.Angle")))


# SPECTRA DATA
df_list <- list.files(file.path(""), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)

# merge and plotting!
df_PSS_full <- bind_rows(
  df_collapsed_MEZCAL,
  df_collapsed_MBON
) %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))


# isolines every 0.1 across your slope range
df<- df_PSS_full[(df_PSS_full$Transect_ID=='LP' & df_PSS_full$season=="Fall"),]


# --- data prep (rename depth column here if yours is "Depth")
df <- df %>%
  filter(!is.na(lon_mid), !is.na(depth_mid), !is.na(slope)) %>%
  select(lon_mid, depth_mid, slope)

# 0.1-step breaks across your slope range
brks <- seq(floor(min(df$slope, na.rm=TRUE)*10)/10,
            ceiling(max(df$slope, na.rm=TRUE)*10)/10, by = 0.1)

# Interpolate to a grid for contouring
surf <- MBA::mba.surf(as.matrix(df), no.X = 150, no.Y = 150, extend = TRUE)$xyz.est
grid <- expand.grid(lon_mid = surf$x, depth_mid = surf$y)
grid$slope <- as.vector(surf$z)

ggplot() +
  # colored points (original data)
  geom_point(data = df, aes(lon_mid, depth_mid, color = slope), alpha = 0.7, size = 1) +
  # isolines at 0.1 increments, colored by their level
  geom_contour(data = grid, aes(lon_mid, depth_mid, z = slope, color = after_stat(level)),
               breaks = brks, linewidth = 0.4, alpha = 0.9) +
  scale_color_viridis_c(name = "Slope") +
  scale_y_reverse() +
  labs(x = "Longitude", y = "Depth") +
  theme_classic(base_size = 12)



# Coarse spatiotemporal data
df_plot <- df_PSS_full %>%
  filter(!is.na(season), !is.na(shelf_loc_200), !is.na(slope), !is.na(depth_mid)) %>%
  mutate(
    season = factor(season, levels = c("Spring","Summer","Fall","Winter")),
    # put on/off shelf in a nice order if both exist
    shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf","off_shelf")),
    depth_val = as.numeric(depth_mid)   # use continuous depth coloring
  )

ggplot(df_plot, aes(x = season, y = slope)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  geom_point(aes(color = depth_val),
             position = position_jitter(width = 0.15, height = 0, seed = 1),
             alpha = 0.6, size = 1.3) +
  scale_color_distiller(palette = "Blues", direction = -1, name = "Depth (m)") +
  labs(x = "Season", y = "Slope") +
  facet_grid(rows = vars(shelf_loc_200), drop = FALSE) +
  theme_classic(base_size = 12)


#tricky plot for showing relative contributions of each taxon to total biomass per size class:
df_PSS_taxa_MEZCAL<- read.csv("/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_MEZCAL/data/size_spectra/by_taxa/ISIIS_NBSS_by_taxa_chunk1.csv")

# 1) Pick first combo and filter (your code)
first_combo <- df_PSS_taxa_MEZCAL %>%
  distinct(Transect_name, season, depth_mid, lat_mid, lon_mid) %>% slice(1)

df_sub <- df_PSS_taxa_MEZCAL %>%
  filter(
    Transect_name == first_combo$Transect_name,
    season        == first_combo$season,
    depth_mid     == first_combo$depth_mid,
    lat_mid       == first_combo$lat_mid,
    lon_mid       == first_combo$lon_mid
  )

# 2) Sum *NB* (linear space) per size bin & taxon
plot_df <- df_sub %>%
  group_by(biovol_um3_mid, predicted_taxon) %>%
  summarise(NB_sum = sum(NB, na.rm = TRUE), .groups = "drop")

x_levels <- sort(unique(plot_df$biovol_um3_mid))

# Base plot (stacked totals)
min_pos <- min(plot_df$NB_sum[plot_df$NB_sum > 0], na.rm = TRUE)
max_pos <- max(plot_df$NB_sum, na.rm = TRUE)
brks <- 10^seq(floor(log10(min_pos)), ceiling(log10(max_pos)), by = 1)

ggplot(plot_df,
            aes(x = factor(biovol_um3_mid, levels = x_levels),
                y = NB_sum, fill = predicted_taxon)) +
  geom_col(width = 0.8) +
  labs(
    x = expression(paste("Biovolume (", mu, m^3, ") bin mid")),
    y = "Total NB",
    fill = "Predicted taxon"
  ) +
  scale_y_reverse() +
  scale_y_log10() +
  theme_classic(base_size = 12) 


NB_by_size <- plot_df %>%
  dplyr::group_by(biovol_um3_mid) %>%
  dplyr::summarise(NB_sum = sum(NB_sum, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(biovol_um3_mid)

NB_by_size
NB_by_size$logNB<- log10(NB_by_size$NB_sum)
