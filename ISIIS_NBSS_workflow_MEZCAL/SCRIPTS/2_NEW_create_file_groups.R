library(dplyr)
library(magrittr)
library(tidyr)
library(NLP)
library(here)
library(hms)
library(lubridate)
library(tidyverse)
library(scales)
library(fs)
library(stringr)
library(optparse)
library(vroom)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-c", "--combinedDir"), type = "character", help = "Path to combined ISIIS data of the merged plus size data ", metavar = "character"),
  make_option(c("--Cruise"), type = "character", default = NULL, help = "MEZCAL cruise code, e.g. s18, w18, s19, w19"),
  make_option(c("--Transect"), type = "character", default = NULL, help = "MEZCAL transect code, e.g. NH or TR"),
  make_option(c("--project"), type = "character", default = "SPECTRA", help = "Project type: SPECTRA or MEZCAL")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# Print and validate
print(opt)

if (is.null(opt$combinedDir)) {
  stop("Error: --combinedDir is required.")
}

#extract arguments


combined_dir <- opt$combinedDir
project <- toupper(opt$project)

extract_mezcal_tow_id <- function(path) {
  parent <- basename(dirname(path))
  parent %>%
    stringr::str_remove("_MEZCAL_merge$") %>%
    stringr::str_remove(paste0("^", Cruise, "_", Transect, "_")) %>%
    stringr::str_remove("[A-Za-z]$")
}



if (project == "MEZCAL") {

  if (is.null(opt$Cruise) || is.null(opt$Transect)) {
    stop("For MEZCAL, --Cruise and --Transect are required.")
  }

  Cruise <- opt$Cruise
  Transect <- opt$Transect

  mezcal_dirs <- list.dirs(combined_dir, recursive = FALSE, full.names = TRUE)

  mezcal_dirs <- mezcal_dirs[
    grepl(
      paste0("^", Cruise, "_", Transect, "_.*_MEZCAL_merge$"),
      basename(mezcal_dirs)
    )
  ]

  if (length(mezcal_dirs) == 0L) {
    stop("No MEZCAL directories found for ", Cruise, " ", Transect)
  }

  message("Found ", length(mezcal_dirs), " MEZCAL directories for ", Cruise, " ", Transect)

  df_list <- list.files(
    mezcal_dirs,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

} else {

  df_list <- list.files(
    file.path(combined_dir),
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  Transect <- basename(combined_dir)
  Cruise <- basename(dirname(combined_dir))
}

df_list <- df_list[ !grepl("cast_id_depth_file_list", df_list, fixed = TRUE) ]

# Keep only files that have at least one data row beyond the header
has_data_rows <- function(path) {
  line_count <- data.table::fread(
    cmd = paste("wc -l <", shQuote(path)),
    header = FALSE,
    showProgress = FALSE
  )[[1]]
  isTRUE(line_count > 1)
}

df_list_all <- df_list
df_list <- df_list[vapply(df_list, has_data_rows, logical(1))]

message("Found ", length(df_list_all), " CSV files total")
message("Keeping ", length(df_list), " CSV files with at least one data row")
message("Skipping ", length(df_list_all) - length(df_list), " empty/header-only CSV files")

#for testing

run_df_groups<- function(df_list, Cruise){
    df_groups<-data.frame()
    for (d in df_list){
           message("Processing: ", basename(d))
           #read and apply 
           combined <- read_csv(d, show_col_types = FALSE) 
           message(paste(names(combined), collapse = ", "))
           combined<- combined %>%
                           rename_with(~ str_replace_all(., " ", "_")) %>%
                           mutate(
                                  Latitude = if ("Lat" %in% names(.)) Lat else Latitude,
                                  Longitude = if ("Long" %in% names(.)) Long else Longitude,
                                  time_Pacific = lubridate::with_tz(
                                    as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"),
                                    tzone = "America/Los_Angeles"),
                                  tar_time_id = if ("tar_id" %in% names(.)) {
                                    str_replace(as.character(tar_id), "\\.[0-9]+$", "")
                                  } else {
                                    tools::file_path_sans_ext(basename(d))
                                  },

                                  frame_id = if ("img" %in% names(.)) {
                                    str_remove(as.character(img), "_crop_.*$")
                                  } else {
                                    tar_time_id
                                  }) %>%
                           mutate(across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA)))   
 
          if (Cruise=="w22"){
              message("Cruise is w22 — removing rows where Depth == 27.435850242")
              combined <- combined %>% dplyr::filter(Depth != 27.435850242)
          }   
    
          if (nrow(combined) == 0) {
              message("No rows remaining after filtering for faulty Depth values in w22 — skipping to next iteration.")
          next
          }   
    
          combined$file_id <- basename(d)
          combined$file_path <- normalizePath(d, mustWork = FALSE)   
	  combined$tow_id <- if (project == "MEZCAL") {
  				extract_mezcal_tow_id(d)
			     } else {
  				paste(Cruise, Transect, sep = "_")
   		  	     }
          
 
          n_frames_key <- combined %>%
                           depth_bin_calc(., depth_bin_width=1) %>%
                           group_by(time_Pacific, depth_mid, Latitude, Longitude, file_id, file_path) %>%
                           repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
                           group_by(time_Pacific, depth_mid, Latitude, Longitude, file_id, file_path, tow_id) %>%
                           summarise( n_frame_id = n_distinct(frame_id), .groups = "drop") %>%
                           repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
                           mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths))) %>%
                           group_by(time_Pacific, depth_mid, depth_run_id, Latitude, Longitude, file_id, file_path, tow_id) %>%
                           summarise(n_frames = sum(n_frame_id), .groups = "drop") %>%
    
                           # Group by depth_run_id to get secondse
                           group_by(depth_run_id) %>%
                           mutate(seconds_in_depth_bin = as.numeric(difftime(max(time_Pacific), min(time_Pacific), units = "secs"))) %>%
                           group_by(depth_run_id) %>%
                           summarize(
                             min_time_Pacific = min(time_Pacific, na.rm = TRUE),
                             max_time_Pacific = max(time_Pacific, na.rm=TRUE),
                             mean_time_Pacific = mean(time_Pacific, na.rm = TRUE),
                             mean_longitude = mean(Longitude, na.rm = TRUE),
                             mean_latitude= mean(Latitude, na.rm = TRUE),
                             depth_mid = unique(depth_mid)[1],
                             file_id = unique(file_id)[1],
                             file_path = unique(file_path)[1],
			     tow_id = unique(tow_id)[1],
                             n_frames = sum(n_frames),
                             seconds_in_depth_bin = unique(seconds_in_depth_bin)[1], .groups = "drop") %>%
                           ungroup()
    
    
          n_frames_key <- n_frames_key %>% repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "mean_time_Pacific") %>% 
                                           mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths))) %>%
                                            group_by(depth_run_id) %>%
                                             summarize(
                                               min_time_Pacific = min(min_time_Pacific, na.rm = TRUE),
					       max_time_Pacific = max(max_time_Pacific, na.rm=TRUE),
                                               mean_time_Pacific = mean(mean_time_Pacific, na.rm = TRUE),
                                               mean_longitude = mean(mean_longitude, na.rm = TRUE),
                                               mean_latitude= mean(mean_latitude, na.rm = TRUE),
                                               depth_mid = unique(depth_mid)[1],
                                               file_id = unique(file_id)[1],
					       file_path = unique(file_path)[1],
					       tow_id = unique(tow_id)[1],
                                               n_frames = sum(n_frames),
                                               seconds_in_depth_bin = sum(seconds_in_depth_bin), 
                                               file_id = unique(file_id)[1], .groups = 'drop') 
          # free memory from 'combined' now that it's no longer needed
          rm(combined)
          gc()
    
            # append result
          df_groups <- dplyr::bind_rows(df_groups, n_frames_key)
          # optionally free intermediate too
          rm(n_frames_key)
          gc()
    }
    
    return(df_groups)
}

df_groups <- run_df_groups(df_list, Cruise = Cruise)

df_groups <- df_groups %>%
  group_by(tow_id) %>%
  group_modify(~ assign_cast_id(
    .x,
    depth_col = "depth_mid",
    time_col = "max_time_Pacific"
  )) %>%
  ungroup() %>%
  mutate(
    cast_id = paste0("tow_", tow_id, "_", cast_id)
  )


df_groups <- df_groups %>% relocate(depth_mid, .after = cast_id)  

# Per-cast overall cutoff time (optional, still useful to keep)
cast_cutoff <- df_groups %>%
  group_by(tow_id, cast_id) %>%
  summarise(cutoff_time = max(max_time_Pacific, na.rm = TRUE), .groups = "drop")

# Per-cast per-file time window
cast_file_windows <- df_groups %>%
  group_by(tow_id, cast_id, file_id, file_path) %>%
  summarise(
    file_start_time = min(min_time_Pacific, na.rm = TRUE),
    file_end_time   = max(max_time_Pacific, na.rm = TRUE),
    .groups = "drop"
  )


# Final long index: each row has a file for a cast + its own time window
cast_file_index <- cast_file_windows %>%
  left_join(cast_cutoff, by = c('tow_id', "cast_id")) %>%
  mutate(
    file_start_iso = format(file_start_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    file_end_iso   = format(file_end_time,   "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    file_start_unix = as.numeric(as.POSIXct(file_start_time, tz = "UTC")),
    file_end_unix   = as.numeric(as.POSIXct(file_end_time,   tz = "UTC")),
    cutoff_iso      = format(cutoff_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    cutoff_unix     = as.numeric(as.POSIXct(cutoff_time, tz = "UTC"))
  ) %>%
  select(
    tow_id, cast_id, file_id,file_path,
    file_start_iso, file_end_iso,
    file_start_unix, file_end_unix,
    cutoff_iso, cutoff_unix
  ) %>%
  arrange(tow_id, cast_id, file_start_unix, file_id)

cast_index_dir <- file.path(base_dir, "data", "cast_id_file_lists")
dir_create(cast_index_dir, recurse = TRUE)

df_groups %>%
  arrange(max_time_Pacific) %>%
  select(tow_id, cast_id, depth_mid, min_time_Pacific, max_time_Pacific, mean_longitude, file_id) %>%
  write_csv(
    file.path(
      cast_index_dir,
      paste("debug", Cruise, Transect, "df_groups_cast_sequence.csv", sep = "_")
    )
  )


p <- ggplot(df_groups,
            aes(max_time_Pacific, depth_mid, color = factor(cast_id))) +
  geom_point(size = 1) +
  scale_y_reverse() +
  theme_bw()



ggsave(
  filename = file.path(
    cast_index_dir,
    paste("debug", Cruise, Transect,
          "df_groups_cast_sequence_time.png",
          sep = "_")
  ),
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

p_l <-ggplot(df_groups,
       aes(mean_longitude, depth_mid, color=factor(cast_id), group=cast_id)) +
    geom_path() +
    geom_point(size=.5) +
    scale_y_reverse()

ggsave(
  filename = file.path(
    cast_index_dir,
    paste("debug", Cruise, Transect,
          "df_groups_cast_sequence_lon.png",
          sep = "_")
  ),  
  plot = p_l,
  width = 10, 
  height = 6,
  dpi = 300 
)


write_csv(
  cast_file_index,
  file.path(
    cast_index_dir,
    paste(Cruise, Transect, "cast_id_depth_file_list.csv", sep = "_")
  )
)
