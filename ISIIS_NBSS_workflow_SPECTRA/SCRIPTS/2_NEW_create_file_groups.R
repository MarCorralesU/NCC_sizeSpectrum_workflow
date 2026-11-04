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
  make_option(c("-c", "--combinedDir"), type = "character", help = "Path to combined ISIIS data of the merged plus size data ", metavar = "character")
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
combined_dir <- opt$combinedDir #use this for testing combined_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/combined_data_2Cams/w23/NH"

df_list <- list.files(file.path(combined_dir), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)

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


 
# Extract Transect (last folder)
Transect <- basename(combined_dir)
# Extract Cruise (the folder above)
Cruise <- basename(dirname(combined_dir))

#for testing

run_df_groups<- function(df_list, Cruise){
    df_groups<-data.frame()
    for (d in df_list){
           message("Processing: ", basename(d))
           #read and apply 
           combined <- read_csv(d) %>%
                                   mutate(tar_time_id = str_replace(str_replace(tar_id, "^(CamTop|CamBottom)-", ""), "\\.[0-9]+$", ""),
                                         time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") - hours(7),
                                         frame_id = str_remove(img, "_crop_.*$")) %>%
                                   mutate(across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA))) %>%
                                   rename(time_Pacific = time) %>% # this is needed for the interplation
                                       rename_with(~ str_replace_all(., " ", "_")) 
    
          if (Cruise=="w22"){
              message("Cruise is w22 — removing rows where Depth == 27.435850242")
              combined <- combined %>% dplyr::filter(Depth != 27.435850242)
          }   
    
          if (nrow(combined) == 0) {
              message("No rows remaining after filtering for faulty Depth values in w22 — skipping to next iteration.")
          next
          }   
    
          combined$file_id<-basename(d)
    
          n_frames_key <- combined %>%
                           depth_bin_calc(., depth_bin_width=1) %>%
                           group_by(time_Pacific, depth_mid, Latitude, Longitude, file_id) %>%
                           repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
                           group_by(time_Pacific, depth_mid, Latitude, Longitude, file_id) %>%
                           summarise( n_frame_id = n_distinct(frame_id), .groups = "drop") %>%
                           repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
                           mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths))) %>%
                           group_by(time_Pacific, depth_mid, depth_run_id, Latitude, Longitude, file_id) %>%
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

df_groups <- assign_cast_id(df_groups, depth_col = "depth_mid",
                            time_col = "max_time_Pacific")

df_groups <- df_groups %>% relocate(depth_mid, .after = cast_id)  

# Per-cast overall cutoff time (optional, still useful to keep)
cast_cutoff <- df_groups %>%
  group_by(cast_id) %>%
  summarise(cutoff_time = max(max_time_Pacific, na.rm = TRUE), .groups = "drop")

# Per-cast per-file time window
cast_file_windows <- df_groups %>%
  group_by(cast_id, file_id) %>%
  summarise(
    file_start_time = min(min_time_Pacific, na.rm = TRUE),
    file_end_time   = max(max_time_Pacific, na.rm = TRUE),
    .groups = "drop"
  )

# Final long index: each row has a file for a cast + its own time window
cast_file_index <- cast_file_windows %>%
  left_join(cast_cutoff, by = "cast_id") %>%
  mutate(
    file_start_iso = format(file_start_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    file_end_iso   = format(file_end_time,   "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    file_start_unix = as.numeric(as.POSIXct(file_start_time, tz = "UTC")),
    file_end_unix   = as.numeric(as.POSIXct(file_end_time,   tz = "UTC")),
    cutoff_iso      = format(cutoff_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    cutoff_unix     = as.numeric(as.POSIXct(cutoff_time, tz = "UTC"))
  ) %>%
  select(
    cast_id, file_id,
    file_start_iso, file_end_iso,
    file_start_unix, file_end_unix,
    cutoff_iso, cutoff_unix
  ) %>%
  arrange(cast_id, file_start_unix, file_id)

 write_csv(cast_file_index, file.path(combined_dir, paste(Cruise, Transect, 'cast_id_depth_file_list_UPDATED.csv', sep="_")))
