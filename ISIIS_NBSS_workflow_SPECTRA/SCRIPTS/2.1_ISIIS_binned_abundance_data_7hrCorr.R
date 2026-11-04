# Goal:get 1-m depth binned plankton abundances.

#COMPARISONS WITH  ISIIS-2 AND FIRST RUNS OF  ISIIS-3:
# this script combines the processin in SPECTRA_PL_2_0_PROCESS_FT_AND_GROUPING.R and SPECTRA_PL_3_BIN_DATA.R
# its core functioning is to get plankton abundance estimates (ind m^3) in variable depth bins (user defined, default value= 1 m)
# This script is ISIIS-3 specific, since it contains a process to combine data from the top and bottom cameras:
#PROCESSING FOLLOWS THESE STEPS:
# 1) reclassify plankton to unknown according to loess probability thresholding (filter_low_clas_scores) and discard
# 2) group data by frame id, get counts by frame id (based on time stamps)
# 3) Combine counts from CamBottom and CamTop
# 4) Get plankton concentration considering info from both camerars (total volume per time step is vol camera*2))
# 5) Apply correction factors 
# 6) assign other relevant variables to group data ( currently on_shelf/off_shelf and day/night/crepuscular)

# THIS SCRIPT CURRENTLY DOES NOT:
# get distance from shore estimates, or standardizes distances from shore in any way
# assigns cast number, or labels of 'down cast' and 'up cast'



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
  make_option(c("-i", "--ISIISver"), type = "integer", default = 2, help = "ISIIS version, either 2 or 3", metavar = "integer"),
  make_option(c("-v", "--dataVersion"), type = "character", default = "MEZCAL", help = "the version of the data, either MEZCAL (use this for MBON) or SPECTRA", metavar = "character"),
  make_option(c("-t", "--transectLats"), type = "character", default = "ancillary/transect_lats.csv", help = "Path to CSV with reference transect latitudes", metavar = "character"),
  make_option(c("-b", "--shelfBreak"), type = "character", default = "ancillary/shelf_break_200.csv", help = "Path to CSV with shelf break longitudes per transect", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"),
  make_option(c("-a", "--addCrepuscular"), action = "store_true", default = FALSE, help = "add crepuscular label to data [TRUE if present]"),
  make_option(c("-s", "--dataForSizeSpectra"), action  = "store_true", default = FALSE, help = "crepuscular label to data [TRUE if present]"),
  make_option(c("-f", "--intermediateFilesExist"), action = "store_true", default = FALSE, help = "whether intermediate files have already been generated [TRUE if present]")
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
combined_dir <- opt$combinedDir #use this for testing combined_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/combined_data_2Cams/w22/NH"
file_name <- gsub("-combined_data", "", tools::file_path_sans_ext(basename(combined_dir)))
ISIISver <- opt$ISIISver
project <- toupper(opt$dataVersion)  # Make it case-insensitive
transect_lat_path <- file.path(base_dir, opt$transectLats) # transect_lat_path <- "/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/ancillary/transect_lats.csv"
shelf_break_path <- file.path(base_dir, opt$shelfBreak) # shelf_break_path <- "/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/ancillary/shelf_break_200.csv"
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- "/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/ancillary/transect_shelf_land.csv"
add_crepuscular <- opt$addCrepuscular
data_sizeSpectra <- opt$dataForSizeSpectra
intermediate_files <- opt$intermediateFilesExist

message("intermediateFilesExist parsed as: ", capture.output(str(opt$intermediateFilesExist)))

df_list <- list.files(file.path(combined_dir), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)
df_list <- df_list[ !grepl("cast_id_depth_file_list", df_list, fixed = TRUE) ]
print(df_list)
# Extract Transect (last folder)
Transect <- basename(combined_dir)
# Extract Cruise (the folder above)
Cruise <- basename(dirname(combined_dir))

output_dir_intermediate <- file.path(base_dir, "data", "intermediate_summaries", "binned_abundance",  Cruise, Transect)

if (!dir_exists(output_dir_intermediate)){
  dir_create(output_dir_intermediate, recurse = TRUE)
  message("Created directory: ", output_dir_intermediate)
}

chunk_size <- 1e6
#message("Processing: ", df_list[1:50])

if (!intermediate_files){
    for (d in df_list){
       message("Processing: ", basename(d))
       #read and apply 
       combined <- read_csv(d) %>% 
    			       mutate(tar_time_id = str_replace(str_replace(tar_id, "^(CamTop|CamBottom)-", ""), "\\.[0-9]+$", ""),
    				     # time = with_tz(as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC"), tzone = "America/Los_Angeles"),
                                     time = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") - hours(7),				      
                                     frame_id = str_remove(img, "_crop_.*$"),
                                     biovolume_um3 = biovol_from_pixels_func(., size_var= 'contour', ISIIS_ver=ISIISver)) %>%
    			       filter_low_class_scores(., base_dir,  project=project) %>%
    			       apply_calibration_factors(.) %>%
                               mutate(across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA))) %>%
    			       rename(time_Pacific = time) %>% # this is needed for the interplation
                                   rename_with(~ str_replace_all(., " ", "_")) %>%
                                   interpolate_missing_values(., c("Temperature", "Conductivity", "Pressure", "Depth", "Salinity","Sound_Velocity", "Horizontal_Speed", "Vertical_Speed", 
    								"Chlorophyll_Wavelength", "Oxygen_ml_l", "pH", "PAR", "Latitude", "Longitude", "chl_a_ug_l")) 
       
       if (Cruise=='w22'){
	   message("Cruise is 'w22' — removing rows where Depth == 27.435850242")
           combined <- combined %>% dplyr::filter(Depth != 27.435850242)
       }    
      
       if (nrow(combined) == 0) {
           message("No rows remaining after filtering for faulty Depth values in w22 — skipping to next iteration.")
       next
       }

 
       if (data_sizeSpectra){
           d_for_binning <- gsub("combined_data_2Cams", "combined_data_2Cams_forBinning", d)
           
           # Extract directory path
           target_dir <- dirname(d_for_binning)
           
           # Create directory recursively if it doesn't exist
           if (!dir.exists(target_dir)) {
             dir.create(target_dir, recursive = TRUE)
             }
           
           # Now save your modified data
           readr::write_csv(combined, d_for_binning)
           
           message("Saved to: ", d_for_binning)
           }
    
       # create a key for the number of frames per depth bin and tar time id. Needs to be done separately from the rest of the summarization 
       n_frames_key <- combined %>%
                       depth_bin_calc(., depth_bin_width=1) %>%
                       group_by(time_Pacific, depth_mid) %>%
                       repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
                       group_by(time_Pacific, depth_mid) %>%
                       summarise( n_frame_id = n_distinct(frame_id), .groups = "drop") %>%
                       repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
                       mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths))) %>%
                       group_by(time_Pacific, depth_mid, depth_run_id) %>%
                       summarise(n_frames = sum(n_frame_id), .groups = "drop") %>%

                       # Group by depth_run_id to get secondse
                       group_by(depth_run_id) %>%
                       mutate(seconds_in_depth_bin = as.numeric(difftime(max(time_Pacific), min(time_Pacific), units = "secs"))) %>%
                       group_by(depth_run_id) %>%
                       summarize(
			 mean_time_Pacific = time_Pacific[ceiling(n()/2)], #mean(time_Pacific, na.rm = TRUE),
			 depth_mid = unique(depth_mid)[1],
			 n_frames = sum(n_frames),
			 seconds_in_depth_bin = unique(seconds_in_depth_bin)[1], .groups = 'drop') %>%
                       ungroup()
                     
    
       # same summarization, but for the environmental variables and counts per taxa (if we get to count the number of frames here, if there were zero organisms for a given  predicted taxon, that frame will not be counted, and abundances will get artificially increased by dividing the counts by a lower sampling volume)
       combined_summarized <- combined %>% 
                               depth_bin_calc(., depth_bin_width=1) %>%
                               group_by(time_Pacific, depth_mid, group) %>%
         			   summarise(
       		                     Cruise_name = unique(Cruise_name)[1],
       			   	     #mean_area = mean(area, na.rm = TRUE),
       			   	     #mean_major = mean(major, na.rm = TRUE),
       			   	     #mean_minor = mean(minor, na.rm = TRUE),
       			   	     #mean_perimeter = mean(perimeter, na.rm = TRUE),
       			   	     #mean_x = mean(x, na.rm = TRUE),
       			   	     #mean_y = mean(y, na.rm = TRUE),
       			   	     #mean_mean = mean(mean, na.rm = TRUE),
       			   	     #mean_height = mean(height, na.rm = TRUE),
       			   	     #mean_probability = mean(probability, na.rm = TRUE),
       			   	     mean_Temperature = mean(Temperature, na.rm = TRUE),
       			   	     mean_Conductivity = mean(Conductivity, na.rm = TRUE),
       			   	     mean_Pressure = mean(Pressure, na.rm = TRUE),
       			   	     mean_Depth = mean(Depth, na.rm = TRUE),
       			   	     mean_Salinity = mean(Salinity, na.rm = TRUE),
       			   	     mean_Sound_Velocity = mean(Sound_Velocity, na.rm = TRUE),
       			   	     mean_Horizontal_Speed = mean(Horizontal_Speed, na.rm = TRUE),
       			   	     mean_Vertical_Speed = mean(Vertical_Speed, na.rm = TRUE),
       			   	     mean_Chlorophyll_Wavelength = mean(Chlorophyll_Wavelength, na.rm = TRUE),
       			   	     mean_Oxygen_ml_l = mean(Oxygen_ml_l, na.rm = TRUE),
       			   	     mean_pH = mean(pH, na.rm = TRUE),
       			   	     mean_PAR = mean(PAR, na.rm = TRUE),
       			   	     mean_Latitude = median(Latitude, na.rm = TRUE),
       			   	     mean_Longitude = median(Longitude, na.rm = TRUE),
       			   	     #mean_total_pixels = mean(total_pixels, na.rm = TRUE),
       			   	     #mean_light = mean(light, na.rm = TRUE),
  				     sum_biovolume_um3 = sum(biovolume_um3, na.rm = TRUE),
  				     mean_biovolume_um3 = mean(biovolume_um3, na.rm = TRUE),
       			   	     mean_dark = mean(dark, na.rm = TRUE),
       			   	     mean_contour = mean(contour, na.rm = TRUE),
       			   	     mean_major_ellipse_axis = mean(major_ellipse_axis, na.rm = TRUE),
       			   	     mean_minor_ellipse_axis = mean(minor_ellipse_axis, na.rm = TRUE),
       			   	     mean_ellipse_area = mean(ellipse_area, na.rm = TRUE),
       			   	     #filter_prob = unique(filter_prob)[1],
       			   	     mean_chl_a_ug_l = mean(chl_a_ug_l, na.rm = TRUE),
                                         n_obs = n(), # THIS IS WHERE THE COUNTS PER TAXA ARE OBTAINED
       			   	     .groups = "drop") 
   
      
        # Step 1: Select only the columns you want group-specific
        group_specific <- combined_summarized %>%
          select(time_Pacific, depth_mid, group,
                 n_obs,
                 sum_biovolume_um3,
                 mean_biovolume_um3,
                 mean_dark,
                 mean_contour,
                 mean_major_ellipse_axis,
                 mean_minor_ellipse_axis,
                 mean_ellipse_area)
        
        # Step 2: Pivot to wide format
        wide_summarized <- group_specific %>%
          pivot_wider(
            id_cols = c(time_Pacific, depth_mid),
            names_from = group,
            names_sep = "_",
            names_glue = "{.value}_{group}",  # e.g., mean_dark_appendicularian
            values_from = c(n_obs, sum_biovolume_um3, mean_biovolume_um3, mean_dark, mean_contour,
                            mean_major_ellipse_axis, mean_minor_ellipse_axis, mean_ellipse_area),
            values_fill = 0
          )
        
         # 2. Environmental summary (just one row per time/depth)
        env_summary <- combined_summarized %>%
          select(-group, -n_obs, -sum_biovolume_um3) %>%
          group_by(time_Pacific, depth_mid) %>%
          summarize(
            Cruise_name = unique(Cruise_name)[1],  # safe pick
            across(-Cruise_name, ~ mean(.x, na.rm = TRUE)),
            .groups = "drop"
          )
        
        
        
        # 3. Merge them back
        wide_summarized <- left_join(env_summary, wide_summarized, by = c("time_Pacific", "depth_mid"))
        
        # Apply the consecutive depth identifier:
        wide_summarized <- wide_summarized %>%
          arrange(time_Pacific) %>%
          repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") %>%
          mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths))) %>%
          relocate(depth_run_id, .before = 1)
        
        
        # Capture the column names
        all_col_names <- names(wide_summarized)
        cols_unique <- c("depth_mid", "Cruise_name")
        cols_range <- all_col_names[5:24]
        cols_pattern <- grep("dark|contour|ellipse", all_col_names, value = TRUE)
        cols_mean <- unique(c("time_Pacific", cols_range, cols_pattern))
        cols_sum <- setdiff(all_col_names, c(cols_unique, cols_mean, "depth_run_id"))
        
        # apply summarization with the collected column names:
        wide_summarized_grouped <- wide_summarized %>%
          group_by(depth_run_id) %>%
          summarize(
            across(all_of(cols_unique), ~ unique(.x)[1]),
        
            mean_time_Pacific = time_Pacific[ceiling(n() / 2)],
        
            across(
              all_of(setdiff(cols_mean, "time_Pacific")),
              ~ mean(.x, na.rm = TRUE)
            ),
        
            across(all_of(cols_sum), ~ sum(.x, na.rm = TRUE)),
            .groups = "drop"
          )

 
      #wide_summarized_grouped <- n_frames_key  %>% left_join( wide_summarized_grouped, by = c("depth_run_id", "mean_time_Pacific", "depth_mid"))
    
      wide_summarized_grouped <- n_frames_key %>% left_join(wide_summarized_grouped, by = c("depth_run_id", "mean_time_Pacific", "depth_mid"), suffix = c("_frames", "_env"))

      # Save per-file summary as an intermediate product
      intermediate_path <- file.path(output_dir_intermediate, paste0(gsub("\\.csv$", "", basename(d)), "_summary.csv"))
      write_csv(wide_summarized_grouped, intermediate_path)
    
      rm(combined, combined_summarized)
      gc()
    }
}

# when mixing intermediate files, make sure to sort them by mean_Time_Pacific
csv_files_intermediate <- Sys.glob(file.path(output_dir_intermediate, "*.csv"))

# get column names across files
all_colnames <- map(csv_files_intermediate, ~ names(vroom(.x, n_max = 1, show_col_types = FALSE))) %>% reduce(union) 

#use custom function since not all taxa are represented in all the files
#df_full <- map_dfr(csv_files_intermediate, ~ read_all_taxa_cols(.x, all_colnames)) %>%
#                                                rename_with(~ "depth_run_id_timeChunk", .cols = "depth_run_id") %>% 
#						collapse_group_columns(., "crustacean_shrimp_molt") %>%
#						arrange(mean_time_Pacific) %>%
#                                               interpolate_missing_values(c("mean_Temperature", "mean_Conductivity", "mean_Pressure", "mean_Depth", 
#                                                                            "mean_Salinity", "mean_Sound_Velocity", "mean_Horizontal_Speed", "mean_Vertical_Speed",
#                                                                            "mean_Chlorophyll_Wavelength", "mean_Oxygen_ml_l", "mean_pH", "mean_PAR", 
#                                                                             "mean_Latitude", "mean_Longitude", "mean_chl_a_ug_l"), time_col='mean_time_Pacific') %>%
#                                                repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "mean_time_Pacific") %>%
#                                                mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths)),
#						       Cruise_name = paste(project, Cruise, sep = "_"),
#                                                       Transect = Transect) %>%
#                                                relocate(Cruise_name, Transect, depth_run_id, .before = 1)


df_full <- map_dfr(csv_files_intermediate, \(path) { 
						    read_all_taxa_cols(path, all_colnames) %>%
						    rename_with(~ "depth_run_id_timeChunk", .cols = any_of( "depth_run_id")) %>%
						    mutate(across(any_of("depth_run_id_timeChunk"), ~ as.integer(as.character(.x))),
                                                    # >>> Normalize mean_time_Pacific to POSIXct in America/Los_Angeles
                                                #    across(any_of("mean_time_Pacific"), ~ as.POSIXct(.x, tz = "UTC")),
                                                          across(any_of("mean_time_Pacific"), ~ {
                                                            if (inherits(.x, "POSIXt")) {
                                                              # keep the SAME clock time, just set tz to LA (no hour shift)
                                                              force_tz(.x, tzone = "America/Los_Angeles")
                                                            } else {
                                                              # parse character as LA local time (no shift from the string)
                                                              as.POSIXct(.x, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Los_Angeles")
							      }}),
                                                          across(any_of(c("Cruise_name", "Transect")), as.character)) %>%
                                                     mutate(across(.cols = -any_of(c("mean_time_Pacific", "Cruise_name", "Transect")),
								    .fns  = ~ suppressWarnings(as.numeric(.x)))
                                                     )}) %>%

                                               collapse_group_columns(., "crustacean_shrimp_molt") %>% 
                                               arrange(mean_time_Pacific)

                                               # ---- Diagnostic BEFORE interpolation ----
                                               if (any(is.na(df_full$mean_Longitude))) {
                                                 warning(
                                                   sum(is.na(df_full$mean_Longitude)),
                                                   " rows have NA mean_Longitude BEFORE interpolation."
                                                 ) 
                                               
                                               df_full %>%
                                                  filter(is.na(mean_Longitude)) %>%
                                                  select(any_of(c("mean_time_Pacific","depth_mid", "depth_run_id_timeChunk",
                                                    "mean_Latitude",
                                                    "mean_Longitude"
                                                  ))) %>%
                                                  print(n = Inf)
                                              }
                                              
                                              if (any(is.na(df_full$mean_Latitude))) {
                                                warning(
                                                  sum(is.na(df_full$mean_Latitude)),
                                                  " rows have NA mean_Latitude BEFORE interpolation."
                                                )
                                              
                                                df_full %>%
                                                  filter(is.na(mean_Latitude)) %>%
                                                  select(any_of(c(
                                                    "mean_time_Pacific",
                                                    "depth_mid",
                                                    "depth_run_id_timeChunk",
                                                    "mean_Latitude",
                                                    "mean_Longitude"
                                                  ))) %>%
                                                  print(n = Inf)
                                              }


                                                df_full<- df_full %>% interpolate_missing_values(c("mean_Temperature", "mean_Conductivity", "mean_Pressure", "mean_Depth", 
                                                                            "mean_Salinity", "mean_Sound_Velocity","mean_Horizontal_Speed", "mean_Vertical_Speed",
                                                                            "mean_Chlorophyll_Wavelength", "mean_Oxygen_ml_l", "mean_pH", "mean_PAR", 
                                                                             "mean_Latitude", "mean_Longitude", "mean_chl_a_ug_l"), time_col='mean_time_Pacific') %>% 
                                                repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "mean_time_Pacific") %>% 
                                                mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths)),
                                                       Cruise_name = paste(project, Cruise, sep = "_"),
                                                       Transect = Transect) %>% 
                                                relocate(Cruise_name, Transect, depth_run_id, .before = 1)


# get lists of columns to be summarized in different methods CONTINUE CHECKING HERE
all_col_names <- names(df_full)
cols_unique <- c("depth_mid", "Cruise_name", "Transect")
cols_mean <- grep("mean|dark|contour|ellipse", all_col_names, value = TRUE)
cols_sum <- setdiff(all_col_names, c(cols_unique, cols_mean, "depth_run_id"))

df_full[setdiff(cols_mean, 'mean_time_Pacific')] <- lapply(df_full[setdiff(cols_mean, 'mean_time_Pacific')], \(x) replace(x, is.na(x), 0))
df_full[cols_sum]  <- lapply(df_full[cols_sum],  \(x) replace(x, is.na(x), 0))

df_full <- df_full %>%
    group_by(depth_run_id) %>%
    arrange(mean_time_Pacific, .by_group = TRUE) %>% 
    summarize(
      mean_time_Pacific = mean_time_Pacific[ceiling(n() / 2)],
      across(all_of(cols_unique), ~ unique(.x)[1]),
      #across(all_of(cols_mean), ~ mean(.x, na.rm = TRUE)),
      across(all_of(setdiff(cols_mean, "mean_time_Pacific")), ~ mean(.x, na.rm = TRUE)),
      across(all_of(cols_sum), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop") 

#final round of despiking
df_full <- df_full %>% repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "mean_time_Pacific") %>% 
			mutate(depth_run_id = with(rle(depth_mid), rep(seq_along(values), lengths)))

df_full <- df_full %>% 
    group_by(depth_run_id) %>%
    arrange(mean_time_Pacific, .by_group = TRUE) %>% 
    summarize(
      mean_time_Pacific = mean_time_Pacific[ceiling(n() / 2)],
      across(all_of(cols_unique), ~ unique(.x)[1]),
      #across(all_of(cols_mean), ~ mean(.x, na.rm = TRUE)),
      across(all_of(setdiff(cols_mean, "mean_time_Pacific")), ~ mean(.x, na.rm = TRUE)),
      across(all_of(cols_sum), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop") %>%
    relocate(n_frames, seconds_in_depth_bin, .after = 4)

df_full <- df_full %>%
  dplyr::as_tibble() %>%
  dplyr::arrange(mean_time_Pacific) %>%
  dplyr::mutate(cast_id = assign_cast_id_simple(depth_mid)) %>%
  dplyr::relocate(cast_id, .after = depth_run_id)

#in cubic meters
if (ISIISver == 2){
frame_vol <- 0.00845 
} else if (ISIISver == 3){
frame_vol <- 0.003 # 10 cm * 10 cm * 30 cm 
}  

# now finally, to get the ind per cubic meter:
df_full <- df_full %>%
  # Step 1–2: Create and rename ind_m3_* columns
  mutate(across(contains("n_obs"),
                ~ .x / (n_frames * frame_vol),
                .names = "ind_m3_{.col}")) %>%
  rename_with(~ gsub("^ind_m3_n_obs_", "ind_m3_", .x), starts_with("ind_m3_n_obs_")) %>%
  
  # Step 3: Drop original n_obs_* columns
  select(-matches("^n_obs")) 

#also get biovolume per cubic meter 
df_full <- df_full %>% 
  # Step 1–2: Create and rename biovol_m3_* columns
  mutate(across(contains("sum_biovolume_um3"),
                ~ .x / (n_frames * frame_vol),
                .names = "biovol_m3_{.col}")) %>% 
  rename_with(~ gsub("^biovol_m3_sum_biovolume_um3_", "biovol_m3_", .x), starts_with("biovol_m3_sum_biovolume_um3_")) %>% 
  
  # Step 3: Drop original biovolume_um3_sum_* columns
  select(-matches("^sum_biovolume_um3")) 



#apply wide format correction factors
df_full <- df_full %>% apply_correction_factors_wide(., base_dir,  project=project, binning = 'abundance')
df_full <- df_full %>% apply_correction_factors_wide(., base_dir,  project=project, binning = 'biovol')

# change time format from UTC to PST
df_full <- df_full %>% 
  dplyr::mutate(
    mean_time_Pacific = format(
      lubridate::with_tz(mean_time_Pacific, "America/Los_Angeles"),
      "%Y-%m-%d %H:%M:%S"
    )   
  ) %>% 
  dplyr::arrange(mean_time_Pacific)



## Create save path
output_dir <- file.path(base_dir, "data", "binned_abundance")
if (!dir_exists(output_dir)) {
  dir_create(output_dir, recurse = TRUE)
  message("Created directory: ", output_dir)
}

# save in chunks:
chunk_size <- 1e6
n_rows <- nrow(df_full)
n_chunks <- ceiling(n_rows / chunk_size)

for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, n_rows)
  chunk <- df_full[start:end, ]
  
  chunk_file <- file.path(output_dir, paste0(Transect, "_", Cruise, "_", "binned_abundance_chunk_", i, ".csv"))
  write_csv(chunk, chunk_file)
  message("Saved to: ", chunk_file)
}


#delete intermediate directory:

#if (grepl("/intermediate_summaries/binned_abundance/", output_dir_intermediate) &&
#    dir.exists(output_dir_intermediate)) {
#  
#  unlink(output_dir_intermediate, recursive = TRUE, force = TRUE)
#  message("Deleted intermediate Cruise/Transect folder: ", output_dir_intermediate)
#} else {
#  warning("Skipping deletion: Path check failed for ", output_dir_intermediate)
#}                                                  
