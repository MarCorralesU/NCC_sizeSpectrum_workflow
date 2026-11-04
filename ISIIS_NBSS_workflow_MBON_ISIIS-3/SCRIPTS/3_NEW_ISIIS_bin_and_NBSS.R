
#Goal: first step toward estimating size spectra from ISIIS, This will be based on Mo's and Elena's previous scripts
library(data.table)
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
library(optparse)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-c", "--combinedDir"), type = "character", help = "Path to combined ISIIS data of the merged plus size data ", metavar = "character"),
  
  make_option(c("-i", "--ISIISver"), type = "integer", default = 2, help = "ISIIS version, either 2 or 3", metavar = "integer"),
  make_option(c("-v", "--dataVersion"), type = "character", default = "MEZCAL", help = "the version of the data, either MEZCAL (use this for MBON) or SPECTRA", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "time_Pacific,Latitude,Longitude,Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-g", "--depthBin"), type = "integer", default = 2, help = "depth bin size to group data", metavar = "integer"),
  make_option(c("-x", "--depthBinTaxa"), type = "integer", default = 5, help = "depth bin size to group data, since data is parsed by taxa, depth bin needs to be bigger to get good counts", metavar = "integer"),
  make_option(c("-p", "--countThres"), type = "double", default = 0.3, help = "maximum count unceratinty allowed following a poisson distribution", metavar = "NUM"),
  make_option(c("-t", "--transectLats"), type = "character", default = "ancillary/transect_lats.csv", help = "Path to CSV with reference transect latitudes", metavar = "character"),
  make_option(c("-b", "--shelfBreak"), type = "character", default = "ancillary/shelf_break_200.csv", help = "Path to CSV with shelf break longitudes per transect", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"),
  make_option(c("-r", "--thresholding"), type = "character", default = "TRUE", help = "whether the size spectra is going to be thresholded, TRUE or FALSE", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$combinedDir)) {
  stop("Error: --combinedDir is required.")
}

# Variables for data management
combined_dir <- normalizePath(opt$combinedDir, mustWork=TRUE) # #use this for testing combined_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/combined_data_2Cams/w22/GH"
df_list <- list.files(file.path(combined_dir), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)
df_list <- df_list[ !grepl("cast_id_depth_file_list", df_list, fixed = TRUE) ]
Tow_ID <- basename(combined_dir)
Transect <- basename(dirname(combined_dir))
Cruise <- basename(dirname(dirname(combined_dir)))
cast_id_index <- file.path(combined_dir, paste(Cruise, Transect, Tow_ID, "cast_id_depth_file_list.csv", sep = "_"))

#checks
if (!Tow_ID %in% c("D", "N")) stop("Unexpected Tow_ID extracted from combinedDir: ", Tow_ID)
if (!file.exists(cast_id_index)) stop("Cast ID index not found: ", cast_id_index)
if (length(df_list) == 0L) stop("No input CSV files found in: ", combined_dir)
message("Processing Cruise = ", Cruise, " | Transect = ", Transect, " | Tow_ID = ", Tow_ID)
message("Cast index: ", cast_id_index)


#processing variables
ISIISver<- opt$ISIISver
project<- toupper(opt$dataVersion)  # Make it case-insensitive
size_var<- opt$sizeVar
env_vars <- strsplit(opt$envVars, ",")[[1]]
depth_bin_width<- opt$depthBin
depth_bin_width_taxa<- opt$depthBinTaxa
count_thres<- opt$countThres
thresholding <- as.logical(opt$thresholding)

if (thresholding == TRUE){
   thres_string<- 'thresholded' 
} else if (thresholding == FALSE){
   thres_string<- 'unthresholded'
}


#use this block for debugging:



#Variables for post-processing
transect_lat_path <- file.path(base_dir, opt$transectLats) # transect_lat_path <- file.path(base_dir, "ancillary/transect_lats.csv")
shelf_break_path <- file.path(base_dir, opt$shelfBreak) # shelf_break_path <- file.path(base_dir, "ancillary/shelf_break_200.csv")
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- "/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/ancillary/transect_shelf_land.csv"


#load cast identifiers
df_castID<- read_csv(cast_id_index)

#checks
required_cast_cols <- c("cast_id", "file_id", "file_start_iso", "file_end_iso")
missing_cast_cols <- setdiff(required_cast_cols, names(df_castID))
if (length(missing_cast_cols) > 0L) stop("Cast index is missing required columns: ", paste(missing_cast_cols, collapse = ", "))

#casts ids
casts <- unique(df_castID$cast_id)

#function to skip NB products if they have no data after all the filtering
process_nb <- function(tbl, kind, cid) {
  if (nrow(tbl) == 0L) {
    message("Cast ", cid, ": ", kind, " is empty after filtering/thresholds — skipping.")
    return(FALSE)
  }
  TRUE
}

# Column list to guide and reduce the amount of data that is loaded
needed_cols <- c(
  "time","img", "tar_id", "Class", "predicted_taxon", "probability",
  "Temperature", "Conductivity", "Pressure", "Depth",  "Salinity",
  "Raw Oxygen [V]", "Raw PAR [V]", "Raw pH [V]", "Chlorophyll", # vars needed for calibration
  "Latitude","Longitude", size_var)


#create dir to save size spectra data
output_dir <- file.path(base_dir, "data", paste("size_spectra_plankton_5-10m-bins", thres_string, sep='_'), size_var, Cruise, Transect, Tow_ID)
if (!dir_exists(output_dir)) {
  dir_create(output_dir, recurse = TRUE)
  message("Created directory: ", output_dir)
}

#processing loop, iterates over each cast
for (cid in casts) { # for debugging: tail(casts, -14) or cid<- casts[18]
  message("Processing ", cid, " ...")
  
  cast_index <- df_castID %>%
  filter(cast_id == cid) %>%
  mutate(
    file_start_time = as.POSIXct(file_start_iso, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    file_end_time   = as.POSIXct(file_end_iso,   format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  )

  ids <- unique(cast_index$file_id)
  files <- df_list[basename(df_list) %in% ids]

  if (length(files) == 0) {
    warning("No files found for ", cid)
  next
  }

  # lookup table from file_id to cast-specific window
  file_windows <- cast_index %>%
     select(file_id, file_start_time, file_end_time)
  
   

   #updated, read one file, filter AND THEN combine
   dat <- map_dfr(files, function(f) { 
				      fname <- basename(f)
                                      win <- file_windows %>% filter(file_id == fname)

                                      if (nrow(win) != 1) {
                                        stop("Expected exactly one window for file ", fname, " in cast ", cid,
                                             ", found ", nrow(win))
                                      }
                                      
                                      file_start <- win$file_start_time[[1]]
                                      file_end   <- win$file_end_time[[1]]

                                      df<-read_csv(f, col_select = any_of(needed_cols), show_col_types = FALSE) %>%
                                          mutate(
                                                 time_Pacific = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") - hours(7),
                                                 tar_time_id  = str_replace(str_replace(tar_id, "^(CamTop|CamBottom)-", ""), "\\.[0-9]+$", ""),
                                                 frame_id     = str_remove(img, "_crop_.*$")
                                                 ) %>%
                                          arrange(time_Pacific) %>%
					  filter(time_Pacific >= file_start, time_Pacific <= file_end)

                                          # --- remove empty dataframe ---

                                          if (nrow(df) == 0) {
                                              message("Skipping empty file: ", basename(f))
                                          return(NULL)
                                          }

					  # --- enforce types so bind_rows won't fail + log NA stats ---
  					  if (!"img" %in% names(df)) {
  					    df$img <- NA_character_
  					    message("File ", fname, ": column 'img' missing — created as NA_character_.")
  					  }
  					  if (!"tar_id" %in% names(df)) {
  					    df$tar_id <- NA_character_
  					    message("File ", fname, ": column 'tar_id' missing — created as NA_character_.")
  					  }

  					  df <- df %>%
  					    mutate(
  					      img             = as.character(img),
  					      tar_id          = as.character(tar_id),
  					      predicted_taxon = as.character(predicted_taxon)
  					    )

  					  # report NA counts (characters only have NA, not NaN)
  					  n_img_na    <- sum(is.na(df$img))
  					  n_tarid_na  <- sum(is.na(df$tar_id))
  					  n_taxa_na   <- sum(is.na(df$predicted_taxon))

  					  if (n_img_na   > 0) message("File ", fname, ": 'img' has ", n_img_na, " NA values.")
  					  if (n_tarid_na > 0) message("File ", fname, ": 'tar_id' has ", n_tarid_na, " NA values.")
  					  if (n_taxa_na  > 0) message("File ", fname, ": 'predicted_taxon' has ", n_taxa_na, " NA values.")
  					  # ------------------------------------------------------------

				          if (nrow(df) == 0L) return(df)  # <-- safe early exit for this file
        				
					  df <- filter_low_class_scores(df, base_dir, project = project)
				          if (nrow(df) == 0L) return(df)  # <-- safe early exit for this file
                                          return(
                                          df %>% filter(!is.na(predicted_taxon),
                                                       !str_detect(predicted_taxon, regex("art.*fact|detritus|carcass|molt", ignore_case = TRUE))) %>% # add |detritus to the regex expression to exclude detritus
                                                 apply_calibration_factors() %>%
                                                mutate(across(where(is.numeric), ~ replace(.x, is.infinite(.x), NA_real_))))
                                      }, .progress = TRUE) %>% arrange(time_Pacific) %>% interpolate_missing_values(c(
                                            "Temperature","Pressure", "Depth", "Conductivity","Salinity",
                                            "Oxygen_ml_l","pH","PAR","Latitude","Longitude","chl_a_ug_l"))
  
  if (Cruise=='w22'){
      message("Cruise is 'w22' — removing rows where Depth == 27.435850242")
      dat <- dat %>% dplyr::filter(Depth != 27.435850242)
  }   

  if (nrow(dat) == 0) {
      message("No rows remaining after filtering for faulty Depth values in w22 — skipping to next iteration.")
  next
  } 

  print(data.frame(head(dat))) 

  # ensure time sorted, then subset
  cast_start <- min(file_windows$file_start_time, na.rm = TRUE)
  cast_end   <- max(file_windows$file_end_time, na.rm = TRUE)
  
  dat <- dat %>% filter(time_Pacific >= cast_start, time_Pacific <= cast_end)
   
  # IF at this point we have the data ONLY for one cast, then we can go ahead and bin by size and depth. ALSO apply correction factor. Lat/Lon and time will come from averaging by depth bin, so that we can retain maximum resolution 
  binned_full<- dat %>% depth_bin_calc(., depth_bin_width=depth_bin_width) %>% size_binning(., size_var=size_var, depth_bin_width=depth_bin_width, ISIIS_ver=ISIISver)  %>% arrange(depth_mid)# after some exploration, depth_bin_width that allowed to have >4 size classes FOR MOST DEPTHS AFTER GETTING ONLY PLANKTON DATA was 2 meters. 

  #can we pass the NB and size spectra functions?
  
  NB_full<-binned_full %>% select(c(time_Pacific, Latitude, Longitude, depth_mid, img, frame_id, predicted_taxon, group,  Temperature, Salinity, Oxygen_ml_l, pH, PAR, chl_a_ug_l, biovol, biovol_binned, biovol_um3_mid, biovol_um3_width)) %>%
                     filter(!is.na(predicted_taxon)) %>%
                     filter(!str_detect(predicted_taxon, regex("art.*fact|unknown|detritus|carcass|molt", ignore_case = TRUE))) %>% #|detritus
                     group_abundance_func(., grouping_vars=c('depth_mid'), env_vars=env_vars, group_by_taxa =FALSE, ISIIS_ver=ISIISver, thresholding = thresholding, count_thres_prob=count_thres) %>%
                     lin_fit_func(., filter_min_size_class=TRUE, n_min_size_class=3)

 
  #now, do the per taxa. It is likely that we will need larger depth integration to get good counts per size class
  binned_taxa<- dat %>% depth_bin_calc(., depth_bin_width=depth_bin_width_taxa) %>% size_binning(., size_var=size_var, depth_bin_width=depth_bin_width_taxa, ISIIS_ver=ISIISver) %>% arrange(depth_mid)# depth bin for taxa product needs to be bigger, to get more counts  

  #can we pass the NB and size spectra functions?

  NB_taxa<-binned_taxa %>% select(c(time_Pacific, Latitude, Longitude, depth_mid, img, frame_id, predicted_taxon, group, Temperature, Salinity, Oxygen_ml_l, pH, PAR, chl_a_ug_l, biovol, biovol_binned, biovol_um3_mid, biovol_um3_width)) %>%
                     filter(!is.na(predicted_taxon)) %>%
                     filter(!str_detect(predicted_taxon, regex("art.*fact|unknown|detritus|carcass|molt", ignore_case = TRUE))) %>%#|detritus
                     group_abundance_func(., grouping_vars=c('depth_mid'), env_vars=env_vars, group_by_taxa =TRUE, ISIIS_ver=ISIISver, thresholding = thresholding, count_thres_prob=count_thres) %>%
                     lin_fit_func(., filter_min_size_class=TRUE, n_min_size_class=3)


  #add relevant columns to keep track of data origin and other spatiotemporal categorical variables for data analysis
  
  if (process_nb(NB_full, "community", cid)) {
      NB_full$Cruise <- Cruise
      NB_full$tow_id <- Tow_ID
      NB_full$cast_id<- cid
      NB_full <- NB_full  %>% 
                                assign_transect_name(transect_lat_path, lat_col_name="Latitude") %>% 
                                assign_shelf_location(shelf_break_path, lon_col_name='Longitude') %>% 
                                assign_season(time_col = "time_Pacific") %>% 
                                assign_day_night(crepuscular = FALSE, time_col= "time_Pacific", lat_col= "Latitude", long_col= "Longitude", tz = "America/Los_Angeles") %>%
                                assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude') %>%
 				relocate(Cruise, Transect_ID, tow_id, cast_id, .before = 1)
  
  write.csv(NB_full, file.path(output_dir, paste(Cruise, Transect, Tow_ID, cid, 'ISIIS_NBSS_by_community.csv', sep='-')), row.names=FALSE)
  }
 
  if (process_nb(NB_taxa, "taxa", cid)) {
      NB_taxa$Cruise <- Cruise
      NB_taxa$tow_id <- Tow_ID
      NB_taxa$cast_id<- cid
      NB_taxa <- NB_taxa  %>%      
                                assign_transect_name(transect_lat_path, lat_col_name="Latitude") %>%
                                assign_shelf_location(shelf_break_path, lon_col_name='Longitude') %>%
                                assign_season(time_col = "time_Pacific") %>%
                                assign_day_night(crepuscular = FALSE, time_col= "time_Pacific", lat_col= "Latitude", long_col= "Longitude", tz = "America/Los_Angeles") %>%
                                assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude') %>%
 				relocate(Cruise, Transect_ID, tow_id, cast_id, .before = 1)
  
  write.csv(NB_taxa, file.path(output_dir, paste(Cruise, Transect, Tow_ID, cid, 'ISIIS_NBSS_by_taxa.csv', sep='-')), row.names=FALSE)
  }
 
  
  # free memory per iteration
  rm(dat, binned_full, binned_taxa, NB_full, NB_taxa)
  gc()
}


#

