# summary size spectra 
setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
library(ggplot2)
library(colorRamps)
library(plyr)
library(dplyr)
library(reshape2)
library(akima)
library(gstat)
library(raster)
library(randomForest)   # not used unless you enable the RF section
library(scales)
library(metR)
library(kriging)
library(ggthemes)
library(optparse)
library(purrr)
library(readr)
library(fuzzyjoin)
base_dir<- getwd()

source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))
source(file.path(base_dir, 'SCRIPTS', 'functions', 'model_functions.R')) 


option_list <- list(
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"),
  make_option(c("-b", "--shelfBreak"), type = "character", default = "ancillary/shelf_break_200.csv", help = "Path to CSV with shelf break longitudes per transect", metavar = "character"),
  make_option(c("-t", "--transectLats"), type = "character", default = "ancillary/transect_lats.csv", help = "Path to CSV with  transect latitudes", metavar = "character")
  )

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)


size_var<- opt$sizeVar # size_var<-'contour'
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "ancillary", "transect_shelf_land.csv")
transect_lat_path <- file.path(base_dir, opt$transectLats) # transect_lat_path<- file.path(base_dir, "ancillary/transect_lats.csv")
shelf_break_path <- file.path(base_dir, opt$shelfBreak) # shelf_break_path <- file.path(base_dir, "ancillary/shelf_break_200.csv")
sizeSpectraDir <- file.path(base_dir,  'data', 'size_spectra_plankton_5-10m-bins_thresholded') 
LISSTdir <-  file.path(base_dir,  'data', "LISST_data", "LISST_size_spectra", "5m_depthBin")
binned_abundanceDir <- file.path(base_dir,  'data', 'binned_abundance') 
merge_env_Dir <- file.path(base_dir,  'data', 'merge_environmental') 



output_dir_figs_model <- file.path(base_dir, "figures", 'size_spectra', 'model')
if (!dir.exists(output_dir_figs_model)) {
  dir.create(output_dir_figs_model, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_model)
}

## Load size spectra data
df_list_ISIIS <- list.files(sizeSpectraDir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_LISST <- list.files(LISSTdir, pattern = "5m-depthBin_LISST_size_spectra\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_abundance <- list.files(binned_abundanceDir, pattern = "binned_abundance_chunk_1\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_env <- list.files(merge_env_Dir, pattern = "merged_enviro\\.csv$", full.names = TRUE, recursive = TRUE)







dat_ISIIS <- dat_ISIIS %>%
  # keep the original for reference
  dplyr::rename(group_id_old = group_id) %>%
  # new globally consecutive id: increments whenever group_id_old changes
  dplyr::mutate(
    group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE))
  )
dat_ISIIS <- dat_ISIIS %>%
  mutate(
    shelf_loc_200 = case_when(
      rel_dist >= 0 ~ "on_shelf",
      rel_dist <  0 ~ "off_shelf",
      TRUE          ~ NA_character_
    ),
    shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf"))
  )





#adjust binned aboundance data to spatial scale of ISIIS data
dat_abundance<- map_dfr( df_list_abundance, ~ read_csv( .x, show_col_types = FALSE)) %>%
  mutate(time_Pacific = mean_time_Pacific - hours(1))  %>%
  dplyr::select(-mean_time_Pacific)

colnames(dat_abundance) <- sub("^mean_", "", colnames(dat_abundance))

dat_abundance<- dat_abundance %>% 
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) %>% 
  assign_season(time_col = 'time_Pacific') %>% 
  assign_transect_name(., transect_lat_path, 'Latitude') %>% 
  assign_shelf_location(., shelf_break_path) %>%
  assign_relative_distance( land_ref_path = shore_path, lat_col  = "Latitude", long_col = "Longitude")

dat_abundance<- dat_abundance %>% depth_bin_calc(., depth_bin_width=5)

#define env vars and taxa cols that can be averaged by depth_mid

wanted_cols <- c(
  "time_Pacific",
  "Latitude",
  "Longitude",
  "Depth",
  "Temperature",
  "Salinity",
  "Oxygen_ml_l",
  "pH",
  "chl_a_ug_l",
  "PAR"
  
)

env_cols_spectra <- intersect(wanted_cols, names(dat_abundance))



taxa_cols_spectra <- names(dat_abundance)[grepl("ind_m3", names(dat_abundance))]

df_abundance_summarized <- dat_abundance %>%
  group_by(depth_run_id, cast_id, Cruise_name, Transect_ID, shelf_loc_200, season, depth_mid) %>%
  summarise(
    across(all_of(env_cols_spectra),  ~ mean(as.numeric(.), na.rm = TRUE)),
    across(all_of(taxa_cols_spectra), ~ sum(as.numeric(.),  na.rm = TRUE)),
    .groups = "drop"
  )  %>% rename(Cruise=Cruise_name)

#remove NA, artifact and detritus counts:
df_abundance_summarized <- df_abundance_summarized %>%
  dplyr::select(
    !matches("NA|art[ie]fact", ignore.case = TRUE)
  )
df_abundance_summarized$Cruise <- sub("^SPECTRA_", "", df_abundance_summarized$Cruise)

df_abundance_summarized <- collapse_taxa_cols(df_abundance_summarized,  start_name = "ind_m3_copepod_calanoid_paracalanidae", end_name   = "ind_m3_crustacean_shrimp_molt", ignore_case = TRUE,
                                              detritus_handling = "remove")
keys <- c("Cruise", "Transect_ID","shelf_loc_200","season","depth_mid","lat_mid","lon_mid")

#standardizing time column
dat_LISST  <- dat_LISST  %>% rename(time_Pacific = mean_time_Pacific)

dat_ISIIS <- dat_ISIIS %>%
  dplyr::select(
    -matches("^time$"),        # remove plain 'time'
    -matches("time_of_day")    # remove potential conflicts
  )

# selecting only columns that are needed:
lisst_trim <- dat_LISST %>%
  dplyr::select(
    Cruise, Transect_ID,
    time_Pacific, depth_mid,
    slope_LISST = slope   # rename for clarity
  )

abund_trim <- df_abundance_summarized %>%
  dplyr::select(
    Cruise, Transect_ID,
    time_Pacific, depth_mid,
    protist, phytoplankton, copepod, crustacean,
    predatory_gelatinous, grazer_gelatinous,
    fish, other
  )

tol_sec   <- 30       # seconds allowed difference in timestamp
tol_depth <- 0.1        # meters allowed difference in depth

#fuzzyjoin function
do_fuzzy <- function(I_sub, X_sub, tol_sec, tol_depth) {
  fuzzyjoin::fuzzy_left_join(
    I_sub,
    X_sub,
    by = c(
      time_Pacific = "time_Pacific",
      depth_mid    = "depth_mid"
    ),
    match_fun = list(
      function(a,b) abs(as.numeric(a - b)) <= tol_sec,
      function(a,b) abs(a - b) <= tol_depth
    )
  ) %>%
    mutate(
      dt    = abs(as.numeric(time_Pacific.x - time_Pacific.y)),
      dd    = abs(depth_mid.x - depth_mid.y),
      score = (dt / tol_sec) + (dd / tol_depth)
    )
}

# ---- Time range diagnostics ----
time_range_str <- function(df) {
  if (!nrow(df)) return("NO DATA")
  rng <- range(df$time_Pacific, na.rm = TRUE)
  paste(format(rng[1]), "→", format(rng[2]))
}

# --- helper: best fuzzy match extractor ---
best_match <- function(cand_df) {
  cand_df %>%
    group_by(time_Pacific.x, depth_mid.x) %>%
    slice_min(score, n = 1, with_ties = FALSE) %>%
    ungroup()
}

# --- build final merged dataset ---
merged_final <- purrr::pmap_dfr(
  test_blocks,
  function(Cruise, Transect_ID) {
    
    message("\n=== Processing block ===")
    message("Cruise: ", Cruise, " | Transect: ", Transect_ID)
    
    # subsets (NO .data$, no !!)
    I_sub <- dplyr::filter(dat_ISIIS,
                           Cruise == Cruise,
                           Transect_ID == Transect_ID)
    
    L_sub <- dplyr::filter(lisst_trim,
                           Cruise == Cruise,
                           Transect_ID == Transect_ID)
    
    A_sub <- dplyr::filter(abund_trim,
                           Cruise == Cruise,
                           Transect_ID == Transect_ID)
    
    message("Rows: ISIIS=",nrow(I_sub),
            " | LISST=",nrow(L_sub),
            " | ABUND=",nrow(A_sub))
    
    cand_LI  <- do_fuzzy(I_sub, L_sub, tol_sec, tol_depth)
    
    best_LI  <- best_match(cand_LI) %>%
      dplyr::select(
        Cruise       = Cruise.x,
        Transect_ID  = Transect_ID.x,
        time_Pacific = time_Pacific.x,
        depth_mid    = depth_mid.x,
        slope_LISST
      )
    
    cand_AI <- do_fuzzy(I_sub, A_sub, tol_sec, tol_depth)
    
    best_AI <- best_match(cand_AI) %>%
      dplyr::select(
        Cruise = Cruise.x,
        Transect_ID = Transect_ID.x,
        time_Pacific = time_Pacific.x,
        depth_mid = depth_mid.x,
        protist, phytoplankton, copepod,
        crustacean, predatory_gelatinous,
        grazer_gelatinous, fish, other
      )
    
    final_block <- I_sub %>%
      dplyr::select(Cruise, Transect_ID, time_Pacific, depth_mid) %>%
      dplyr::left_join(best_LI,
                       by = c("Cruise","Transect_ID","time_Pacific","depth_mid")) %>%
      dplyr::left_join(best_AI,
                       by = c("Cruise","Transect_ID","time_Pacific","depth_mid"))
    
    # -----------------------------------------
    # SAVE to disk (with auto directory create)
    # -----------------------------------------
    out_dir <- file.path(base_dir, "data", "data_modeling")
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
      message("Created directory: ", out_dir)
    }
    
    outfile <- paste0(Cruise, "_", Transect_ID, "_merged_for_model.csv")
    fpath <- file.path(out_dir, outfile)
    
    message("Writing → ", fpath)
    readr::write_csv(final_block, fpath)
    
    # return tiny tibble so pmap_dfr succeeds
    tibble::tibble(
      Cruise = Cruise,
      Transect_ID = Transect_ID,
      outfile = fpath
    )
  }
)

keys <- c("Cruise", "Transect_ID","time_Pacific","depth_mid")
# 
# 
# 
# 
# 
 df_combined_SPECTRA <- df_abundance_summarized_MEZCAL %>%
   mutate(Transect_ID = as.character(Transect_ID)) %>%
   left_join(
     df_unique_PSS_MEZCAL %>% mutate(Transect_ID = as.character(Transect_ID)),
    by = keys
) %>%
drop_na(slope)
# 
# # troubleshoot env vars data
# 
#   
#   # Ensure column exists
#   
#   # Ensure output directory exists
#   
#   # Symbol for tidy evaluation
#   
#   #----------------------------------------------
#   # Compute min/max per season (ignoring NA/Inf)
#   #----------------------------------------------
#       # Build label text
#       # Positions for annotation inside facet
#       # (x = max x value, y = 90th percentile height)
#   
#   #----------------------------------------------
#   # Build histogram with facets
#   #----------------------------------------------
#     
#     #----------------------------------------------
#   # Add min/max text to each facet
#   #----------------------------------------------
#   
#   #----------------------------------------------
#   # Save & return plot
#   #----------------------------------------------
#   
#   
# 
# 
#   
# 
# # apply functions to dat_enviro and dat_abundance so that we can make direct comparisons
# 
# 
# 
# 
# 
# 
# # compare the env data in size spectra products to raw env data and binned_abundance
# 
# 
# 
# 