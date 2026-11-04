setwd('~/SCRIPTS_backups')
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
library(lubridate)
library(zoo)
base_dir<- getwd()

# --- 0) Paths & identifiers (ADAPT: uses your base_dir and w23_NH) ----
# main directory for outputs (pick your SPECTRA workflow root)
main_dir <- file.path(base_dir, "ISIIS_NBSS_workflow_SPECTRA")
source(file.path(main_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

option_list <- list(
  make_option(c("-c", "--Cruise"), type = "character", help = "Cruise: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t", "--Transect"), type = "character", help = "Cruise: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character")
  )

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$Cruise)) {
  stop("Error: --Cruise is required.")
}

#variables to load data:
abundanceDir <- file.path(main_dir, 'data', 'binned_abundance')
message("abundance_dir: ", abundanceDir)
Cruise <- opt$Cruise # Cruise<- 'w23'
Transect <- opt$Transect # Transect<- 'GH'
shore_path <- file.path(main_dir, opt$distFromShore) # shore_path  <- file.path(main_dir, "ancillary", "transect_shelf_land.csv")


output_dir_data <- file.path(main_dir, "data", "bathymetry_polygons")
if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_data)
}

# Build a regex pattern: e.g., "^NH_w23_binned_abundance_test_chunk_.*\\.csv$" and list binned abundance files
pat <- sprintf("^%s_%s_binned_abundance_chunk_.*\\.csv$", Transect, Cruise)

df_list <- list.files(abundanceDir, pattern = pat, full.names = TRUE, recursive  = TRUE)


dat <- map_dfr(
  df_list,
  ~ read_csv(
    .x,
    show_col_types = FALSE,
    col_select = c(cast_id, mean_Latitude, mean_Longitude, 
    )
  ) |> distinct()
)


dat_ord <- dat %>% arrange(mean_Longitude)  # break ties by time if present

# ---------- 8) seafloor polygon ----------
# load enviro data
env_paths <- list.files(file.path(base_dir, 'ISIIS_NBSS_workflow_SPECTRA', 'data', 'merge_environmental'), full.names = TRUE, recursive = TRUE)

env_paths <- env_paths[grepl(Cruise, env_paths)]

env_dat <- map_dfr(
  env_paths,
  ~ read_csv(.x, show_col_types = FALSE, col_select = c(Time, Latitude, Longitude, Depth, Altitude)) %>%
    distinct() %>%
    mutate(Altitude = as.numeric(Altitude)) %>%
    filter(!is.na(Altitude)) %>%
    arrange(Longitude)
)

env_dat$Transect_ID<- Transect
env_dat <- env_dat %>% assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude')


env_dat <- env_dat %>% mutate(time_Pacific = as.POSIXct(Time, tz = "UTC"))

env_dat <- env_dat %>%
  dplyr::filter(
    Latitude >= min(dat_ord$mean_Latitude),
    Latitude <= max(dat_ord$mean_Latitude),
    Longitude >= min(dat_ord$mean_Longitude),
    Longitude <= max(dat_ord$mean_Longitude)
  )

# do a fuzzy join of precise Depth and Altimeter values:

#sort env_dat and then apply filling for NaN
env_dat <- env_dat %>% arrange(time_Pacific) %>% 
  interpolate_missing_values(., c( 'Latitude', 'Longitude', 'Depth', 'Altitude')) 

env_dat$bottomDepth<- env_dat$Depth +env_dat$Altitude

dat_ord <- dat_ord %>% rename(Latitude = mean_Latitude) %>% rename(Longitude = mean_Longitude) %>% arrange(Longitude)
env_dat <- env_dat %>% arrange(Longitude)

# define maximum tolerances
time_tol  <- dseconds(10)   # within 10 seconds (adjust as needed)
latlon_tol <- 0.0001          # roughly ~10 m, tweak depending on resolution

message("Adding precise Depth and Altimeter values to draw bathymetry, please wait...")

# perform nearest join
dat_joined <- fuzzy_left_join(
  dat_ord,
  env_dat,
  by = c(
    "Longitude"    = "Longitude"
  ),
  match_fun = list(
    function(x, y) abs(x - y) <= latlon_tol
  )
)

dat_joined <- dat_joined %>%
  dplyr::select(
    cast_id, rel_dist, 
    Latitude = Latitude.x, Longitude = Longitude.x,
    Depth, Altitude, bottomDepth
  )

dat_joined <- dat_joined %>% arrange(Longitude)
dat_joined <- na.omit(dat_joined)
dat_joined$Altitude[dat_joined$rel_dist < 0] <- NA
dat_joined$bottomDepth[dat_joined$rel_dist < 0] <- NA

dat_max_alt <- dat_joined %>%
  group_by(cast_id) %>%
  slice_max(order_by = bottomDepth, n = 1, with_ties = FALSE) %>%
  ungroup()

dat_max_alt <- dat_max_alt %>% arrange(Longitude)

# since 0 is the shelf break, which is 200 m, any altitude readings where rel_dist<0 are bugs
# Find index of rel_dist closest to 0
dat_max_alt$bottomDepth[dat_max_alt$rel_dist < 0] <- 200



# smoothing andf fill up of inf values
dat_max_alt <- dat_max_alt %>%
  arrange(rel_dist)  # smaller (negative) = offshore

# optional: smooth first to remove sharp jumps
loess_fit <- loess(bottomDepth ~ rel_dist, data = dat_max_alt, span = 0.2)
dat_max_alt$B_smooth <- predict(loess_fit, dat_max_alt$rel_dist)

# Clean and limit bathy_df
bathy_df <- dat_max_alt %>%
  dplyr::ungroup() %>%
  dplyr::select(rel_dist, B_smooth, Longitude) %>%
  dplyr::filter(rel_dist > 0) %>%   # keep only on-shelf region
  tidyr::drop_na() %>%
  as.data.frame()

output_filePath <-file.path(output_dir_data, paste(Cruise, Transect, 'bathymetry.csv', sep='-'))
write.csv(bathy_df, output_filePath, row.names=FALSE)

message("saved: ", output_filePath)
