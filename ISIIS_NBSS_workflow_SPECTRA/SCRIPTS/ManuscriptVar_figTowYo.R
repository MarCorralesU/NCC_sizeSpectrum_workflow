# ============================================================
# ManuscriptVar_figTowYo.R
# Isolated from 4_Kriging_new_size_spectra.R
# Creates and saves only the binned_plot figure.
# ============================================================

# ============================
# ADAPTED KRIGING PIPELINE (your w23_NH "dat")
# keeps original functionality; only maps columns & paths
# ============================
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
library(lubridate)
library(zoo)
library(FNN)
library(colorspace)
base_dir<- getwd()

# --- 0) Paths & identifiers (ADAPT: uses your base_dir and w23_NH) ----
# main directory for outputs (pick your SPECTRA workflow root)
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


option_list <- list(
  make_option(c("-c", "--Cruise"), type = "character", help = "Cruise: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t", "--Transect"), type = "character", help = "Cruise: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-i", "--instrument"), type = "character", default = "ISIIS", help = "LISST-ISIIS, ISIIS or LISST", metavar = "character"),
  make_option(c("-f", "--planktonFraction"), type = "character", default = "allParticles", help = "fraction of the plankton to be processed: allParticles, detritus, plankton", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"),
  make_option(c("-a", "--AkimaExtrapolation"), action = "store_true", default = FALSE, help = "whether to perfom Akima interpolation  [TRUE if present]"))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$Cruise)) {
  stop("Error: --Cruise is required.")
}

#variables to load data:
Cruise <- opt$Cruise # Cruise<- 's22'
Transect <- opt$Transect # Transect<- 'NH'
planktonFraction <- opt$planktonFraction # planktonFraction<- 'plankton'

size_var<- opt$sizeVar # size_var<-'contour'
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "ancillary", "transect_shelf_land.csv")
instrument <- opt$instrument #instrument <- 'ISIIS'

if (instrument=='ISIIS'){
    sizeSpectraDir <- file.path(base_dir, 'data', paste0('size_spectra_', planktonFraction, '_5-10m-bins_thresholded'), size_var, Cruise, Transect) #use this for testing sizeSpectraDir <- file.path(base_dir, 'ISIIS_NBSS_workflow_SPECTRA', 'data', 'size_spectra', 'contour', 'w23', 'NH')
    df_list <- list.files(sizeSpectraDir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
    dat <- map_dfr(
      df_list,
      ~ read_csv(
        .x,
        show_col_types = FALSE,
        col_select = c(
          group_id, total_biovolume, biovol_um3_mid, biovol_um3_width,
          time_Pacific, rel_dist, depth_mid, Latitude, Longitude,
          slope, intercept, dplyr::any_of(env_vars)
        )
      )
    ) %>%
      arrange(time_Pacific) %>%
      dplyr::rename(group_id_old = group_id) %>%
      dplyr::mutate(
        group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
        total_biovolume  = as.numeric(total_biovolume),
        biovol_um3_mid   = as.numeric(biovol_um3_mid),
        biovol_um3_width = as.numeric(biovol_um3_width),
        slope            = as.numeric(slope),
        intercept        = as.numeric(intercept),
        lower_size_edge  = biovol_um3_mid - biovol_um3_width / 2,
        upper_size_edge  = biovol_um3_mid + biovol_um3_width / 2
      ) %>%
      dplyr::group_by(group_id) %>%
      dplyr::summarise(
        time_Pacific = dplyr::first(time_Pacific),
        rel_dist     = dplyr::first(rel_dist),
        depth_mid    = dplyr::first(depth_mid),
        Latitude     = dplyr::first(Latitude),
        Longitude    = dplyr::first(Longitude),
        slope        = dplyr::first(slope),
        intercept    = dplyr::first(intercept),
        total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
        min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
        max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
        dplyr::across(dplyr::any_of(env_vars), dplyr::first),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        total_biovolume_sum = as.numeric(total_biovolume_sum),
        min_size = as.numeric(min_size),
        max_size = as.numeric(max_size),
        slope_anomaly = slope + 1
      ) %>%
      dplyr::filter(!is.na(slope) & slope < 0) %>%
      arrange(time_Pacific)
    
    dat_ord <- dat %>% arrange(Longitude, time_Pacific)  # break ties by time if present

} else if(instrument=='LISST-ISIIS'){
  df_list<- file.path(base_dir, 'data/ISIIS-LISST_mergedNBSS_5m_allParticles' ,paste(Cruise, Transect, 'LISST-ISIIS_NBSS.csv', sep='_'))
  dat <- map_dfr(
    df_list,
    ~ read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        group_id, biovolume_um3, biovol_um3_mid, biovol_um3_width,
        time_Pacific, rel_dist, depth_mid, Latitude, Longitude,
        slope, intercept, dplyr::any_of(env_vars)
      )
    )
  ) %>%
    arrange(time_Pacific) %>%
    dplyr::rename(group_id_old = group_id) %>%
    dplyr::mutate(
      group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
      biovolume_um3  = as.numeric(biovolume_um3),
      biovol_um3_mid   = as.numeric(biovol_um3_mid),
      biovol_um3_width = as.numeric(biovol_um3_width),
      slope            = as.numeric(slope),
      intercept        = as.numeric(intercept),
      lower_size_edge  = biovol_um3_mid - biovol_um3_width / 2,
      upper_size_edge  = biovol_um3_mid + biovol_um3_width / 2
    ) %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(
      time_Pacific = dplyr::first(time_Pacific),
      rel_dist     = dplyr::first(rel_dist),
      depth_mid    = dplyr::first(depth_mid),
      Latitude     = dplyr::first(Latitude),
      Longitude    = dplyr::first(Longitude),
      slope        = dplyr::first(slope),
      intercept    = dplyr::first(intercept),
      total_biovolume_sum = sum(biovolume_um3, na.rm = TRUE),
      min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
      max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
      dplyr::across(dplyr::any_of(env_vars), dplyr::first),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      total_biovolume_sum = as.numeric(total_biovolume_sum),
      min_size = as.numeric(min_size),
      max_size = as.numeric(max_size),
      slope_anomaly = slope + 1
    ) %>%
    dplyr::filter(!is.na(slope) & slope < 0) %>%
    arrange(time_Pacific)
  
  dat_ord <- dat %>% arrange(Longitude, time_Pacific) # break ties by time if present
  
} else if(instrument=='LISST'){
    df_list<- file.path(base_dir, 'data/LISST_data/LISST_size_spectra_thresholded/5m_depthBin' ,paste('SPECTRA', Cruise, Transect, '5m-depthBin_LISST_size_spectra.csv', sep='_'))
    dat <- map_dfr(
      df_list,
      ~ read_csv(
        .x,
        show_col_types = FALSE,
        col_select = c(
          group_id, concentration_um3_m3, size_class_mid_um3, bin_width_um3,
          mean_time_Pacific, depth_mid, Latitude, Longitude,
          slope, intercept, dplyr::any_of(env_vars)
        )
      )
    ) %>%
      dplyr::rename(time_Pacific = mean_time_Pacific,
                    total_biovolume = concentration_um3_m3,
                    biovol_um3_mid = size_class_mid_um3,
                    biovol_um3_width = bin_width_um3) %>%
      arrange(time_Pacific) %>%
      dplyr::rename(group_id_old = group_id) %>%
      dplyr::mutate(
        group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
        total_biovolume  = as.numeric(total_biovolume),
        biovol_um3_mid   = as.numeric(biovol_um3_mid),
        biovol_um3_width = as.numeric(biovol_um3_width),
        slope            = as.numeric(slope),
        intercept        = as.numeric(intercept),
        lower_size_edge  = biovol_um3_mid - biovol_um3_width / 2,
        upper_size_edge  = biovol_um3_mid + biovol_um3_width / 2
      ) %>%
      dplyr::group_by(group_id) %>%
      dplyr::summarise(
        time_Pacific = dplyr::first(time_Pacific),
        depth_mid    = dplyr::first(depth_mid),
        Latitude     = dplyr::first(Latitude),
        Longitude    = dplyr::first(Longitude),
        slope        = dplyr::first(slope),
        intercept    = dplyr::first(intercept),
        total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
        min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
        max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
        dplyr::across(dplyr::any_of(env_vars), dplyr::first),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        total_biovolume_sum = as.numeric(total_biovolume_sum),
        min_size = as.numeric(min_size),
        max_size = as.numeric(max_size),
        slope_anomaly = slope + 1
      ) %>%
      dplyr::filter(!is.na(slope) & slope < 0) %>%
      arrange(time_Pacific)
    
    dat_ord <- dat %>% arrange(Longitude, time_Pacific)
    #filter out LISST data with very high slopes
    dat_ord<- dat_ord[dat_ord$slope<1,]
}

dat_ord$Transect_ID <- Transect

dat_ord <- dat_ord %>% assign_relative_distance(land_ref_path = shore_path, lat_col = "Latitude", long_col = "Longitude")

# create absolute distance for the transect (0 is the easternmost start of the transect)
dat_ord <- dat_ord %>%
  dplyr::mutate(
    distance_from_shelf_signed_km = sign(rel_dist) * distance_sample_km
  )
  


#function to get distance from shore 
get_distance_sample_axis <- function(df) {
  
  axis_dat <- df %>%
    dplyr::filter(
      is.finite(Longitude),
      is.finite(distance_from_shelf_signed_km)
    ) %>%
    dplyr::distinct(
      Longitude,
      distance_from_shelf_signed_km
    ) %>%
    dplyr::mutate(
      distance_from_shore_km =
        max(distance_from_shelf_signed_km, na.rm = TRUE) -
        distance_from_shelf_signed_km
    )
  
  if (nrow(axis_dat) < 2) {
    return(NULL)
  }
  
  dist_fit <- stats::lm(
    distance_from_shore_km ~ Longitude,
    data = axis_dat
  )
  
  dist_intercept <- unname(stats::coef(dist_fit)[1])
  
  dist_slope <- unname(stats::coef(dist_fit)[2])
  
  ggplot2::sec_axis(
    transform = ~ dist_intercept + dist_slope * .,
    name = "km"
  )
}


binned_plot <- ggplot(dat_ord, aes(x = Longitude, y = depth_mid)) +
  
  geom_path(
    aes(group = 1),
    color = "black",
    linewidth = 0.6,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  
  geom_point(
    aes(color = slope),
    alpha = 0.8,
    size = 2
  ) +
  
  scale_y_reverse() +
  
  scale_x_continuous(
    sec.axis = get_distance_sample_axis(dat_ord)
  ) +
  
  colorspace::scale_color_continuous_divergingx(
    palette = "RdBu",
    mid = -1,
    rev = TRUE,
    limits = c(-2, 0),
    name = "Slope"
  ) +
  
  labs(x = "Longitude", y = "Depth (m)", color = "Slope") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x.top = element_text(size = 12),
    axis.ticks.x.top = element_line(),
    axis.title.x.top = element_text(
      size = 12,
      hjust = 1
    )
  )



ggsave(filename = file.path(base_dir, "figures", paste(Cruise, Transect, size_var, instrument, planktonFraction, "binned_data_plot.png", sep = "-")),
       plot = binned_plot, width  = 12, height = 6, units  = "in", dpi = 300)
