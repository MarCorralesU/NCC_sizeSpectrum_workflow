# ============================
# ADAPTED KRIGING PIPELINE (your w23_NH "dat")
# keeps original functionality; only maps columns & paths
# ============================
setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(akima)
library(gstat)
library(scales)
library(metR)
library(optparse)
library(purrr)
library(readr)
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
Cruise <- opt$Cruise # Cruise<- 'w22'
Transect <- opt$Transect # Transect<- 'GH'
planktonFraction <- opt$planktonFraction # planktonFraction<- 'allParticles'

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
          group_id, total_biovolume, logNB, biovol_um3_mid, biovol_um3_width,
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
        logNB             = as.numeric(logNB),
        biovolume_um3_m3  = as.numeric((10^logNB) * biovol_um3_width),
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
        total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE),
        total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
        min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
        max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
        dplyr::across(dplyr::any_of(env_vars), dplyr::first),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        total_biovolume_um3_m3= as.numeric(total_biovolume_um3_m3),
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
        group_id, biovolume_um3, logNB, biovol_um3_mid, biovol_um3_width,
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
      logNB             = as.numeric(logNB),
      biovolume_um3_m3  = as.numeric((10^logNB) * biovol_um3_width),
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
      total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE),
      min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
      max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
      dplyr::across(dplyr::any_of(env_vars), dplyr::first),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      total_biovolume_um3_m3= as.numeric(total_biovolume_um3_m3),
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
          group_id, concentration_um3_m3,  logNB, size_class_mid_um3, bin_width_um3,
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
        logNB             = as.numeric(logNB),
        biovolume_um3_m3  = as.numeric((10^logNB) * biovol_um3_width),
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
        total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE),
        min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
        max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
        dplyr::across(dplyr::any_of(env_vars), dplyr::first),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        total_biovolume_sum = as.numeric(total_biovolume_sum),
        total_biovolume_um3_m3= as.numeric(total_biovolume_um3_m3),
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

# calculate integrated biovolume
x1 <- 12.78849       # minimum range contemplated by LISST
x2 <- 4.547442e13    # maximum range contemplated by ISIIS

dat_ord <- dat_ord %>%
  dplyr::mutate(
    slope = as.numeric(slope),
    intercept = as.numeric(intercept),
    min_size = as.numeric(min_size),
    max_size = as.numeric(max_size),
    b_lin = 10^intercept,
    
    # integral across the full common instrument range
    integrated_biovol_all = dplyr::if_else(
      is.finite(slope) & is.finite(b_lin),
      ifelse(
        abs(slope + 1) < 1e-12,
        b_lin * log(x2 / x1),
        (b_lin / (slope + 1)) * (x2^(slope + 1) - x1^(slope + 1))
      ),
      NA_real_
    ),
    
    # integral across the sample-specific size range for that group_id
    integrated_biovol_sample = dplyr::if_else(
      is.finite(slope) & is.finite(b_lin) &
        is.finite(min_size) & is.finite(max_size) &
        min_size > 0 & max_size > 0 & max_size >= min_size,
      ifelse(
        abs(slope + 1) < 1e-12,
        b_lin * log(max_size / min_size),
        (b_lin / (slope + 1)) * (max_size^(slope + 1) - min_size^(slope + 1))
      ),
      NA_real_
    )
  )

dat_ord <- dat_ord %>% dplyr::mutate(biovol_diff = integrated_biovol_sample - total_biovolume_um3_m3)

# transform integrated biovolume due to large differences 
dat_ord <- dat_ord %>% mutate(log10_integrated_biovol_all = log10(integrated_biovol_all),
                              log10_total_biovolume_um3_m3 = log10(total_biovolume_um3_m3),
                              log10_integrated_biovol_sample = log10(integrated_biovol_sample),
                              log10_total_biovolume_sum = log10(total_biovolume_sum))



# Create output folders 
interp_tag <- if (opt$AkimaExtrapolation) "interp_and_krige_dist10" else "krige-only" # interp_tag <-"interp_and_krige_dist10"

output_dir_data <- if (instrument == "ISIIS") {
  file.path(base_dir, "data", interp_tag, instrument, planktonFraction, "size_spectra_coeffs")
} else {
  file.path(base_dir, "data", interp_tag, instrument, "size_spectra_coeffs")
}

if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_data)
}

output_dir_figs_normal <- file.path(base_dir, "figures", interp_tag, "size_spectra_coeffs", instrument, Cruise, Transect)

output_dir_figs_normal_diagnostics <- file.path(base_dir, "figures", interp_tag, "size_spectra_coeffs", instrument, Cruise, Transect, "diagnostic_figsAndTables")

if (!dir.exists(output_dir_figs_normal)) {
  dir.create(output_dir_figs_normal, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_normal)
}

if (!dir.exists(output_dir_figs_normal_diagnostics)) {
  dir.create(output_dir_figs_normal_diagnostics, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_normal_diagnostics)
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
    legend.position = "top"
  )



ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, instrument, planktonFraction, "binned_data_plot.png", sep = "-")),
       plot = binned_plot, width  = 12, height = 6, units  = "in", dpi = 300)

if (Cruise == 'w22' && Transect == 'CM' && instrument != 'LISST') {
  #  # Drop all rows where rel_dist < -1
  message('w22 CM data has issues, subsetting')
  dat_ord <- dat_ord %>% 
    dplyr::filter(Longitude >= -124.3) %>% 
    dplyr::filter(Longitude <= -124.1)
}

# ---------- 1) Build myDf from your `dat` (use Longitude -> distance km) ----------
# Expect: dat has at least: time_Pacific, rel_dist, depth_mid, Latitude, Longitude, slope

# We'll compute distance in km from the westernmost longitude (start of transect):
# This preserves the original pipeline's functionality (distance in km along x).

# Use dplyr:: to avoid plyr clashes
myDf <- dat_ord %>%
  dplyr::mutate(
    Transect_ID = Cruise,
    HD_ID       = Transect,
    Depth       = as.numeric(depth_mid),
    Longitude   = as.numeric(Longitude),
    Latitude    = as.numeric(Latitude),
    slope       = as.numeric(slope),
    slope_anomaly = as.numeric(slope_anomaly),
    intercept   = as.numeric(intercept),
    log10_integrated_biovol_all = as.numeric(log10_integrated_biovol_all),
    log10_total_biovolume_um3_m3 = as.numeric(log10_total_biovolume_um3_m3),
    log10_integrated_biovol_sample = as.numeric(log10_integrated_biovol_sample),
    log10_total_biovolume_sum = as.numeric(log10_total_biovolume_sum),
    biovol_diff = as.numeric(biovol_diff)
  ) %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, Latitude, Longitude,
    dplyr::any_of(env_vars),
    slope, slope_anomaly, intercept,
    log10_integrated_biovol_all, log10_total_biovolume_um3_m3, log10_integrated_biovol_sample,
    log10_total_biovolume_sum, biovol_diff, time_Pacific
  )

# Horizontal distance along transect:
# 0 km = easternmost observation; positive values increase westward
lat0 <- mean(myDf$Latitude, na.rm = TRUE)
km_per_deg_lon <- 111.320 * cos(lat0 * pi / 180)

lon0 <- max(myDf$Longitude, na.rm = TRUE)

myDf <- myDf %>%
  dplyr::mutate(
    GPS.dist.from.start.km = (lon0 - Longitude) * km_per_deg_lon
  ) %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, GPS.dist.from.start.km,
    dplyr::any_of(env_vars),
    slope, slope_anomaly, intercept,
    log10_total_biovolume_sum, log10_total_biovolume_um3_m3,
    log10_integrated_biovol_all, log10_integrated_biovol_sample,
    biovol_diff, Latitude, Longitude, time_Pacific
  )

# Drop odd columns if they exist (same as original script)
if ("NA" %in% names(myDf)) myDf$`NA` <- NULL
if ("NA_mean_area"  %in% names(myDf)) myDf$NA_mean_area  <- NULL
if ("NA_mean_major" %in% names(myDf)) myDf$NA_mean_major <- NULL
if ("NA_mean_minor" %in% names(myDf)) myDf$NA_mean_minor <- NULL
gc()



# ------------------------------------------------------------
# Fixed anisotropic spatial coordinates
# 1 km horizontal = 10 m vertical
# ------------------------------------------------------------

anisotropy_m_per_km <- 5 # 10

myDf <- myDf %>%
  dplyr::mutate(
    x_km = GPS.dist.from.start.km,
    z_m = Depth,
    x_scaled = x_km,
    z_scaled = z_m / anisotropy_m_per_km
  )

# ---------- 2) Long variable table (variables to krige; here: slope) ----------
biophys.10 <- myDf %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, GPS.dist.from.start.km,
    dplyr::any_of(c(
      "slope", "slope_anomaly", "intercept",
      "log10_integrated_biovol_all", "log10_total_biovolume_um3_m3",
      "log10_integrated_biovol_sample",
      "log10_total_biovolume_sum", "biovol_diff"
    ))
  )
message("selected columns for krige")

biophys.11 <- reshape2::melt(
  biophys.10,
  id.vars = c("Transect_ID","HD_ID","Depth","GPS.dist.from.start.km")
)

biophys.11 <- biophys.11 %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(PCT = percent_rank(value), ntile = ntile(value, 100))

krige_complete <- NULL
variable1 <- unique(biophys.11$variable)

for (subvar in variable1) {
  di <- plyr::ddply(
    .data = biophys.11[biophys.11$variable == subvar, ],
    .variables = ~ HD_ID,
    function(x) {
      mean_ntil97 <- mean(x$value[x$ntile == 98], na.rm = TRUE)
      x$value[x$PCT > 0.9959] <- as.numeric(mean_ntil97)
      x
    }
  )
  krige_complete <- rbind(krige_complete, di)
}
message("reset outliers")


# ---------- 3) Environmental data (optional; same behavior) ----------
env_keep <- env_vars[env_vars %in% names(myDf)]

biophys.12 <- myDf %>%
  dplyr::select(dplyr::any_of(c("Transect_ID","HD_ID","Depth","GPS.dist.from.start.km", env_keep)))

biophys.13 <- reshape2::melt(
  biophys.12,
  id.vars = c("Transect_ID","Depth","HD_ID","GPS.dist.from.start.km")
)

krige_complete.b <- krige_complete %>%
  dplyr::select(-PCT, -ntile) %>%
  dplyr::mutate(value = as.numeric(value))

biophys.14 <- dplyr::bind_rows(krige_complete.b, biophys.13)

#exclude vars

# # ---------- 4) Coarse interpolation (akima) per variable (xo 0.5 km, yo 5 m) ----------
 variable <- as.character(unique(biophys.14$variable))
 biophys.14$HD_ID <- as.factor(biophys.14$HD_ID)

 exclude_var <- c("camera","Chlorophyll Wavelength","height","horizontal.vel","Sound Velocity","x","y")
 variable <- variable[!variable %in% exclude_var]
 biophys.14 <- biophys.14[biophys.14$variable %in% variable, ]

 
 # ------------------------------------------------------------
 # Spatial interpolation settings
 # ------------------------------------------------------------
 
 krige_dx_km <- 0.5
 krige_dz_m <- 2.5
 akima_dx_km <- 0.5
 akima_dz_m <- 5
 
 
 # ------------------------------------------------------------
 # Local ordinary kriging settings
 # ------------------------------------------------------------
 
 krige_nmax <- 30
 krige_nmin <- 5
 krige_maxdist <- 5 #10
 
 # ------------------------------------------------------------
 # Variogram fitting distances to test
 # ------------------------------------------------------------
 
 variogram_cutoff <- 10
 
 # ------------------------------------------------------------
 # Optional Akima interpolation before kriging
 # ------------------------------------------------------------
 
 biophys_observed <- biophys.14
 
 if (opt$AkimaExtrapolation) {
   
   message("Akima interpolation ENABLED — two-stage interpolation.")
   
   interpol_complete.a <- list()
   
   for (subvar in variable) {
     
     df <- biophys_observed %>%
       dplyr::filter(variable == subvar,
                     is.finite(GPS.dist.from.start.km),
                     is.finite(Depth),
                     is.finite(value))
     
     if (nrow(df) < 3) next
     
     x0 <- seq(min(df$GPS.dist.from.start.km), max(df$GPS.dist.from.start.km), by = akima_dx_km)
     y0 <- seq(min(df$Depth), max(df$Depth), by = akima_dz_m)
     
     if (length(x0) < 2 || length(y0) < 2) next
     
     pli <- tryCatch(
       akima::interp(
         x = df$GPS.dist.from.start.km,
         y = df$Depth,
         z = df$value,
         xo = x0,
         yo = y0,
         duplicate = "mean",
         extrap = FALSE,
         linear = TRUE
       ),
       error = function(e) {
         warning(paste("Akima interpolation failed for", subvar, ":", conditionMessage(e)))
         NULL
       }
     )
     
     if (is.null(pli)) next
     
     interpol_complete.a[[subvar]] <- expand.grid(
       GPS.dist.from.start.km = pli$x,
       Depth = pli$y
     ) %>%
       dplyr::mutate(
         value = as.vector(pli$z),
         variable = subvar,
         HD_ID = Transect,
         Transect_ID = Cruise
       ) %>%
       dplyr::filter(is.finite(value))
   }
   
   biophys.14 <- dplyr::bind_rows(interpol_complete.a)
   
   if (nrow(biophys.14) == 0) {
     stop("Akima interpolation produced no valid observations.")
   }
   
 } else {
   
   message("Akima interpolation DISABLED — using original observations.")
   
   biophys.14 <- biophys_observed
 }
 

 # Get coordinates of observed points
 obs <- myDf %>%
   dplyr::select(GPS.dist.from.start.km, Depth) %>%
   na.omit()


 
 grid <- expand.grid(
   GPS.dist.from.start.km = seq(
     min(myDf$GPS.dist.from.start.km, na.rm = TRUE),
     max(myDf$GPS.dist.from.start.km, na.rm = TRUE),
     by = krige_dx_km
   ),
   Depth = seq(
     0,
     max(myDf$Depth, na.rm = TRUE),
     by = krige_dz_m
   )
 )

 ## Compute nearest-neighbor distances (in grid units)

 ## Choose a reasonable threshold — e.g. 1 km horizontally or 10 m vertically equivalent
 ## Here, since Akima works in the same units for x and y, use ~0.6 to 1.0 as a start

 
 
 # Compute nearest-neighbor distances after scaling depth
 # Here: 10 m vertical distance is treated as roughly equivalent to 1 km horizontal distance
 obs_scaled <- obs %>%
   dplyr::transmute(
     x_scaled = GPS.dist.from.start.km,
     y_scaled = Depth / anisotropy_m_per_km
   )
 
 grid_scaled <- grid %>%
   dplyr::transmute(
     x_scaled = GPS.dist.from.start.km,
     y_scaled = Depth / anisotropy_m_per_km
   )
 
 nn <- FNN::get.knnx(obs_scaled, grid_scaled, k = 1)
 
 grid$dist_to_nearest <- as.numeric(nn$nn.dist[, 1])
 
 # Smaller threshold = more conservative / wider mask
 thr <- 5
 
 grid$no_data <- grid$dist_to_nearest > thr
 

 # prepare  no data mask for plotting:
 mask_slope <- grid %>%
   dplyr::filter(no_data) %>%
   mutate(
     Transect_ID = Transect,
     Longitude = lon0 - GPS.dist.from.start.km / km_per_deg_lon,
     Latitude = mean(myDf$Latitude)
   ) %>%
   assign_relative_distance(land_ref_path = shore_path, lat_col = "Latitude", long_col = "Longitude")

 # save no data mask to reuse later
 mask_out <- mask_slope %>%
   dplyr::transmute(
     Cruise = Cruise,
     Transect_ID = Transect,
     dist_kriged_km = GPS.dist.from.start.km,
     depth_kriged_m = Depth,
     Longitude,
     Latitude,
     rel_dist,
     no_data
   )
 
 print(vapply(mask_out, typeof, character(1)))
 print(vapply(mask_out, is.matrix, logical(1)))
 
 readr::write_csv(
   mask_out,
   file.path(output_dir_data, paste(Cruise, Transect, size_var, instrument, paste0("no_data_mask_", interp_tag, ".csv"), sep = "-"))
 )
 
 
# ---------- 5) Kriging (fine grid; Exp, maxdist=10) ----------
krige_biophys1 <- NULL
variable2 <- unique(biophys.14$variable)

# ------------------------------------------------------------
# Ordinary kriging with fixed anisotropic coordinates
# ------------------------------------------------------------

cv_summary_list <- list()
cv_predictions_list <- list()

for (subvar in variable2) {
  
  x <- biophys.14 %>%
    dplyr::filter(variable == subvar, is.finite(value),
                  is.finite(GPS.dist.from.start.km), is.finite(Depth)) %>%
    dplyr::mutate(
      x_scaled = GPS.dist.from.start.km,
      z_scaled = Depth / anisotropy_m_per_km
    )
  
  if (nrow(x) < 10) next
  
  # Collapse coincident coordinates to avoid singular kriging systems
  x <- x %>%
    dplyr::group_by(x_scaled, z_scaled) %>%
    dplyr::summarise(value = mean(value), .groups = "drop")
  
  if (nrow(x) < 10) next
  
  # Prediction grid: 500 m horizontally, 2.5 m vertically
  x0 <- seq(min(x$x_scaled), max(x$x_scaled), by = krige_dx_km)
  z0 <- seq(0, max(x$z_scaled) * anisotropy_m_per_km, by = krige_dz_m)
  
  d0 <- expand.grid(
    x = x0,
    y = z0 / anisotropy_m_per_km
  )
  
  sp::coordinates(d0) <- ~ x + y
  sp::coordinates(x) <- ~ x_scaled + z_scaled
  
 
  
  # Empirical variogram and fitted exponential model
  x_var <- gstat::variogram(
    value ~ 1,
    x,
    cutoff = variogram_cutoff,
    width = 1
  )
  
  x_fit <- tryCatch(
    gstat::fit.variogram(x_var, model = gstat::vgm("Exp")),
    error = function(e) NULL
  )
  
  if (is.null(x_fit) || any(!is.finite(x_fit$psill)) ||
      any(!is.finite(x_fit$range)) || any(x_fit$psill < 0) ||
      any(x_fit$range[-1] <= 0)) {
    warning(paste("Invalid variogram fit:", subvar))
    next
  }
  
  # ------------------------------------------------------------
  # Save empirical and fitted exponential variogram
  # ------------------------------------------------------------
  
  variogram_dir <- file.path(
    output_dir_figs_normal_diagnostics,
    "variograms"
  )
  
  dir.create(
    variogram_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Empirical variogram values
  readr::write_csv(
    as.data.frame(x_var),
    file.path(variogram_dir, paste0(subvar, "_empirical_variogram.csv"))
  )
  
  # Fitted variogram parameters
  readr::write_csv(
    as.data.frame(x_fit),
    file.path(variogram_dir, paste0(subvar, "_fitted_variogram.csv"))
  )
  
  # Empirical variogram points and fitted exponential curve
  p_variogram <- ggplot2::ggplot(
    as.data.frame(x_var),
    ggplot2::aes(x = dist, y = gamma)
  ) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_line(
      data = as.data.frame(
        gstat::variogramLine(
          x_fit,
          maxdist = max(x_var$dist, na.rm = TRUE),
          n = 300
        )
      ),
      ggplot2::aes(x = dist, y = gamma),
      color = "red",
      linewidth = 1
    ) +
    ggplot2::labs(
      title = paste(Cruise, Transect, subvar),
      subtitle = "Empirical variogram and fitted exponential model",
      x = "Separation distance (anisotropically scaled units)",
      y = "Semivariance"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  ggplot2::ggsave(
    filename = file.path(
      variogram_dir,
      paste0(subvar, "_variogram.png")
    ),
    plot = p_variogram,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300
  )
  
  # Ordinary kriging
  x_kriged <- tryCatch(
    gstat::krige(
      value ~ 1, x, d0,
      model = x_fit,
      nmax = krige_nmax,
      nmin = krige_nmin,
      maxdist = krige_maxdist
    ),
    error = function(e) {
      warning(paste("Kriging failed for", subvar, ":", conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(x_kriged)) next
  
  x_kriged_df <- as.data.frame(x_kriged)
  
  # Restore physical coordinates for plotting
  x_kriged_df$Depth <- x_kriged_df$y * anisotropy_m_per_km
  x_kriged_df$GPS.dist.from.start.km <- x_kriged_df$x
  x_kriged_df$Longitude <- lon0 - x_kriged_df$x / km_per_deg_lon
  
  x_kriged_df$variable <- subvar
  x_kriged_df$HD_ID <- Transect
  x_kriged_df$Transect_ID <- Cruise
  
  krige_biophys1 <- rbind(krige_biophys1, x_kriged_df)
  
  # ----------------------------------------------------------
  # Cross-validation: leave-one-observation-out
  # ----------------------------------------------------------
  
  cv <- tryCatch(
    gstat::krige.cv(
      value ~ 1, x,
      model = x_fit,
      nfold = nrow(x),
      nmax = krige_nmax,
      nmin = krige_nmin,
      maxdist = krige_maxdist
    ),
    error = function(e) {
      warning(paste("Cross-validation failed for", subvar, ":", conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(cv)) next
  
  cv_df <- as.data.frame(cv) %>%
    dplyr::mutate(
      variable = subvar,
      observed = observed,
      predicted = var1.pred,
      residual = observed - predicted
    ) %>%
    dplyr::filter(is.finite(observed), is.finite(predicted))
  
  if (nrow(cv_df) == 0) next
  
  cv_predictions_list[[subvar]] <- cv_df
  
  cv_summary_list[[subvar]] <- cv_df %>%
    dplyr::summarise(
      Cruise = Cruise,
      Transect_ID = Transect,
      variable = subvar,
      n = dplyr::n(),
      ME = mean(residual),
      MAE = mean(abs(residual)),
      RMSE = sqrt(mean(residual^2)),
      R2 = if (dplyr::n() > 1 && stats::sd(observed) > 0 &&
               stats::sd(predicted) > 0)
        stats::cor(observed, predicted)^2 else NA_real_,
      R2_pred = if (dplyr::n() > 1 &&
               sum((observed - mean(observed))^2) > 0)
        1 - sum(residual^2) / sum((observed - mean(observed))^2)
      else NA_real_,
      nmax = krige_nmax,
      nmin = krige_nmin,
      maxdist = krige_maxdist,
      anisotropy_m_per_km = anisotropy_m_per_km,
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Save cross-validation diagnostics
# ------------------------------------------------------------

cv_summary <- dplyr::bind_rows(cv_summary_list)
cv_predictions <- dplyr::bind_rows(cv_predictions_list)

readr::write_csv(
  cv_summary,
  file.path(output_dir_figs_normal_diagnostics, "kriging_cv_summary.csv")
)

readr::write_csv(
  cv_predictions,
  file.path(output_dir_figs_normal_diagnostics, "kriging_cv_predictions.csv")
)

# ------------------------------------------------------------
# Kriging validation figure
# ------------------------------------------------------------

if (nrow(cv_predictions) > 0) {
  
  p_cv <- ggplot2::ggplot(
    cv_predictions,
    ggplot2::aes(x = observed, y = predicted)
  ) +
    ggplot2::geom_point(alpha = 0.3, size = 0.8) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::labs(
      title = paste(Cruise, Transect, "ordinary kriging cross-validation"),
      x = "Observed",
      y = "Cross-validated prediction"
    ) +
    ggplot2::theme_classic()
  
  ggplot2::ggsave(
    filename = file.path(output_dir_figs_normal_diagnostics,
                         paste(Cruise, Transect, "kriging_cross_validation.png", sep = "-")),
    plot = p_cv, width = 12, height = 8, units = "in", dpi = 300
  )
}

message("Kriged")


scaled1 <- krige_biophys1

# ---------- 7) Save kriged, rescaled grids (includes Longitude) and assign relative distance for plotting ----------
scaled1 <- scaled1 %>%
  dplyr::rename(
    Cruise                  = Transect_ID,
    Transect_ID             = HD_ID,
    dist_kriged_km          = x,
    depth_kriged_m          = Depth,
    predicted_kriged_value  = var1.pred,
    krige_variance          = var1.var
  ) %>% mutate(Latitude = mean(dat$Latitude, na.rm = TRUE)) %>% assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude')

# ---------- 7.1)  ----------


scaled1 <- scaled1 %>%
  mutate(
    predicted_kriged_value = case_when(
      variable == "log10_integrated_biovol_all" ~ 10^predicted_kriged_value,
      variable == "log10_total_biovolume_um3_m3" ~ 10^predicted_kriged_value,
      variable == "log10_total_biovolume_sum" ~ 10^predicted_kriged_value,
      variable == "log10_integrated_biovol_sample" ~ 10^predicted_kriged_value,
      TRUE ~ predicted_kriged_value
    )
  )

scaled1 <- scaled1 %>%
  mutate(
    variable = case_when(
      variable == "log10_integrated_biovol_all" ~ "integrated_biovol_all",
      variable == "log10_total_biovolume_um3_m3" ~ "total_biovolume_um3_m3",
      variable == "log10_total_biovolume_sum"   ~ "total_biovolume_sum",
      variable == "log10_integrated_biovol_sample"   ~ "integrated_biovol_sample",
      TRUE ~ variable
    )
  )

# Save the final kriged predictions
readr::write_csv(
  scaled1,
  file.path(
    output_dir_data,
    paste(
      Cruise, Transect, size_var, instrument,
      paste0("size-spectrum_envs_", interp_tag, ".csv"),
      sep = "-"
    )
  )
)
# ---------- 8) seafloor polygon ----------
# load bathy data 
bathy_dir<- file.path(base_dir, 'data', 'bathymetry_polygons', paste(Cruise, Transect, 'bathymetry.csv', sep='-'))
bathy_df <- read_csv(bathy_dir, col_names = TRUE)

# Restrict polygon to kriged depth range
max_depth <- max(scaled1$depth_kriged_m, na.rm = TRUE)
min_x     <- min(scaled1$rel_dist, na.rm = TRUE)
max_x     <- max(scaled1$rel_dist, na.rm = TRUE)

lon_offset<- max(bathy_df$rel_dist) - max_x 

bathy_df$rel_dist<-bathy_df$rel_dist - lon_offset
# ---------- Build bathymetry polygon for plotting ----------

# Build the two edge rows with identical column structure
edge_top <- data.frame(
  rel_dist  = max(bathy_df$rel_dist, na.rm = TRUE),
  B_smooth  = max(bathy_df$B_smooth, na.rm = TRUE) + 10,
  Longitude = NA_real_
)

edge_bottom <- data.frame(
  rel_dist  = min(bathy_df$rel_dist, na.rm = TRUE),
  B_smooth  = max(bathy_df$B_smooth, na.rm = TRUE) + 10,
  Longitude = NA_real_
)

# Combine into a single polygon df
bathy_poly <- rbind(bathy_df, edge_top, edge_bottom)

  # keep only depths within plotted range
  # keep only rel_dist values inside slope grid window
  # order polygon points shoreward → offshore



# add "closing" vertices so polygon can fill
bathy_poly <- rbind(
  bathy_poly,
  data.frame(rel_dist = max(bathy_poly$rel_dist), B_smooth = max_depth, Longitude = NA),
  data.frame(rel_dist = min(bathy_poly$rel_dist), B_smooth = max_depth, Longitude = NA)
)

# Ensure offshore is smaller rel_dist, shore is larger rel_dist
bathy_df <- bathy_df %>%
  arrange(desc(rel_dist)) %>%  # start from shore → offshore
  
  # enforce monotonic deepening
  mutate(B_smooth_monotonic = cummax(B_smooth)) %>%
  
  # rename it for convenience
  rename(B_smooth_original = B_smooth) %>%
  rename(B_smooth = B_smooth_monotonic)

bathy_poly <- bathy_df %>%
  dplyr::select(rel_dist, B_smooth, Longitude) %>%
  dplyr::filter(!is.na(B_smooth), B_smooth <= max(scaled1$depth_kriged_m, na.rm = TRUE)) %>%
  arrange(rel_dist) %>%
  rbind(
    data.frame(
      rel_dist = max(.$rel_dist),
      B_smooth = max(scaled1$depth_kriged_m, na.rm = TRUE),
      Longitude = NA_real_
    ),
    data.frame(
      rel_dist = min(.$rel_dist),
      B_smooth = max(scaled1$depth_kriged_m, na.rm = TRUE),
      Longitude = NA_real_
    )
  )

# observation dataset for plotting
# ------------------------------------------------------------
# Original observations for kriging diagnostic overlays
# ------------------------------------------------------------

observed_points <- myDf %>%
  dplyr::mutate(
    Transect_ID = Transect
  ) %>%
  assign_relative_distance(
    land_ref_path = shore_path,
    lat_col = "Latitude",
    long_col = "Longitude"
  ) %>%
  dplyr::filter(
    is.finite(rel_dist),
    is.finite(Depth)
  )



# ---------- 9) Plots (same behavior as original) ----------
# base (filled) data for slope
fmt <- scales::label_scientific(digits = 2)

s_df <- scaled1 %>% dplyr::filter(variable == "slope") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

a_df <- scaled1 %>% filter(variable == "slope_anomaly") %>% drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

i_df <- scaled1 %>% dplyr::filter(variable == "intercept") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

bts_df <- scaled1 %>% dplyr::filter(variable == "total_biovolume_sum") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

bt_df <- scaled1 %>% dplyr::filter(variable == "total_biovolume_um3_m3") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

bia_df <- scaled1 %>% dplyr::filter(variable == "integrated_biovol_all") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

bis_df <- scaled1 %>% dplyr::filter(variable == "integrated_biovol_sample") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

bd_df <- scaled1 %>% dplyr::filter(variable == "biovol_diff") %>% tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

# overlay data: Temperature + Chl-a
temp_df <- scaled1 %>%
  dplyr::filter(variable == "Temperature") %>%
  transmute(rel_dist, depth_kriged_m, temp = predicted_kriged_value) %>%
  tidyr::drop_na()

chla_df <- scaled1 %>%
  dplyr::filter(variable == "chl_a_ug_l") %>%
  transmute(rel_dist, depth_kriged_m, chla = predicted_kriged_value) %>%
  tidyr::drop_na()

# nice breaks for the contour lines (adjust n as you like)
brk_temp <- pretty(range(temp_df$temp, na.rm = TRUE), n = 6)
brk_chla <- pretty(range(chla_df$chla, na.rm = TRUE), n = 6)

# ----- choose breaks (robust: fallback if range is tiny) -----
rng_t <- range(temp_df$temp, na.rm = TRUE)
rng_c <- range(chla_df$chla, na.rm = TRUE)

brk_temp <- pretty(rng_t, n = 7); if (length(unique(brk_temp)) < 2) brk_temp <- seq(rng_t[1], rng_t[2], length.out = 7)
brk_chla <- pretty(rng_c, n = 7); if (length(unique(brk_chla)) < 2) brk_chla <- seq(rng_c[1], rng_c[2], length.out = 7)

if (instrument=='ISIIS'){
   slope_breaks <- c(-Inf, seq(-1.8, 0, by = 0.25), Inf)
} else if (instrument =='LISST-ISIIS'){
  slope_breaks <- c(-Inf, seq(-1.1, -0.2, by = 0.1), Inf)
} else if (instrument == 'LISST'){
  slope_breaks <- c(-Inf, seq(-1.1, -0.2, by = 0.1), Inf)
}

p_slope <- ggplot(s_df, aes(x = rel_dist, y = depth_kriged_m)) +
  # filled contours of SLOPE
  geom_contour_filled(aes(z = predicted_kriged_value), breaks = slope_breaks, na.rm=TRUE) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  # overlay: TEMPERATURE contour lines
  geom_contour(
    data = temp_df,
    aes(z = temp, color = "Temperature"),
    bins = 8, linewidth = 0.4
  ) +
  # overlay: CHL-A contour lines (dotted)
  geom_contour(
    data = chla_df,
    aes(z = chla, color = "Chl-a"),
    bins = 8, linewidth = 0.4, linetype = "dotted"
  ) +
  # axes tight to data
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0), limits = c(max(s_df$depth_kriged_m, na.rm = TRUE), 0)) +
  
  scale_fill_manual(
    drop = FALSE,
    values = viridis::viridis(length(slope_breaks) - 1),
    labels = paste(
      head(slope_breaks, -1),
      slope_breaks[-1],
      sep = " – "
    ),
    name = expression(paste("Slope (", mu, "m"^{-3}, " ", m^{-3}, ")"))
  )+

  # line legend for overlays
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
    labels = c(
      "Temperature" = expression(paste("Temperature (", degree*C, ")")),
      "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
    )
  ) +
  labs(
    title = paste0("Slope (+ Temp & Chl-a contours) — ", unique(s_df$Cruise), " / ", unique(s_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text()
  )


p_slope<- p_slope +
  # temperature isolines
  geom_contour(
    data = temp_df,
    aes(z = temp, color = "Temperature"),
    breaks = brk_temp, linewidth = 0.4
  ) +
  # temperature labels
  metR::geom_text_contour(
    data = temp_df,
    aes(z = temp, label = after_stat(level), color = "Temperature"),
    breaks = brk_temp,
    size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  # chl-a isolines
  geom_contour(
    data = chla_df,
    aes(z = chla, color = "Chl-a"),
    breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
  ) +
  # chl-a labels (formatted)
  metR::geom_text_contour(
    data = chla_df,
    aes(
      z = chla,
      label = scales::number(after_stat(level), accuracy = 0.1),
      color = "Chl-a"
    ),
    breaks = brk_chla,
    size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  )



# semi-transparent white polygon overlay for extrapolated regions

p_slope <- p_slope +
  geom_tile(
    data = mask_slope,
    aes(x = rel_dist, y = Depth),
    fill = "white", alpha = 0.5
  )

p_slope <- p_slope+
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, planktonFraction, instrument, "slope_interp.png", sep = "-")),
       plot = p_slope, width  = 12, height = 6, units  = "in", dpi = 300)

#slope anomaly
an_max <- max(abs(range(s_df$predicted_kriged_value, na.rm = TRUE)))
an_breaks <- seq(-an_max, an_max, length.out = 9)

p_anomaly <- ggplot(s_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_raster(aes(fill = predicted_kriged_value)) +
  
  geom_contour(aes(z = predicted_kriged_value),
               color = "black", alpha = 0.3, bins = 5) +
  metR::geom_text_contour(
    aes(z = predicted_kriged_value,
        label = scales::number(after_stat(level), accuracy = 0.1)),
    stroke = 0.2,
    size = 5,
    check_overlap = TRUE,
    skip = 0.5,    # reduce label density (adjust as needed)
    min.size = 20 
  )+
  
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", linewidth = 2) +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  colorspace::scale_fill_continuous_divergingx(
    palette = "RdBu",
    mid = -1,
    rev = TRUE,
    limits = c(-2, 0),
    name = "Slope"
  ) +
  
  labs(
    title = paste0("Slope (relative to -1) — ",
                   unique(a_df$Cruise), " / ", unique(a_df$Transect_ID)),
    x = "Relative distance",
    y = "Depth (m)"
  ) +
  
  theme_classic(base_size = 12)



p_anomaly <- p_anomaly +
  geom_tile(data = mask_slope,
            aes(x = rel_dist, y = Depth),
            fill = "white", alpha = 0.5) +
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

p_anomaly <- p_anomaly +
  theme(
    # Axis titles
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    
    # Tick labels (numbers)
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    
    # Tick marks themselves
    axis.ticks.length = unit(6, "pt"),
    axis.ticks = element_line(linewidth = 3)
  )

ggsave(
  filename = file.path(output_dir_figs_normal,
                       paste(Cruise, Transect, size_var, planktonFraction, instrument, "slopeANOM_interp.png", sep = "-")),
  plot = p_anomaly, width = 12, height = 6, units = "in", dpi = 300
)


p_anomaly_obs <- p_anomaly +
  ggplot2::geom_point(
    data = observed_points %>%
      dplyr::filter(is.finite(slope)),
    mapping = ggplot2::aes(
      x = rel_dist,
      y = Depth
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 0.4,
    alpha = 0.35
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir_figs_normal_diagnostics,
    paste(Cruise, Transect, "slope_with_observations.png", sep = "-")
  ),
  plot = p_anomaly_obs,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)

#intercept plot
if (instrument=='ISIIS'){
  intercept_breaks <- c(-Inf, seq(5, 20, by = 1.6), Inf) 
} else if (instrument =='LISST-ISIIS'){
  intercept_breaks <- c(-Inf, seq(6.5, 11.5, by = 1), Inf)
} else if (instrument == 'LISST'){
  intercept_breaks <- c(-Inf, seq(6.5, 11.5, by = 1), Inf)
}


p_intercept <- ggplot(i_df, aes(x = rel_dist, y = depth_kriged_m))  +
  # filled contours for INTERCEPT
  geom_contour_filled(aes(z = predicted_kriged_value), breaks = intercept_breaks, na.rm = TRUE) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  # --- overlays: Temperature and Chl-a isolines + labels (same as p_slope) ---
  geom_contour(
    data = temp_df, aes(z = temp, color = "Temperature"),
    breaks = brk_temp, linewidth = 0.4
  ) +
  metR::geom_text_contour(
    data = temp_df,
    aes(z = temp, label = after_stat(level), color = "Temperature"),
    breaks = brk_temp, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  geom_contour(
    data = chla_df, aes(z = chla, color = "Chl-a"),
    breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
  ) +
  metR::geom_text_contour(
    data = chla_df,
    aes(z = chla, label = scales::number(after_stat(level), accuracy = 0.1), color = "Chl-a"),
    breaks = brk_chla, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  # axes tight
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  scale_fill_manual(
    drop = FALSE,
    values = viridis::viridis(length(intercept_breaks) - 1),
    labels = paste(
      head(intercept_breaks, -1),
      intercept_breaks[-1],
      sep = " – "
    ),
    name = expression(paste("Intercept (", mu, "m"^3, " ", mu, "m"^{-3}, " ", m^{-3}, ")"))
  )+

  # line legend for overlays (same as p_slope)
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
    labels = c(
      "Temperature" = expression(paste("Temperature (", degree*C, ")")),
      "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
    )
  ) +
  labs(
    title = paste0("Intercept (+ Temp & Chl-a contours) — ",
                   unique(i_df$Cruise), " / ", unique(i_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")


# semi-transparent white polygon overlay for extrapolated regions

p_intercept <- p_intercept +
  geom_tile(
    data = mask_slope,
    aes(x = rel_dist, y = Depth),
    fill = "white", alpha = 0.5
  )

p_intercept <- p_intercept + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, planktonFraction, "intercept_interp0.png", sep = "-")),
       plot = p_intercept, width  = 12, height = 6, units  = "in", dpi = 300)

p_intercept_obs <- p_intercept +
  ggplot2::geom_point(
    data = observed_points %>%
      dplyr::filter(is.finite(intercept)),
    mapping = ggplot2::aes(
      x = rel_dist,
      y = Depth
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 0.4,
    alpha = 0.35
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir_figs_normal_diagnostics,
    paste(Cruise, Transect, "intercept_with_observations.png", sep = "-")
  ),
  plot = p_intercept_obs,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)


#biovol plots
if (instrument=='ISIIS'){
  biovol_breaks <- c(-Inf, 10^(9:20), Inf)
} else if (instrument =='LISST-ISIIS'){
  biovol_breaks <- c(-Inf, 10^(seq(11, 15, by = 0.3)), Inf)
} else if (instrument == 'LISST'){
  biovol_breaks <- c(-Inf, 10^(9:20), Inf) 
}

#total biovolume
p_biovol_ts <- ggplot(bts_df, aes(x = rel_dist, y = depth_kriged_m))  +
  # filled contours for INTERCEPT
  geom_contour_filled(aes(z = predicted_kriged_value), breaks = biovol_breaks,  na.rm = TRUE) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  # --- overlays: Temperature and Chl-a isolines + labels (same as p_slope) ---
  geom_contour(
    data = temp_df, aes(z = temp, color = "Temperature"),
    breaks = brk_temp, linewidth = 0.4
  ) +
  metR::geom_text_contour(
    data = temp_df,
    aes(z = temp, label = after_stat(level), color = "Temperature"),
    breaks = brk_temp, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  geom_contour(
    data = chla_df, aes(z = chla, color = "Chl-a"),
    breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
  ) +
  metR::geom_text_contour(
    data = chla_df,
    aes(z = chla, label = scales::number(after_stat(level), accuracy = 0.1), color = "Chl-a"),
    breaks = brk_chla, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  # axes tight
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  scale_fill_manual(
    drop = FALSE,
    values = viridis::viridis(length(biovol_breaks) - 1),
    labels = paste(fmt(head(biovol_breaks, -1)),
                   fmt(biovol_breaks[-1]),
                   sep = " – "),
    name = expression(paste("Biovolume (", mu, "m"^3, " ", m^{-3}, ")"))
  )+
  
  # line legend for overlays (same as p_slope)
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
    labels = c(
      "Temperature" = expression(paste("Temperature (", degree*C, ")")),
      "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
    )
  ) +
  labs(
    title = paste0(" Biovolume (+ Temp & Chl-a contours) — ",
                   unique(i_df$Cruise), " / ", unique(i_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")


# semi-transparent white polygon overlay for extrapolated regions

p_biovol_ts <- p_biovol_ts +
  geom_tile(
    data = mask_slope,
    aes(x = rel_dist, y = Depth),
    fill = "white", alpha = 0.5
  )

p_biovol_ts <- p_biovol_ts + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, planktonFraction, "biovol_totalSum_interp0.png", sep = "-")),
       plot = p_biovol_ts, width  = 12, height = 6, units  = "in", dpi = 300)



#total biovolume per cubic meter
p_biovol_total <- ggplot(bt_df, aes(x = rel_dist, y = depth_kriged_m))  +
  # filled contours for INTERCEPT
  geom_contour_filled(aes(z = predicted_kriged_value), breaks = biovol_breaks,  na.rm = TRUE) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  # --- overlays: Temperature and Chl-a isolines + labels (same as p_slope) ---
  geom_contour(
    data = temp_df, aes(z = temp, color = "Temperature"),
    breaks = brk_temp, linewidth = 0.4
  ) +
  metR::geom_text_contour(
    data = temp_df,
    aes(z = temp, label = after_stat(level), color = "Temperature"),
    breaks = brk_temp, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  geom_contour(
    data = chla_df, aes(z = chla, color = "Chl-a"),
    breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
  ) +
  metR::geom_text_contour(
    data = chla_df,
    aes(z = chla, label = scales::number(after_stat(level), accuracy = 0.1), color = "Chl-a"),
    breaks = brk_chla, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  # axes tight
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  scale_fill_manual(
    drop = FALSE,
    values = viridis::viridis(length(biovol_breaks) - 1),
    labels = paste(fmt(head(biovol_breaks, -1)),
                   fmt(biovol_breaks[-1]),
                   sep = " – "),
    name = expression(paste("Biovolume (", mu, "m"^3, " ", m^{-3}, ")"))
  )+
  
  # line legend for overlays (same as p_slope)
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
    labels = c(
      "Temperature" = expression(paste("Temperature (", degree*C, ")")),
      "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
    )
  ) +
  labs(
    title = paste0(" Total Biovolume (+ Temp & Chl-a contours) — ",
                   unique(bt_df$Cruise), " / ", unique(bt_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")


# semi-transparent white polygon overlay for extrapolated regions

p_biovol_total <- p_biovol_total +
  geom_tile(
    data = mask_slope,
    aes(x = rel_dist, y = Depth),
    fill = "white", alpha = 0.5
  )

p_biovol_total <- p_biovol_total + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, planktonFraction, "biovol_total_um3_m3_interp0.png", sep = "-")),
       plot = p_biovol_total, width  = 12, height = 6, units  = "in", dpi = 300)






# intregrated biovolume for full size range
p_biovol_ia <- ggplot(bia_df, aes(x = rel_dist, y = depth_kriged_m))  +
  # filled contours for INTERCEPT
  geom_contour_filled(aes(z = predicted_kriged_value), breaks = biovol_breaks,  na.rm = TRUE) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  # --- overlays: Temperature and Chl-a isolines + labels (same as p_slope) ---
  geom_contour(
    data = temp_df, aes(z = temp, color = "Temperature"),
    breaks = brk_temp, linewidth = 0.4
  ) +
  metR::geom_text_contour(
    data = temp_df,
    aes(z = temp, label = after_stat(level), color = "Temperature"),
    breaks = brk_temp, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  geom_contour(
    data = chla_df, aes(z = chla, color = "Chl-a"),
    breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
  ) +
  metR::geom_text_contour(
    data = chla_df,
    aes(z = chla, label = scales::number(after_stat(level), accuracy = 0.1), color = "Chl-a"),
    breaks = brk_chla, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  # axes tight
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  scale_fill_manual(
    drop = FALSE,
    values = viridis::viridis(length(biovol_breaks) - 1),
    labels = paste(fmt(head(biovol_breaks, -1)),
                   fmt(biovol_breaks[-1]),
                   sep = " – "),
    name = expression(paste("Biovolume (", mu, "m"^3, " ", m^{-3}, ")"))
  )+
  
  # line legend for overlays (same as p_slope)
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
    labels = c(
      "Temperature" = expression(paste("Temperature (", degree*C, ")")),
      "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
    )
  ) +
  labs(
    title = paste0("int. Biovolume all sizes (+ Temp & Chl-a contours) — ",
                   unique(i_df$Cruise), " / ", unique(i_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")


# semi-transparent white polygon overlay for extrapolated regions

p_biovol_ia <- p_biovol_ia +
  geom_tile(
    data = mask_slope,
    aes(x = rel_dist, y = Depth),
    fill = "white", alpha = 0.5
  )

p_biovol_ia <- p_biovol_ia + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, planktonFraction, "biovol_integratedAll_interp0.png", sep = "-")),
       plot = p_biovol_ia, width  = 12, height = 6, units  = "in", dpi = 300)

#integrated biovolume per sample:

p_biovol_is <- ggplot(bis_df, aes(x = rel_dist, y = depth_kriged_m))  +
  # filled contours for INTERCEPT
  geom_contour_filled(aes(z = predicted_kriged_value), breaks = biovol_breaks,  na.rm = TRUE) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  # --- overlays: Temperature and Chl-a isolines + labels (same as p_slope) ---
  geom_contour(
    data = temp_df, aes(z = temp, color = "Temperature"),
    breaks = brk_temp, linewidth = 0.4
  ) +
  metR::geom_text_contour(
    data = temp_df,
    aes(z = temp, label = after_stat(level), color = "Temperature"),
    breaks = brk_temp, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  geom_contour(
    data = chla_df, aes(z = chla, color = "Chl-a"),
    breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
  ) +
  metR::geom_text_contour(
    data = chla_df,
    aes(z = chla, label = scales::number(after_stat(level), accuracy = 0.1), color = "Chl-a"),
    breaks = brk_chla, size = 3, stroke = 0.2, check_overlap = TRUE, show.legend = FALSE
  ) +
  # axes tight
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  scale_fill_manual(
    drop = FALSE,
    values = viridis::viridis(length(biovol_breaks) - 1),
    labels = paste(fmt(head(biovol_breaks, -1)),
                   fmt(biovol_breaks[-1]),
                   sep = " – "),
    name = expression(paste("Biovolume (", mu, "m"^3, " ", m^{-3}, ")"))
  )+
  
  # line legend for overlays (same as p_slope)
  scale_color_manual(
    name = NULL,
    values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
    labels = c(
      "Temperature" = expression(paste("Temperature (", degree*C, ")")),
      "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
    )
  ) +
  labs(
    title = paste0(" int. Biovolume per sample (+ Temp & Chl-a contours) — ",
                   unique(i_df$Cruise), " / ", unique(i_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")


# semi-transparent white polygon overlay for extrapolated regions

p_biovol_is <- p_biovol_is +
  geom_tile(
    data = mask_slope,
    aes(x = rel_dist, y = Depth),
    fill = "white", alpha = 0.5
  )

p_biovol_is <- p_biovol_is + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, planktonFraction, "biovol_integratedSample_interp0.png", sep = "-")),
       plot = p_biovol_is, width  = 12, height = 6, units  = "in", dpi = 300)

#percentage diff between total biovolume sampled and integrated biovolume per sample

p_biovol_df <- ggplot(bd_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_raster(aes(fill = predicted_kriged_value)) +
  
  geom_contour(aes(z = predicted_kriged_value),
               color = "black", alpha = 0.3, bins = 5) +
  metR::geom_text_contour(
    aes(z = predicted_kriged_value,
        label = scales::number(after_stat(level), accuracy = 0.1)),
    stroke = 0.2,
    size = 5,
    check_overlap = TRUE,
    skip = 0.5,    # reduce label density (adjust as needed)
    min.size = 20 
  )+
  
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "black", linewidth = 2) +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  
  colorspace::scale_fill_continuous_divergingx(
    palette = "RdBu",
    mid = 0,
    rev = TRUE,
    limits = c(-150, 150),
    name = "Slope"
  ) +
  
  labs(
    title = paste0("percentage diff. biovol — ",
                   unique(a_df$Cruise), " / ", unique(a_df$Transect_ID)),
    x = "Relative distance",
    y = "Depth (m)"
  ) +
  
  theme_classic(base_size = 12)



p_biovol_df <- p_biovol_df +
  geom_tile(data = mask_slope,
            aes(x = rel_dist, y = Depth),
            fill = "white", alpha = 0.5) +
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )

p_biovol_df <-p_biovol_df +
  theme(
    # Axis titles
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    
    # Tick labels (numbers)
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    
    # Tick marks themselves
    axis.ticks.length = unit(6, "pt"),
    axis.ticks = element_line(linewidth = 3)
  )

ggsave(
  filename = file.path(output_dir_figs_normal,
                       paste(Cruise, Transect, size_var, planktonFraction, instrument, "Biovol_diffSample_interp.png", sep = "-")),
  plot = p_biovol_df, width = 12, height = 6, units = "in", dpi = 300
)




chl_df <- scaled1 %>%
  dplyr::filter(variable == "chl_a_ug_l") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

p_chla <- ggplot(chl_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_viridis_d(
    name = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
  ) +
  labs(
    title = paste0("Chlorophyll-a — ", unique(chl_df$Cruise), " / ", unique(chl_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_chla <- p_chla + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "chla_interp.png", sep = "-")),
       plot = p_chla, width  = 12, height = 6, units  = "in", dpi = 300)

t_df <- scaled1 %>%
  dplyr::filter(variable == "Temperature") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

p_temp <- ggplot(t_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  # Blue (cold) -> Red (warm)
  scale_fill_brewer(
    palette = "RdYlBu", direction = -1,
    name = expression(paste("Temperature (", degree*C, ")"))
  ) +
  labs(
    title = paste0("Temperature — ", unique(t_df$Cruise), " / ", unique(t_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_temp <- p_temp + 
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "temperature_interp.png", sep = "-")),
       plot = p_temp, width  = 12, height = 6, units  = "in", dpi = 300)

sali_df <- scaled1 %>%
  dplyr::filter(variable == "Salinity") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

p_sali <- ggplot(sali_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_brewer(
    palette = "RdYlBu", direction = -1,
    name = "Salinity (PSU)"
  ) +
  labs(
    title = paste0("Salinity — ", unique(sali_df$Cruise), " / ", unique(sali_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_sali+
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "salinity_interp.png", sep = "-")),
       plot = p_sali, width  = 12, height = 6, units  = "in", dpi = 300)

o2_df <- scaled1 %>%
  dplyr::filter(variable == "Oxygen_ml_l") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

p_o2 <- ggplot(o2_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_brewer(
    palette = "RdYlBu", direction = -1,
    name = expression(paste("Oxygen (mL ", L^{-1}, ")"))
  ) +
  labs(
    title = paste0("Oxygen — ", unique(o2_df$Cruise), " / ", unique(o2_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_o2+
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "oxygen_interp.png", sep = "-")),
       plot = p_o2, width  = 12, height = 6, units  = "in", dpi = 300)

par_df <- scaled1 %>%
  dplyr::filter(variable == "PAR") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

p_par <- ggplot(par_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_brewer(
    palette = "RdYlBu", direction = -1,
    name = "PAR (unitless)"
  ) +
  labs(
    title = paste0("PAR — ", unique(par_df$Cruise), " / ", unique(par_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_par+
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "PAR_interp.png", sep = "-")),
       plot = p_par, width  = 12, height = 6, units  = "in", dpi = 300)

ph_df <- scaled1 %>%
  dplyr::filter(variable == "pH") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

p_ph <- ggplot(ph_df, aes(x = rel_dist, y = depth_kriged_m)) +
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_brewer(
    palette = "RdYlBu", direction = -1,
    name = "pH"
  ) +
  labs(
    title = paste0("pH — ", unique(ph_df$Cruise), " / ", unique(ph_df$Transect_ID)),
    x = "Relative distance (dimensionless)",
    y = "Depth (m)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_ph+
  # add bathymetry
  geom_polygon(
    data = bathy_poly,
    aes(x = rel_dist, y = B_smooth),
    fill = "gray40", alpha = 1
  )
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "ph_interp.png", sep = "-")),
       plot = p_ph, width  = 12, height = 6, units  = "in", dpi = 300)



