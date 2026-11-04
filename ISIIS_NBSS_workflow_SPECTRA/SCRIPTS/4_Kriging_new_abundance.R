setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
library(ggplot2)
library(colorRamps)
library(plyr)
library(dplyr)
library(reshape2)
library(akima)
library(scales)
library(kriging)
library(ggthemes)
library(forcats)
library(optparse)
library(readr)
base_dir<- getwd()

# --- 0) Paths & identifiers (ADAPT: uses your base_dir and w23_NH) ----
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

option_list <- list(
  make_option(c("-c", "--Cruise"), type = "character", help = "Cruise: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t", "--Transect"), type = "character", help = "Cruise: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)


#variables to load data:
abundanceDir <- file.path(base_dir, 'data', 'binned_abundance')
message("abundance_dir: ", abundanceDir)
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "ancillary", "transect_shelf_land.csv")

Cruise <- opt$Cruise # Cruise<- 'w23'
Transect <- opt$Transect # Transect<- 'NH'


# checking if directories extist and create them if not
output_dir_data <- file.path(base_dir, "data", "interp_and_krige_dist10", 'abundance_data')
if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_data)
}

output_dir_figs_normal <- file.path(base_dir, "figures", "interp_and_krige_dist10", 'abundance', 'normal', Cruise, Transect)
if (!dir.exists(output_dir_figs_normal)) {
  dir.create(output_dir_figs_normal, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_normal)
}

output_dir_figs_log <- file.path(base_dir, "figures", "interp_and_krige_dist10", 'abundance', 'log', Cruise, Transect)
if (!dir.exists(output_dir_figs_log)) {
  dir.create(output_dir_figs_log, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_log)
}


# Build a regex pattern: e.g., "^NH_w23_binned_abundance_test_chunk_.*\\.csv$" and list binned abundance files
pat <- sprintf("^%s_%s_binned_abundance_chunk_.*\\.csv$", Transect, Cruise)

df_list <- list.files(abundanceDir, pattern = pat, full.names = TRUE, recursive  = TRUE)

print(df_list)
#create list of env variables with the mean suffix for indexing
mean_env_vars <- paste0("mean_", env_vars)


#load and clean up data
dat <- map_dfr(df_list, ~ readr::read_csv( .x, show_col_types = FALSE,
              col_select = c( mean_time_Pacific, depth_mid, mean_Latitude, mean_Longitude,
      dplyr::any_of(mean_env_vars),
      dplyr::contains("ind_m3", ignore.case = TRUE)
    )
  ) |> dplyr::distinct()
) |>
  dplyr::rename_with(~ sub("^mean_", "", .x), .cols = dplyr::starts_with("mean_")) |>
  dplyr::arrange(time_Pacific) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    ind_m3_total = sum(dplyr::c_across(dplyr::contains("ind_m3") & !dplyr::matches("^ind_m3_total$")),
                       na.rm = TRUE)
  ) |>
  dplyr::ungroup()

dat$Cruise<-Cruise
dat$Transect<-Transect


binned_plot<- ggplot(dat, aes(x = Longitude, y = depth_mid)) +
  geom_point(aes(color = ind_m3_total), alpha = 0.6, size = 1) +
  scale_y_reverse() +
  scale_color_viridis_c(option = "viridis") +
  labs(x = "Longitude", y = "Depth (m)", color = "ind_m3_total") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect,  "binned_data_plot.png", sep = "-")),
       plot = binned_plot, width  = 12, height = 6, units  = "in", dpi = 300)

# to get rough resolution in Km:


# ---------- 1) Build myDf from your `dat` (use Longitude -> distance km) ----------
# Expect: dat has at least: time_Pacific, rel_dist, depth_mid, Latitude, Longitude, slope

# We'll compute distance in km from the westernmost longitude (start of transect):
# This preserves the original pipeline's functionality (distance in km along x).


# second round of standardization and cleanup for kriging
# Use dplyr:: to avoid plyr clashes
myDf <- dat %>%
  dplyr::mutate(
    Transect_ID = Cruise,
    HD_ID       = Transect,
    Depth       = as.numeric(depth_mid),   # vertical axis (m)
    Longitude   = as.numeric(Longitude),
    Latitude    = as.numeric(Latitude)
  ) %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, Latitude, Longitude,
    dplyr::any_of(env_vars),              # env vars (already de-prefixed earlier)
    time_Pacific,
    dplyr::contains("ind_m3", ignore.case = TRUE)  # all abundance columns
  )

# Compute km-per-degree at mean latitude, and define the start longitude (easternmost) 
lat0 <- mean(myDf$Latitude, na.rm = TRUE)
km_per_deg_lon <- 111.320 * cos(lat0 * pi/180)

lon0 <- max(myDf$Longitude, na.rm = TRUE)   # start = easternmost point

# Distance grows westward (more negative longitudes):
#   east -> west : GPS.dist.from.start.km increases
myDf <- myDf %>%
  dplyr::mutate(GPS.dist.from.start.km = (lon0 - Longitude) * km_per_deg_lon)

# Add along-transect distance (km) from start longitude; keep Longitude in the table (this chunk might be redundant) 
myDf <- myDf %>%
  dplyr::mutate(
    GPS.dist.from.start.km = (lon0 - Longitude) * km_per_deg_lon
  ) %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, GPS.dist.from.start.km,
    dplyr::any_of(env_vars),          # env vars (already de-prefixed earlier)
    Latitude, Longitude, time_Pacific,
    dplyr::contains("ind_m3", ignore.case = TRUE)  # all abundance columns
  )

# Drop odd columns if they exist (same as original script) ( also redundant, columns that we need were already selected)
if ("NA" %in% names(myDf)) myDf$`NA` <- NULL
if ("NA_mean_area"  %in% names(myDf)) myDf$NA_mean_area  <- NULL
if ("NA_mean_major" %in% names(myDf)) myDf$NA_mean_major <- NULL
if ("NA_mean_minor" %in% names(myDf)) myDf$NA_mean_minor <- NULL
gc()


# ---------- 2) Long variable table (variables to krige; here: abundances ind_m3) ----------

abund_cols <- names(myDf)[grepl("ind_m3", names(myDf), ignore.case = TRUE)]

#this is a dataset with only abundance data and spatial variables
biophys.10 <- myDf %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, GPS.dist.from.start.km,
    dplyr::all_of(abund_cols)
  ) %>%
  # make sure everything we will melt is numeric
  dplyr::mutate(
    dplyr::across(
      -c(Transect_ID, HD_ID, Depth, GPS.dist.from.start.km),
      ~ suppressWarnings(as.numeric(.))
    )
  )

message("selected columns for krige: all *ind_m3*")

# ---------- long table ----------

# long versiion of biophys10 
biophys.11 <- reshape2::melt(
  biophys.10,
  id.vars = c("Transect_ID","HD_ID","Depth","GPS.dist.from.start.km")
)

# ranks/ntiles per variable
biophys.11 <- biophys.11 %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(
    PCT   = percent_rank(value),
    ntile = ntile(value, 100)
  )

# ---------- cap extreme upper tail per variable & HD_ID ----------
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
message("reset outliers for  all *ind_m3* variables")



# ---------- 3) Environmental data (optional; same behavior) ----------
env_keep <- env_vars[env_vars %in% names(myDf)]

# similar to biophys 10 and 11, wide to long for env variables.
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


# ---------- 4) Coarse interpolation (akima) per variable (xo 0.5 km, yo 5 m) ---------- This and step 5 were generating strange patterns
# 
# 
# 
#   
#   
#   
#   
# 
# 
# 


 # ---------- 5) Kriging (fine grid; Exp, maxdist=10) ---------- # old kriging
 krige_biophys1 <- NULL
 variable2 <- unique(biophys.14$variable) # variable2 <- c('crab_zoea_corrected', ...)
 
 
 for (subvar in variable2) {
   x <- biophys.14[biophys.14$variable == subvar, ] %>% stats::na.omit()
   if (nrow(x) == 0) next
   
   print(unique(x$variable))
   
   x0 <- seq(min(x$GPS.dist.from.start.km, na.rm = TRUE),
             max(x$GPS.dist.from.start.km, na.rm = TRUE), by = 0.5) # edit distance grid here, OG 0.5
   y0 <- seq(0, max(x$Depth, na.rm = TRUE), by = 2) #edit depth grid here, OG 2
   
   d0 <- expand.grid(x = x0, y = y0)
   sp::coordinates(d0) <- ~ x + y
   
   sp::coordinates(x) <- ~ GPS.dist.from.start.km + Depth
   
   x.var <- gstat::variogram(value ~ 1, x)
   possibleError <- tryCatch(gstat::fit.variogram(x.var, model = gstat::vgm("Exp")),
                             error = function(e) e)
   if (inherits(possibleError, "error")) next
   
   x_fit <- gstat::fit.variogram(x.var, model = gstat::vgm("Exp"))
   
   possibleError2 <- tryCatch(gstat::krige(value ~ 1, x, d0, model = x_fit, maxdist = 10),
                              error = function(e) e)
   if (inherits(possibleError2, "error")) next
   
   x_kriged <- gstat::krige(value ~ 1, x, d0, model = x_fit, maxdist = 10)
   x_kriged_df <- as.data.frame(x_kriged)
   
   # --- Important: add Longitude back to the grid using the inverse transform ---
   # x_kriged_df$x is distance (km); convert to degrees and add lon0
   x_kriged_df$Longitude <- lon0 - (x_kriged_df$x / km_per_deg_lon)
   
   x_kriged_df$variable    <- subvar
   x_kriged_df$HD_ID       <- Transect
   x_kriged_df$Transect_ID <- Cruise
   
   krige_biophys1 <- rbind(krige_biophys1, x_kriged_df)
 }
 
 message("Kriged")

# # ---------- 6) Rescale predictions back to original [min, max] per variable ----------
 scaled1 <- NULL
 scale.variab <- unique(krige_biophys1$variable)
 
 for (subvar in scale.variab) {
   di <- plyr::ddply(
     .data = krige_biophys1[krige_biophys1$variable == subvar, ],
     .variables = ~ HD_ID,
     function(x) {
       values.select1 <- biophys.14$value[biophys.14$variable == subvar] %>% stats::na.omit()
       minv <- min(values.select1); maxv <- max(values.select1)
       x$var1.pred <- scales::rescale(x$var1.pred, to = c(minv, maxv))
       x
     }
   )
   scaled1 <- rbind(scaled1, di)
 }
 message("scaled back to original value range")
 gc()
 
 # Also keep Longitude in the final, rescaled set (present from step 5)
 # (If for any reason it's missing, recompute it here)
 if (!"Longitude" %in% names(scaled1)) {
   if ("dist_kriged_km" %in% names(scaled1)) {
     scaled1$Longitude <- lon0 - scaled1$dist_kriged_km / km_per_deg_lon
   } else if ("x" %in% names(scaled1)) {
     scaled1$Longitude <- lon0 - scaled1$x / km_per_deg_lon
   }
 }
 
 # ---------- 7) Save kriged, rescaled grids (includes Longitude) and assign relative distance for plotting ----------
 scaled1 <- scaled1 %>%
   dplyr::rename(
     Cruise                  = Transect_ID,
     Transect_ID             = HD_ID,
     dist_kriged_km          = x,
     depth_kriged_m          = y,
     predicted_kriged_value  = var1.pred,
     krige_variance          = var1.var
  ) %>% mutate(Latitude = mean(dat$Latitude, na.rm = TRUE)) %>% assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude')
 
 write.csv(scaled1, file.path(output_dir_data, paste(Cruise, Transect, 'plankton-abundance_envs_interp_and_krige_dist10.csv', sep='-')), row.names=FALSE)


# # ---------- 5) Kriging (fine grid; Spherical, block kriging, log-handling for abund) ----------
# 
# # Grid & neighborhood controls (tweak as desired)
# 
# 
#   # data for this variable
#   
#   
#   # (a) LOG-HANDLING for abundance variables (ind_m3*)
#     # small epsilon per variable to handle zeros
#   
#   # OPTIONAL: thin duplicates within small cells to reduce speckle
#   
#   # (b) RESCALE DEPTH to km for kriging geometry
#   
#   # Prediction grid in (km, km)
#   
#   
#   # set coordinates
#   
#   # Variogram and fit (on transformed scale if abundance)
#   
#   # (c) BLOCK KRIGING: block size in the same units (km, km)
#   
#   
#   
#   # (a, cont.) BACK-TRANSFORM abundance predictions with Laurent correction
#     # mean on natural scale for lognormal:
#     #   E[Y] = exp(mu + 0.5 * sigma^2) - eps
#     # (var1.var remains variance on log scale; keep if useful diagnostically)
#   
#   # Convert y back to meters for downstream plotting convenience
#   
#   # (d) LONGITUDE back-transform from km to degrees (east->west increasing)
#   
#   # Tag metadata
#   
# 
# 
# # ---------- 6) Prepare final kriged grid (NO rescale) ----------
# # We now trust var1.pred to already be on the correct scale:
# #  - abundance variables are back-transformed from log with Laurent correction
# #  - environmental variables were kriged on their native scale
# 
# # Ensure Longitude exists (Step 5 should have created it, but guard just in case)
#   # x is along-transect distance in km measured westward from lon0
# 
# # ---------- 7) Rename columns for plotting + add rel_dist; then save ----------
#   # constant latitude (needed by assign_relative_distance); use mean of input data
#   # compute rel_dist using your shoreline reference file
# 
# # Write CSV with everything needed downstream
# 
# 


##PLOTTING


## ---------- overlay data (compute once) ---------- ( this part can be edited to add other env variables)
temp_df <- scaled1 %>%
  filter(variable == "Temperature") %>%
  transmute(rel_dist, depth_kriged_m, temp = predicted_kriged_value) %>%
  tidyr::drop_na()

chla_df <- scaled1 %>%
  filter(variable == "chl_a_ug_l") %>%
  transmute(rel_dist, depth_kriged_m, chla = predicted_kriged_value) %>%
  tidyr::drop_na()

# robust breaks for overlays
rng_t <- range(temp_df$temp, na.rm = TRUE)
rng_c <- range(chla_df$chla, na.rm = TRUE)

brk_temp <- pretty(rng_t, n = 7); if (length(unique(brk_temp)) < 2) brk_temp <- seq(rng_t[1], rng_t[2], length.out = 7)
brk_chla <- pretty(rng_c, n = 7); if (length(unique(brk_chla)) < 2) brk_chla <- seq(rng_c[1], rng_c[2], length.out = 7)

## ---------- which abundance variables to loop over ----------
ind_vars <- scaled1 %>%
  filter(grepl("ind_m3", variable, ignore.case = TRUE)) %>%
  distinct(variable) %>%
  pull(variable) %>%
  sort()

# helpers for titles/filenames
cruise_id   <- unique(scaled1$Cruise)[1]
transect_id <- unique(scaled1$Transect_ID)[1]


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
  filter(!is.na(B_smooth), B_smooth <= max(scaled1$depth_kriged_m, na.rm = TRUE)) %>%
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


## ---------- plot builder (normal/log) ----------
make_abund_plot <- function(df, var_label, use_log = FALSE, bins_fill = 6) {
  
  # --- prepare z variable ---
  vals <- df$predicted_kriged_value
  eps  <- if (any(vals <= 0, na.rm = TRUE)) {
    pos <- vals[vals > 0]
    if (length(pos)) max(min(pos, na.rm = TRUE)/10, 1e-8) else 1e-8
  } else 0
  
  if (!use_log) {
    df <- df %>% dplyr::mutate(zplot = predicted_kriged_value)
    nonzero_vals <- df$zplot[df$zplot > 0 & is.finite(df$zplot)]
    fill_name <- expression(paste("Abundance (Ind. ", m^{-3}, ")"))
  } else {
    df <- df %>% dplyr::mutate(zplot = log10(predicted_kriged_value + eps))
    nonzero_vals <- df$zplot[df$predicted_kriged_value > 0 & is.finite(df$zplot)]
    fill_name <- expression(paste("log"[10], "(Ind. ", m^{-3}, ")"))
  }
  
  # --- compute quartiles (on appropriate scale) ---
  qtiles <- quantile(nonzero_vals, probs = c(0.25, 0.5, 0.75, 1.0), na.rm = TRUE)
  
  # --- define contour breaks with conditional zero-bin ---
  if (!use_log) {
    if (any(df$predicted_kriged_value == 0, na.rm = TRUE)) {
      breaks <- c(-1e-8, 0, qtiles, Inf)
      label_ranges <- c(
        "0",
        paste0(formatC(qtiles[1], digits = 2, format = "fg"), " – ",
               formatC(qtiles[2], digits = 2, format = "fg")),
        paste0(formatC(qtiles[2], digits = 2, format = "fg"), " – ",
               formatC(qtiles[3], digits = 2, format = "fg")),
        paste0(formatC(qtiles[3], digits = 2, format = "fg"), " – ",
               formatC(qtiles[4], digits = 2, format = "fg")),
        paste0("> ", formatC(qtiles[4], digits = 2, format = "fg"))
      )
    } else {
      breaks <- c(min(vals), qtiles, Inf)
      label_ranges <- c(
        paste0(formatC(qtiles[1], digits = 2, format = "fg"), " – ",
               formatC(qtiles[2], digits = 2, format = "fg")),
        paste0(formatC(qtiles[2], digits = 2, format = "fg"), " – ",
               formatC(qtiles[3], digits = 2, format = "fg")),
        paste0(formatC(qtiles[3], digits = 2, format = "fg"), " – ",
               formatC(qtiles[4], digits = 2, format = "fg")),
        paste0("> ", formatC(qtiles[4], digits = 2, format = "fg"))
      )
    }
  } else {
    if (any(df$predicted_kriged_value == 0, na.rm = TRUE)) {
      breaks <- c(log10(eps), log10(eps * 1.0001), qtiles, Inf)  # zero bin + quartiles
      label_ranges <- c(
        "0",
        paste0(formatC(10^qtiles[1], digits = 2, format = "fg"), " – ",
               formatC(10^qtiles[2], digits = 2, format = "fg")),
        paste0(formatC(10^qtiles[2], digits = 2, format = "fg"), " – ",
               formatC(10^qtiles[3], digits = 2, format = "fg")),
        paste0(formatC(10^qtiles[3], digits = 2, format = "fg"), " – ",
               formatC(10^qtiles[4], digits = 2, format = "fg")),
        paste0("> ", formatC(10^qtiles[4], digits = 2, format = "fg"))
      )
    } else {
      breaks <- c(min(df$zplot, na.rm = TRUE), qtiles, Inf)       # no zero bin
      label_ranges <- c(
        paste0(formatC(10^qtiles[1], digits = 2, format = "fg"), " – ",
               formatC(10^qtiles[2], digits = 2, format = "fg")),
        paste0(formatC(10^qtiles[2], digits = 2, format = "fg"), " – ",
               formatC(10^qtiles[3], digits = 2, format = "fg")),
        paste0(formatC(10^qtiles[3], digits = 2, format = "fg"), " – ",
               formatC(10^qtiles[4], digits = 2, format = "fg")),
        paste0("> ", formatC(10^qtiles[4], digits = 2, format = "fg"))
      )
    }
  }
  
  # --- build the contour plot first to inspect the levels ---
  tmp <- ggplot_build(
    ggplot(df, aes(x = rel_dist, y = depth_kriged_m)) +
      geom_contour_filled(aes(z = zplot), breaks = breaks)
  )
  n_lvls <- length(unique(tmp$data[[1]]$level))
  rm(tmp)
  
  # --- color palette adjusted dynamically ---
  vir_cols <- viridis::viridis(n_lvls, option = "D")
  vir_cols[1] <- "#08306B"  # dark blue for zero bin
  
  # --- ensure labels vector length matches levels ---
  fill_labels <- label_ranges
  if (length(fill_labels) < n_lvls) {
    fill_labels <- c(fill_labels, rep("", n_lvls - length(fill_labels)))
  }
  
  ggplot(df, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_contour_filled(aes(z = zplot), breaks = breaks) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    
    # overlays
    geom_contour(
      data = temp_df, aes(z = temp, color = "Temperature"),
      breaks = brk_temp, linewidth = 0.4
    ) +
    metR::geom_text_contour(
      data = temp_df,
      aes(z = temp, label = after_stat(level), color = "Temperature"),
      breaks = brk_temp, size = 3, stroke = 0.2,
      check_overlap = TRUE, show.legend = FALSE
    ) +
    geom_contour(
      data = chla_df, aes(z = chla, color = "Chl-a"),
      breaks = brk_chla, linewidth = 0.4, linetype = "dotted"
    ) +
    metR::geom_text_contour(
      data = chla_df,
      aes(z = chla, label = scales::number(after_stat(level), accuracy = 0.1),
          color = "Chl-a"),
      breaks = brk_chla, size = 3, stroke = 0.2,
      check_overlap = TRUE, show.legend = FALSE
    ) +
    
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_manual(
      name = fill_name,
      values = vir_cols,
      labels = fill_labels
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Temperature" = "red3", "Chl-a" = "darkgreen"),
      labels = c(
        "Temperature" = expression(paste("Temperature (", degree*C, ")")),
        "Chl-a"       = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")"))
      )
    ) +
    labs(
      title = paste0(var_label, " — ", cruise_id, " / ", transect_id),
      x = "Relative distance (dimensionless)",
      y = "Depth (m)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      legend.title = element_text(),
      legend.text  = element_text()
    ) +
    geom_polygon(
      data = bathy_poly,
      aes(x = rel_dist, y = B_smooth),
      fill = "gray40", alpha = 1
    )
}


## ---------- loop, build + save ----------
for (v in ind_vars) {
  message("Plotting: ", v)
  v_df <- scaled1 %>%
    filter(variable == v) %>%
    tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)
  
  if (nrow(v_df) < 10) {
    warning("Skipping ", v, " (too few points after NA drop).")
    next
  }
  
  # normal scale
  p_norm <- make_abund_plot(v_df, var_label = v, use_log = FALSE, bins_fill = 6)
  ggsave(
    filename = file.path(
      output_dir_figs_normal,
      paste(cruise_id, transect_id,  paste0(v, "_interp_and_krige_dist10.png"), sep = "-")
    ),
    plot   = p_norm, width = 12, height = 6, units = "in", dpi = 300
  )
  
  # log10 scale
  p_log <- make_abund_plot(v_df, var_label = paste0(v, " (log10)"), use_log = TRUE, bins_fill = 6)
  ggsave(
    filename = file.path(
      output_dir_figs_log,
      paste(cruise_id, transect_id,  paste0(v, "_log10_interp_and_krige_dist10.png"), sep = "-")
    ),
    plot   = p_log, width = 12, height = 6, units = "in", dpi = 300
  )
}

## after this is code to assess variance, might not be as useful


# plot to assess krigin variance:
ind_vars <- scaled1 %>%
  dplyr::filter(grepl("ind_m3", variable, ignore.case = TRUE)) %>%
  dplyr::distinct(variable) %>%
  dplyr::pull(variable)

env_vars_present <- intersect(env_vars, unique(scaled1$variable))

exclude_vars <- c("ind_m3_detritus", "ind_m3_NA", "ind_m3_artifact")
vars_to_plot <- sort(unique(c(ind_vars, env_vars_present)))
vars_to_plot <- setdiff(vars_to_plot, exclude_vars)

# 2) Per-variable epsilon (make sure `variable` stays in the table)
eps_tbl <- scaled1 %>%
  dplyr::filter(variable %in% vars_to_plot) %>%
  dplyr::mutate(variable = as.character(variable)) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    eps = {
      pos <- krige_variance[is.finite(krige_variance) & krige_variance > 0]
      if (length(pos)) max(min(pos, na.rm = TRUE)/10, 1e-12) else 1e-12
    },
    .groups = "drop"
  ) %>%
  dplyr::select(variable, eps)

# 3) Safe log transform per row, then summarise in log-space
var_summary_log <- scaled1 %>%
  dplyr::filter(variable %in% vars_to_plot) %>%
  dplyr::mutate(
    variable = as.character(variable),
    krige_variance = pmax(krige_variance, 0)     # clamp tiny negatives
  ) %>%
  dplyr::inner_join(eps_tbl, by = "variable") %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(log_var = log10(krige_variance + eps)) %>%
  dplyr::filter(is.finite(log_var)) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    mean_log = mean(log_var),
    sd_log   = sd(log_var),
    se_log   = sd_log / sqrt(n),
    gmean    = 10^mean_log,
    ci95_low = 10^(mean_log - 1.96 * se_log),
    ci95_high= 10^(mean_log + 1.96 * se_log),
    .groups  = "drop"
  ) %>%
  # tiny whisker if completely flat
  dplyr::mutate(
    ci95_low  = ifelse(sd_log == 0, gmean / 1.05, ci95_low),
    ci95_high = ifelse(sd_log == 0, gmean * 1.05, ci95_high),
    type      = ifelse(grepl("ind_m3", variable, ignore.case = TRUE), "Abundance", "Environment")
  ) %>%
  # order by geometric mean and lock that order in factor levels
  dplyr::arrange(gmean) %>%
  dplyr::mutate(variable = factor(variable, levels = variable))

# 4) Plot (log x-axis) with whiskers
p_var_summary <- ggplot(var_summary_log, aes(x = gmean, y = variable, color = type)) +
  geom_errorbarh(aes(xmin = ci95_low, xmax = ci95_high), height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.3) +
  scale_color_manual(values = c("Abundance" = "black", "Environment" = "grey40")) +
  scale_x_log10() +
  labs(
    title = paste0(
      "Mean kriging variance (geometric mean) — ",
      unique(scaled1$Cruise)[1], " / ", unique(scaled1$Transect_ID)[1]
    ),
    x = "Kriging variance (log scale)",
    y = NULL, color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

p_var_summary
ggsave( filename = file.path( output_dir_figs_normal, paste(cruise_id, transect_id, "krige_variance_interp_and_krige_dist10.png", sep="_")), plot = p_var_summary,
  width  = 6, height = 12, units = "in", dpi = 300
)

ggsave( filename = file.path( output_dir_figs_log, paste(cruise_id, transect_id, "krige_variance_interp_and_krige_dist10.png", sep="_")), plot = p_var_summary,
        width  = 6, height = 12, units = "in", dpi = 300
)
##----- Plot Variance-- not needed but interesting to see krigin effects
cruise_id   <- unique(scaled1$Cruise)[1]
transect_id <- unique(scaled1$Transect_ID)[1]

# 1) ind_m3* variables
ind_vars <- scaled1 %>%
  filter(grepl("ind_m3", variable, ignore.case = TRUE)) %>%
  distinct(variable) %>%
  pull(variable)

# 2) env vars that are present in 'scaled1$variable'
env_vars_present <- intersect(env_vars, unique(scaled1$variable))

# all variables to plot
vars_to_plot <- sort(unique(c(ind_vars, env_vars_present)))

# helper: legend label (variance units)
var_variance_label <- function(v) {
  if (grepl("ind_m3", v, ignore.case = TRUE)) {
    return(expression(paste("Kriging variance (Ind. ", m^{-3}, ")^2")))
  }
  switch(v,
         "Temperature"  = expression(paste("Kriging variance (", degree*C, ")^2")),
         "Salinity"     = expression(paste("Kriging variance (PSU)^2")),
         "Oxygen_ml_l"  = expression(paste("Kriging variance (mL ", L^{-1}, ")^2")),
         "pH"           = expression(paste("Kriging variance (unitless)^2")),
         "PAR"          = expression(paste("Kriging variance (unitless)^2")),
         "chl_a_ug_l"   = expression(paste("Kriging variance (", mu, "g ", L^{-1}, ")^2")),
         # default fallback:
         expression("Kriging variance")
  )
}

# --- loop: plot kriging variance for each selected variable ---
for (v in vars_to_plot) {
  message("Plotting variance: ", v)
  v_df <- scaled1 %>%
    filter(variable == v) %>%
    tidyr::drop_na(krige_variance, depth_kriged_m, rel_dist)
  
  # light guard against degenerate fields
  ok_enough <- nrow(v_df) >= 20 &&
    dplyr::n_distinct(v_df$rel_dist) >= 3 &&
    dplyr::n_distinct(v_df$depth_kriged_m) >= 3 &&
    is.finite(stats::sd(v_df$krige_variance, na.rm = TRUE)) &&
    stats::sd(v_df$krige_variance, na.rm = TRUE) > .Machine$double.eps
  if (!ok_enough) {
    warning("Skipping variance for ", v, " (insufficient grid or near-constant variance).")
    next
  }
  
  # discretize for readability
  brks_var <- pretty(range(v_df$krige_variance, na.rm = TRUE), n = 8)
  
  p_var <- ggplot(v_df, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_contour_filled(aes(z = krige_variance), breaks = brks_var) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_viridis_d(name = var_variance_label(v)) +
    labs(
      title = paste0(v, " — Kriging variance — ", cruise_id, " / ", transect_id),
      x = "Relative distance (dimensionless)",
      y = "Depth (m)"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right")
  
}



