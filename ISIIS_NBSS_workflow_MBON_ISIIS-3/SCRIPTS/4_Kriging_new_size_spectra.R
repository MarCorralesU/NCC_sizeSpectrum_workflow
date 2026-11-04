# ============================
# ADAPTED KRIGING PIPELINE (your w23_NH "dat")
# keeps original functionality; only maps columns & paths
# ============================
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
base_dir<- getwd()

# --- 0) Paths & identifiers (ADAPT: uses your base_dir and w23_NH) ----
# main directory for outputs (pick your SPECTRA workflow root)
main_dir <- file.path(base_dir, "ISIIS_NBSS_workflow_SPECTRA")
source(file.path(main_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


option_list <- list(
  make_option(c("-c", "--Cruise"), type = "character", help = "Cruise: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t", "--Transect"), type = "character", help = "Cruise: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$Cruise)) {
  stop("Error: --Cruise is required.")
}

#variables to load data:
Cruise <- opt$Cruise # Cruise<- 'w23'
Transect <- opt$Transect # Transect<- 'NH'

size_var<- opt$sizeVar # size_var<-'contour'
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(main_dir, opt$distFromShore) # shore_path  <- file.path(main_dir, "ancillary", "transect_shelf_land.csv")
sizeSpectraDir <- file.path(base_dir, 'ISIIS_NBSS_workflow_SPECTRA', 'data', 'size_spectra', size_var, Cruise, Transect) #use this for testing sizeSpectraDir <- file.path(base_dir, 'ISIIS_NBSS_workflow_SPECTRA', 'data', 'size_spectra', 'contour', 'w23', 'NH')


df_list <- list.files(sizeSpectraDir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)


# Create output folders 
output_dir_data <- file.path(main_dir, "data", "interp_and_krige_dist10", 'size_spectra_coeffs')
if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_data)
}

output_dir_figs_normal <- file.path(main_dir, "figures", "interp_and_krige_dist10", 'size_spectra_coeffs', Cruise, Transect)
if (!dir.exists(output_dir_figs_normal)) {
  dir.create(output_dir_figs_normal, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_normal)
}





print(sizeSpectraDir)

dat <- map_dfr(
  df_list,
  ~ read_csv(
    .x,
    show_col_types = FALSE,
    col_select = c(
      time_Pacific, rel_dist, depth_mid, Latitude, Longitude,
      slope, intercept, dplyr::any_of(env_vars)
    )
  ) |> distinct()
) %>%
  mutate(slope = as.numeric(slope)) %>%   # ensure numeric for comparison
  filter(!is.na(slope) & slope < 0) %>%   # exclude all rows where slope >= 0
  arrange(time_Pacific)

dat_ord <- dat %>% arrange(Longitude, time_Pacific)  # break ties by time if present

binned_plot<- ggplot(dat_ord, aes(x = Longitude, y = depth_mid)) +
  geom_point(aes(color = slope), alpha = 0.6, size = 1) +
  scale_y_reverse() +
  scale_color_viridis_c(option = "viridis") +
  labs(x = "Longitude", y = "Depth (m)", color = "Slope") +
  theme_classic(base_size = 12)

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "binned_data_plot.png", sep = "-")),
       plot = binned_plot, width  = 12, height = 6, units  = "in", dpi = 300)

# ---------- 1) Build myDf from your `dat` (use Longitude -> distance km) ----------
# Expect: dat has at least: time_Pacific, rel_dist, depth_mid, Latitude, Longitude, slope

# We'll compute distance in km from the westernmost longitude (start of transect):
# This preserves the original pipeline's functionality (distance in km along x).

# Use dplyr:: to avoid plyr clashes
myDf <- dat %>%
  dplyr::mutate(
    Transect_ID = Cruise,
    HD_ID       = Transect,
    Depth       = as.numeric(depth_mid),         # vertical axis (m)
    Longitude   = as.numeric(Longitude),
    Latitude    = as.numeric(Latitude),
    slope       = as.numeric(slope),
    intercept   = as.numeric(intercept)
  ) %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, Latitude, Longitude,
    dplyr::any_of(env_vars),
    slope, intercept,  time_Pacific
  )

# Compute km-per-degree at mean latitude, and define the start longitude (easternmost)
lat0 <- mean(myDf$Latitude, na.rm = TRUE)
km_per_deg_lon <- 111.320 * cos(lat0 * pi/180)

lon0 <- max(myDf$Longitude, na.rm = TRUE)   # start = easternmost point

# Distance grows westward (more negative longitudes):
#   east -> west : GPS.dist.from.start.km increases
myDf <- myDf %>%
  dplyr::mutate(GPS.dist.from.start.km = (lon0 - Longitude) * km_per_deg_lon)

# Add along-transect distance (km) from start longitude; keep Longitude in the table
myDf <- myDf %>%
  dplyr::mutate(GPS.dist.from.start.km = (Longitude - lon0) * km_per_deg_lon) %>%
  dplyr::select(
    Transect_ID, HD_ID, Depth, GPS.dist.from.start.km,
    dplyr::any_of(env_vars),
    slope, intercept, Latitude, Longitude, time_Pacific
  )

# Drop odd columns if they exist (same as original script)
if ("NA" %in% names(myDf)) myDf$`NA` <- NULL
if ("NA_mean_area"  %in% names(myDf)) myDf$NA_mean_area  <- NULL
if ("NA_mean_major" %in% names(myDf)) myDf$NA_mean_major <- NULL
if ("NA_mean_minor" %in% names(myDf)) myDf$NA_mean_minor <- NULL
gc()


# ---------- 2) Long variable table (variables to krige; here: slope) ----------
biophys.10 <- myDf %>%
  dplyr::select(Transect_ID, HD_ID, Depth, GPS.dist.from.start.km, dplyr::any_of(c("slope", "intercept")))

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


# # ---------- 4) Coarse interpolation (akima) per variable (xo 0.5 km, yo 5 m) ----------
# 
# 
# 
#   
#   
#   
#   
# 


# ---------- 5) Kriging (fine grid; Exp, maxdist=10) ----------
krige_biophys1 <- NULL
variable2 <- unique(biophys.14$variable)


for (subvar in variable2) {
  x <- biophys.14[biophys.14$variable == subvar, ] %>% stats::na.omit()
  if (nrow(x) == 0) next
  
  print(unique(x$variable))
  
  x0 <- seq(min(x$GPS.dist.from.start.km, na.rm = TRUE),
            max(x$GPS.dist.from.start.km, na.rm = TRUE), by = 0.5)
  y0 <- seq(0, max(x$Depth, na.rm = TRUE), by = 2)
  
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
  x_kriged_df$Longitude <- (x_kriged_df$x / km_per_deg_lon) + lon0
  
  x_kriged_df$variable    <- subvar
  x_kriged_df$HD_ID       <- Transect
  x_kriged_df$Transect_ID <- Cruise
  
  krige_biophys1 <- rbind(krige_biophys1, x_kriged_df)
}

message("Kriged")


# ---------- 6) Rescale predictions back to original [min, max] per variable ----------
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
  scaled1$Longitude <- (scaled1$x / km_per_deg_lon) + lon0
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

write.csv(scaled1, file.path(output_dir_data, paste(Cruise, Transect, size_var, 'size-spectrum_envs_interp_and_krige_dist10.csv', sep='-')), row.names=FALSE)


# ---------- 8) Plots (same behavior as original) ----------
# base (filled) data for slope
s_df <- scaled1 %>%
  filter(variable == "slope") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)

# overlay data: Temperature + Chl-a
temp_df <- scaled1 %>%
  filter(variable == "Temperature") %>%
  transmute(rel_dist, depth_kriged_m, temp = predicted_kriged_value) %>%
  tidyr::drop_na()

chla_df <- scaled1 %>%
  filter(variable == "chl_a_ug_l") %>%
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


p_slope <- ggplot(s_df, aes(x = rel_dist, y = depth_kriged_m)) +
  # filled contours of SLOPE
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 6) +
  # shelf-break marker
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
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
  scale_y_reverse(expand = c(0, 0)) +
  # slope fill legend (viridis discrete)
  scale_fill_viridis_d(
    name = expression(paste("Slope (", mu, "m"^{-3}, " ", m^{-3}, ")"))
  ) +
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

p_slope
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "slope_interp_and_krige_dist10.png", sep = "-")),
       plot = p_slope, width  = 12, height = 6, units  = "in", dpi = 300)

i_df <- scaled1 %>%
  dplyr::filter(variable == "intercept") %>%
  tidyr::drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)


p_intercept <- ggplot(i_df, aes(x = rel_dist, y = depth_kriged_m)) +
  # filled contours for INTERCEPT
  geom_contour_filled(aes(z = predicted_kriged_value), bins = 6) +
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
  # intercept fill legend
  scale_fill_viridis_d(
    name = expression(paste("Intercept (", mu, "m"^3, " ", mu, "m"^{-3}, " ", m^{-3}, ")"))
  ) +
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

p_intercept

ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "intercept_interp_and_krige_dist10.png", sep = "-")),
       plot = p_intercept, width  = 12, height = 6, units  = "in", dpi = 300)


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

p_chla
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "chla_interp_and_krige_dist10.png", sep = "-")),
       plot = p_chla, width  = 12, height = 6, units  = "in", dpi = 300)

t_df <- scaled1 %>%
  filter(variable == "Temperature") %>%
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

p_temp
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "temperature_interp_and_krige_dist10.png", sep = "-")),
       plot = p_temp, width  = 12, height = 6, units  = "in", dpi = 300)

sali_df <- scaled1 %>%
  filter(variable == "Salinity") %>%
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

p_sali
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "salinity_interp_and_krige_dist10.png", sep = "-")),
       plot = p_sali, width  = 12, height = 6, units  = "in", dpi = 300)

o2_df <- scaled1 %>%
  filter(variable == "Oxygen_ml_l") %>%
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

p_o2
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "oxygen_interp_and_krige_dist10.png", sep = "-")),
       plot = p_o2, width  = 12, height = 6, units  = "in", dpi = 300)

par_df <- scaled1 %>%
  filter(variable == "PAR") %>%
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

p_par
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "PAR_interp_and_krige_dist10.png", sep = "-")),
       plot = p_par, width  = 12, height = 6, units  = "in", dpi = 300)

ph_df <- scaled1 %>%
  filter(variable == "pH") %>%
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

p_ph
ggsave(filename = file.path(output_dir_figs_normal, paste(Cruise, Transect, size_var, "ph_interp_and_krige_dist10.png", sep = "-")),
       plot = p_ph, width  = 12, height = 6, units  = "in", dpi = 300)



