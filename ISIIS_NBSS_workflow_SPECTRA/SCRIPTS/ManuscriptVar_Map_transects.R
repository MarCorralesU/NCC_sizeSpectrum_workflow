# ---- packages
setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)       # for bathymetry / depth rasters
library(viridis)      # nice color palettes
library(ggOceanMaps)
library(terra)
library(marmap)
library(readr)
library(dplyr)
library(purrr)
library(stringr)

base_dir<- getwd()



# grab all CSVs
files <- list.files(
  file.path(base_dir, "data", "binned_abundance"),
  pattern = "\\.csv$",
  full.names = TRUE
) |> 
  (\(x) x[!grepl("OB", x)])()
# filename pattern: e.g. "CM_s22_binned_abundance_test_chunk_1.csv"
pat <- "^([A-Za-z]+)_([sw]\\d{2})_.*\\.csv$"  # capture Transect_ID and Cruise (w22/s23)

# helper to find column by regex (case-insensitive)
pick_col <- function(nms, pattern) {
  hits <- nms[str_detect(nms, regex(pattern, ignore_case = TRUE))]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

summary_df <- map_dfr(files, function(f) {
  bn <- basename(f)
  m <- str_match(bn, pat)
  if (is.na(m[1, 1])) stop("Filename didn't match expected pattern: ", bn)
  
  transect_id <- m[1, 2]
  cruise      <- m[1, 3]
  
  dat <- suppressMessages(read_csv(f, show_col_types = FALSE, progress = FALSE))
  
  lat_col <- pick_col(names(dat), "^mean_?lat(itude)?$")
  lon_col <- pick_col(names(dat), "^mean_?lon(gitude)?$")
  
  if (is.na(lat_col) || is.na(lon_col)) {
    stop("Could not find mean_Latitude/mean_Longitude in: ", bn)
  }
  
  tibble(
    Transect_ID = transect_id,
    Cruise      = cruise,
    mean_Lat    = mean(dat[[lat_col]], na.rm = TRUE),
    min_Lon     = min(dat[[lon_col]],  na.rm = TRUE),
    max_Lon     = max(dat[[lon_col]],  na.rm = TRUE)
  )
}) %>%
  arrange(Transect_ID, Cruise)

summary_by_transect <- summary_df %>%
  dplyr::summarise(
    mean_Lat = mean(mean_Lat, na.rm = TRUE),
    min_Lon  = min(min_Lon,  na.rm = TRUE),
    max_Lon  = max(max_Lon,  na.rm = TRUE),
    .by = Transect_ID
  ) %>%
  dplyr::arrange(Transect_ID)

transect_names <- c(
  LP = 'La Push',
  GH = 'Grays Harbor',
  OB=  'Ocean Beaches',
  CM = "Cape Meares",
  CR = "Columbia River",
  HH = "Heceta Head",
  NH = "Newport Line",
  RR = "Rogue River",
  CC = "Crescent City",
  TR = "Trinidad"
)


lims <- c(-126.5, -122.2, 42, 47.5)  # xmin, xmax, ymin, ymax

b <- marmap::getNOAA.bathy(  lon1 = lims[1], lon2 = lims[2],
                             lat1 = lims[3], lat2 = lims[4],
                             resolution = 0.01)            # 15″
r <- marmap::as.raster(b)# RasterLayer
r_terra <- terra::rast(r)

fine_template <- terra::rast(
  xmin = terra::ext(r_terra)[1],
  xmax = terra::ext(r_terra)[2],
  ymin = terra::ext(r_terra)[3],
  ymax = terra::ext(r_terra)[4],
  resolution = 0.01,
  crs = terra::crs(r_terra)
)

r_smooth <- terra::resample(r_terra, fine_template, method = "bilinear")

outfile <- "~/ggOceanMapsLargeData/etopo15s_PNW.tif"
writeRaster(r_smooth, outfile, overwrite = TRUE)
outfile <- "~/ggOceanMapsLargeData/etopo15s_PNW.tif"
writeRaster(terra::rast(r), outfile, overwrite = TRUE)

# 2) Tell ggOceanMaps where that file is, and plot with the user style
options(ggOceanMaps.userpath = outfile)



# summary_by_transect should already have: Transect_ID, mean_Lat, min_Lon, max_Lon
# build a plotting df + a safe x for labels within the x-limits
x_span <- diff(lims[1:2])
lines_df <- summary_by_transect %>%
  mutate(
    label   = recode(Transect_ID, !!!transect_names),
    x0      = min_Lon,
    x1      = max_Lon,
    y       = mean_Lat,
    x_label = pmin(x1 + 0.04 * x_span, lims[2] - 0.01 * x_span)  # nudge right, keep inside panel
  )

# ---- stepped white→blue palette for your custom breaks
brks <- c(0, 200, 1000, 2000, 3000)  # edges
pal_blues <- grDevices::colorRampPalette(
  c("white", "#deebf7", "#9ecae1", "#3182bd", "#08519c")
)(length(brks) - 1)

# ---- rebuild base map and override its fill scale
p <- basemap(
  limits = lims,
  crs = 4326,
  bathymetry = TRUE,
  bathy.style = "rub",     # using your high-res user raster
  downsample = FALSE,
  projection.grid = FALSE,
  land.col = "grey40",
  land.border.col = "grey20"
)

# If your raster has **positive** depths (0..3000):
p <- p + ggplot2::scale_fill_stepsn(
  name   = "Depth (m)",
  breaks = brks,
  limits = c(min(brks), max(brks)),
  colors = pal_blues,
  oob    = scales::squish,                 # >3000 squished to last bin color
  guide  = ggplot2::guide_colorsteps(show.limits = TRUE)
)

# If your raster has **negative** depths (typical NOAA: -3000..0), use this instead:

# ---- (re)add your transect lines & labels on top

x_span <- diff(lims[1:2])
lines_df <- summary_by_transect |>
  dplyr::mutate(
    label   = dplyr::recode(Transect_ID, !!!transect_names),
    x0      = min_Lon,
    x1      = max_Lon,
    y       = mean_Lat,
    x_label = pmin(x1 + 0.04 * x_span, lims[2] - 0.01 * x_span)
  )

p <- p +
  ggplot2::geom_segment(
    data = lines_df,
    aes(x = x0, xend = x1, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 2,
    lineend = "round",
    color = "black"
  ) +
  ggplot2::geom_text(
    data = lines_df,
    aes(x = x_label, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5,
    size = 6,
    color = "black"
  )

p <- p +
  labs(x = "Lon", y = "Lat") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p

#OPTIONAL: add comparisons with  ESMs spatial resolution

# ---- Add reference squares with true geographic scaling ----
library(ggnewscale)

# Conversion helpers
km_to_latdeg <- function(km) km / 111.32
km_to_londeg <- function(km, lat) km / (111.32 * cos(lat * pi / 180))

# Define a reference point (offshore central Oregon)
center_lat <- 44.3
center_lon <- -125.0

# 1️⃣ 4×4 km "Satellite products" box


# 2️⃣ 1×1 degree "Earth System Model output" cell
half_deg <- 0.5   # half side = 0.5° in both lat and lon
square_1deg <- data.frame(
  lon = c(center_lon - half_deg, center_lon + half_deg,
          center_lon + half_deg, center_lon - half_deg, center_lon - half_deg),
  lat = c(center_lat - half_deg, center_lat - half_deg,
          center_lat + half_deg, center_lat + half_deg, center_lat - half_deg),
  type = "Earth System Model output (1×1°)"
)

squares_df <- rbind(square_1deg) #square_4km 

# ---- Reset fill scale so this doesn't conflict with bathymetry ----
p <- p + ggnewscale::new_scale_fill()

# ---- Add squares to the map ----
p <- p +
  geom_polygon(
    data = squares_df,
    aes(x = lon, y = lat, group = type, fill = type),
    color = "black",
    alpha = 0.4,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Satellite products (4×4 km)" = "green3",
      "Earth System Model output (1×1°)" = "orange"
    ),
    guide = guide_legend(override.aes = list(alpha = 0.6))
  )
p


ggsave(filename = file.path(base_dir, "figures", "map_transects.png"),
       plot = p, width  = 6, height = 12, units  = "in", dpi = 300)


# ============================================================
# Separate high-resolution transect map using geom_contour_filled
# ============================================================

library(maps)

# ---- limits (same as your current map)
lims <- c(-126.3, -122.2, 42, 47.5)  # xmin, xmax, ymin, ymax

# ---- states / land polygon
states <- map_data("state") |>
  dplyr::filter(region %in% c("oregon", "washington", "california"))

# ---- bathymetry
bathy <- marmap::getNOAA.bathy(
  lon1 = lims[1], lon2 = lims[2],
  lat1 = lims[3], lat2 = lims[4],
  resolution = 0.01
)

bathy.xyz <- marmap::as.xyz(bathy)
colnames(bathy.xyz) <- c("x", "y", "z")
bathy.xyz <- bathy.xyz[bathy.xyz$z <= 0, ]   # keep ocean only

# ---- same transect labels / lines df as before
x_span <- diff(lims[1:2])

lines_df <- summary_by_transect |>
  dplyr::mutate(
    label   = dplyr::recode(Transect_ID, !!!transect_names),
    x0      = min_Lon,
    x1      = max_Lon,
    y       = mean_Lat,
    x_label = pmin(x1 + 0.04 * x_span, lims[2] - 0.01 * x_span)
  )

# ---- bathymetry breaks (same idea as before, but negative and increasing)
bathy_breaks <- c(-3500, -2000, -1000, -500, -200, -100, 0)

# ---- palette matching your previous style
contour_levels <- c(
  "(-3500, -2000]",
  "(-2000, -1000]",
  "(-1000, -500]",
  "(-500, -200]",
  "(-200, -100]",
  "(-100, 0]"
)

pal_blues <- c(
  "#08519c",
  "#2171b5",
  "#6baed6",
  "#9ecae1",
  "#deebf7",
  "#f7fbff"
)

# ---- build hi-res contour map
p_hiRes <- ggplot() +
  
  # filled bathymetry
  geom_contour_filled(
    data = bathy.xyz,
    aes(x = x, y = y, z = z),
    breaks = bathy_breaks
  ) +
  scale_fill_manual(
    name = "Depth (m)",
    values = setNames(pal_blues, contour_levels),
    breaks = contour_levels,
    labels = c("3500–2000", "2000–1000", "1000–500", "500–200", "200–100", "100–0"),
    drop = FALSE
  ) +
  
  # optional contour line at -200 m
  geom_contour(
    data = bathy.xyz,
    aes(x = x, y = y, z = z),
    breaks = -200,
    color = "black",
    linewidth = 0.6
  ) +
  
  # land
  geom_polygon(
    data = states,
    aes(x = long, y = lat, group = group),
    fill = "grey40",
    color = "grey20",
    linewidth = 0.3
  ) +
  
  # transect lines
  geom_segment(
    data = lines_df,
    aes(x = x0, xend = x1, y = y, yend = y),
    linewidth = 2,
    lineend = "round",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  # transect labels
  geom_text(
    data = lines_df,
    aes(x = x_label, y = y, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  coord_cartesian(
    xlim = lims[1:2],
    ylim = lims[3:4],
    expand = FALSE
  ) +
  
  scale_x_continuous(
    breaks = c(-126, -125, -124, -123, -122),
    labels = function(x) paste0(abs(x), "°W")
  ) +
  scale_y_continuous(
    breaks = c(42, 43, 44, 45, 46, 47),
    labels = function(y) paste0(y, "°N")
  ) +
  
  labs(x = "Longitude °W", y = "Latitude °N") +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

p_hiRes

ggsave( filename = file.path(base_dir, "figures", "map_transects_hiRes.png"), plot = p_hiRes,
  width = 6, height = 9, units = "in", dpi = 300 )
