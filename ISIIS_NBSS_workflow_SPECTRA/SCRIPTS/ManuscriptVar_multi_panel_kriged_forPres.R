#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(colorspace)
  library(metR)
  library(optparse)
  library(viridis)
  library(patchwork)
  library(cowplot)
})

# ============================================================
# Multi-panel kriged figures assembled from per-panel plots
# This preserves the appearance of the original single-panel
# kriged plots much better than faceting.
#
# FOR-PRESS VERSION
# Data source is fixed to ISIIS plankton.
# Only three multi-panel figures are produced:
#   1) slope
#   2) integrated_biovol_all
#   3) integrated_biovol_sample
#
# All three use continuous viridis fills, interpolated rasters,
# contour isolines, and a horizontal legend above the panel grid.
#
# rows    = GH, CR, CM, NH, HH, RR
# columns = w22, s22, w23, s23
# ============================================================

# ---------------------------
# Paths
# ---------------------------
base_dir <- getwd()

# ---------------------------
# Options
# ---------------------------
# This for-press script is intentionally fixed to ISIIS plankton.
# The only optional argument is an alternate output directory.
option_list <- list(
  make_option(
    c("-o", "--outputDir"),
    type = "character",
    default = NULL,
    help = "Optional output directory"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data_source <- "plankton"
instrument_dir <- "ISIIS"

input_dir <- file.path(
  base_dir, "data", "interp_and_krige_dist10",
  "ISIIS", "plankton", "size_spectra_coeffs"
)

output_dir <- if (is.null(opt$outputDir)) {
  file.path(
    base_dir, "figures", "interp_and_krige_dist10",
    "multi_panel_kriged_preserved_panels", "ISIIS"
  )
} else {
  opt$outputDir
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading kriged files from: ", input_dir)
message("Writing figures to: ", output_dir)

# NOTE: This multi-panel script reads already-kriged outputs and does not contain
# a raw/binned-data plot. Therefore no binned-data color scale is modified here.
# The three for-press multi-panel outputs below all use viridis.

# ---------------------------
# Constants
# ---------------------------
transect_levels <- c("GH", "CR", "CM", "NH", "HH", "RR")
cruise_levels   <- c("w22", "s22", "w23", "s23")

metric_vars <- c(
  "slope",
  "slope_anomaly",
  "intercept",
  "total_biovolume_sum",
  "integrated_biovol_all",
  "integrated_biovol_sample",
  "biovol_diff"
)

env_vars <- c(
  "Temperature",
  "Salinity",
  "Oxygen_ml_l",
  "PAR",
  "pH",
  "chl_a_ug_l"
)

# Only the three for-press variables are built and saved.
target_vars <- c(
  "slope",
  "integrated_biovol_all",
  "integrated_biovol_sample"
)

cruise_label_map <- c(
  w22 = "Winter 2022",
  s22 = "Summer 2022",
  w23 = "Winter 2023",
  s23 = "Summer 2023"
)

transect_label_map <- c(
  GH = "Grays Harbor",
  CR = "Columbia River",
  CM = "Cape Meares",
  NH = "Newport Line",
  HH = "Heceta Head",
  RR = "Rogue River"
)


fmt_sci <- scales::label_scientific(digits = 2)

# ---------------------------
# Read all kriged outputs
# ---------------------------
csv_files <- list.files(
  input_dir,
  pattern = "size-spectrum_envs_interp_and_krige_dist10\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(csv_files) == 0) {
  stop("No kriged CSV files found in: ", input_dir)
}

kriged_all <- purrr::map_dfr(csv_files, function(f) {
  readr::read_csv(f, show_col_types = FALSE)
})

#read the no data masks
mask_files <- list.files(
  input_dir,
  pattern = "no_data_mask_interp_and_krige_dist10\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

#      no_data
#  NULL

mask_all <- if (length(mask_files) > 0) {
  
  purrr::map_dfr(mask_files, function(f) {
    readr::read_csv(
      f,
      show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    )
  }) %>%
    dplyr::mutate(
      Cruise = as.character(Cruise),
      Transect_ID = as.character(Transect_ID),
      dist_kriged_km = as.numeric(dist_kriged_km),
      depth_kriged_m = as.numeric(depth_kriged_m),
      Longitude = as.numeric(Longitude),
      Latitude = as.numeric(Latitude),
      rel_dist = as.numeric(rel_dist),
      no_data = dplyr::case_when(
        no_data %in% c("TRUE", "True", "true", "1") ~ TRUE,
        no_data %in% c("FALSE", "False", "false", "0") ~ FALSE,
        TRUE ~ NA
      )
    ) %>%
    dplyr::filter(
      Cruise %in% cruise_levels,
      Transect_ID %in% transect_levels,
      no_data
    )
  
} else {
  NULL
}

#helper function to add the no data mask
get_mask_panel <- function(mask_df, cruise, transect) {
  if (is.null(mask_df) || nrow(mask_df) == 0) return(NULL)
  
  out <- mask_df %>%
    dplyr::filter(Cruise == cruise, Transect_ID == transect, no_data)
  
  if (nrow(out) == 0) return(NULL)
  out
}

required_cols <- c(
  "Cruise", "Transect_ID", "variable",
  "rel_dist", "depth_kriged_m", "predicted_kriged_value"
)

missing_cols <- setdiff(required_cols, names(kriged_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns in kriged data: ", paste(missing_cols, collapse = ", "))
}

kriged_all <- kriged_all %>%
  mutate(
    Cruise = as.character(Cruise),
    Transect_ID = as.character(Transect_ID),
    rel_dist = as.numeric(rel_dist),
    depth_kriged_m = as.numeric(depth_kriged_m),
    predicted_kriged_value = as.numeric(predicted_kriged_value)
  ) %>%
  filter(
    variable %in% target_vars,
    Cruise %in% cruise_levels,
    Transect_ID %in% transect_levels
  )

# ---------------------------
# Bathymetry helper
# ---------------------------
build_bathy_poly <- function(cruise, transect, max_depth) {
  
  bathy_file <- file.path(
    base_dir, "data", "bathymetry_polygons",
    paste(cruise, transect, "bathymetry.csv", sep = "-")
  )
  
  if (!file.exists(bathy_file)) {
    return(NULL)
  }
  
  bathy_df <- readr::read_csv(bathy_file, show_col_types = FALSE)
  
  if (!all(c("rel_dist", "B_smooth") %in% names(bathy_df))) {
    return(NULL)
  }
  
  bathy_df <- bathy_df %>%
    arrange(desc(rel_dist)) %>%
    mutate(B_smooth_monotonic = cummax(B_smooth)) %>%
    rename(B_smooth_original = B_smooth) %>%
    rename(B_smooth = B_smooth_monotonic)
  
  bathy_poly <- bathy_df %>%
    dplyr::select(rel_dist, B_smooth, dplyr::any_of("Longitude")) %>%
    filter(!is.na(B_smooth), B_smooth <= max_depth) %>%
    arrange(rel_dist)
  
  if (nrow(bathy_poly) == 0) return(NULL)
  
  bathy_poly <- rbind(
    bathy_poly,
    data.frame(
      rel_dist = max(bathy_poly$rel_dist, na.rm = TRUE),
      B_smooth = max_depth,
      Longitude = NA_real_
    ),
    data.frame(
      rel_dist = min(bathy_poly$rel_dist, na.rm = TRUE),
      B_smooth = max_depth,
      Longitude = NA_real_
    )
  )
  
  bathy_poly
}

# ---------------------------
# Breaks
# ---------------------------
get_slope_breaks <- function(data_source) {
  if (data_source == "plankton") {
    c(-Inf, seq(-1.8, 0, by = 0.25), Inf)
  } else {
    c(-Inf, seq(-1.1, -0.2, by = 0.1), Inf)
  }
}

get_intercept_breaks <- function(data_source) {
  if (data_source == "plankton") {
    c(-Inf, seq(5, 20, by = 1.6), Inf)
  } else {
    c(-Inf, seq(6.5, 11.5, by = 1), Inf)
  }
}

get_biovol_breaks <- function(data_source) {
  if (data_source == "plankton") {
    c(-Inf, 10^(9:20), Inf)
  } else {
    c(-Inf, 10^(seq(11, 15, by = 0.3)), Inf)
  }
}


get_env_palette <- function(var_name, n) {
  switch(
    var_name,
    "Temperature" = rev(RColorBrewer::brewer.pal(n, "RdYlBu")),
    "Salinity"    = rev(RColorBrewer::brewer.pal(n, "RdYlBu")),
    "Oxygen_ml_l" = rev(RColorBrewer::brewer.pal(n, "RdYlBu")),
    "PAR"         = rev(RColorBrewer::brewer.pal(n, "RdYlBu")),
    "pH"          = rev(RColorBrewer::brewer.pal(n, "RdYlBu")),
    "chl_a_ug_l"  = viridis::viridis(n)
  )
}


# ---------------------------
# Label helpers
# ---------------------------
make_interval_labels <- function(breaks_vec, fmt_fun = function(x) x) {
  lefts  <- head(breaks_vec, -1)
  rights <- breaks_vec[-1]
  
  lab_left  <- ifelse(is.infinite(lefts),  "-Inf", fmt_fun(lefts))
  lab_right <- ifelse(is.infinite(rights), "Inf",  fmt_fun(rights))
  
  paste(lab_left, lab_right, sep = " – ")
}

bin_values <- function(x, breaks_vec, labels_vec) {
  cut(
    x,
    breaks = breaks_vec,
    include.lowest = TRUE,
    right = TRUE,
    labels = labels_vec
  )
}



# ---------------------------
# Theme helpers
# ---------------------------
panel_theme <- function(show_x = FALSE, show_y = FALSE) {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, hjust = 0, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x  = if (show_x) element_text(size = 14) else element_blank(),
      axis.text.y  = if (show_y) element_text(size = 14) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      legend.position = "right",
      plot.margin = margin(4, 4, 4, 4)#plot.margin = margin(10, 5, 16, 5)
    )
}

blank_panel <- function() {
  ggplot() +
    theme_void() +
    theme(
      plot.margin = margin(5, 5, 5, 5)
    )
}

column_label_plot <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(2, 2, 6, 2)
    )
}

row_label_plot <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 4.6, angle = 270) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(2, 6, 2, 6)
    )
}

get_env_breaks <- function(df_all, var_name, bins = 8) {
  x <- df_all %>%
    dplyr::filter(variable == var_name) %>%
    dplyr::pull(predicted_kriged_value)
  
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NULL)
  
  seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = bins + 1)
}

make_contour_level_labels <- function(breaks_vec, digits = 2) {
  lefts  <- head(breaks_vec, -1)
  rights <- tail(breaks_vec, -1)
  
  fmt <- function(x) formatC(x, format = "f", digits = digits)
  
  paste0("(", fmt(lefts), ", ", fmt(rights), "]")
}

# ---------------------------
# Plot builders for one panel
# ---------------------------
make_discrete_panel <- function(df_panel, cruise, transect,
                                variable, breaks_vec, fill_name,
                                labels_vec, palette_vals,
                                title_text = NULL,
                                show_x = FALSE, show_y = FALSE) {
  
  if (is.null(df_panel) || nrow(df_panel) == 0) {
    return(blank_panel())
  }
  
  max_depth <- max(df_panel$depth_kriged_m, na.rm = TRUE)
  bathy_poly <- build_bathy_poly(cruise, transect, max_depth)
  mask_panel <- get_mask_panel(mask_all, cruise, transect)
  
  df_plot <- df_panel %>%
    mutate(fill_bin = bin_values(predicted_kriged_value, breaks_vec, labels_vec))
  
  p <- ggplot(df_plot, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_raster(aes(fill = fill_bin)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 1.2
    ) +
    scale_x_continuous(breaks = get_rel_dist_breaks(df_plot$rel_dist), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_manual(
      values = stats::setNames(palette_vals, labels_vec),
      limits = labels_vec,
      drop = FALSE,
      name = fill_name,
      na.value = "white"
    ) +
    labs(
      title = title_text,
      x = if (show_x) "Relative distance" else NULL,
      y = if (show_y) "Depth (m)" else NULL
    ) +
    panel_theme(show_x = show_x, show_y = show_y) +
    theme(legend.position = "none")
  
  if (!is.null(mask_panel) && nrow(mask_panel) > 0) {
    p <- p +
      geom_tile(
        data = mask_panel,
        aes(x = rel_dist, y = depth_kriged_m),
        inherit.aes = FALSE,
        fill = "white",
        alpha = 0.9
      )
  }
  
  if (!is.null(bathy_poly) && nrow(bathy_poly) > 0) {
    p <- p +
      geom_polygon(
        data = bathy_poly,
        aes(x = rel_dist, y = B_smooth),
        inherit.aes = FALSE,
        fill = "gray40",
        alpha = 1
      )
  }
  
  p
}

make_continuous_panel <- function(df_panel, cruise, transect,
                                  variable, scale_fun,
                                  title_text = NULL,
                                  show_x = FALSE, show_y = FALSE,
                                  add_contours = FALSE,
                                  contour_bins = 5,
                                  contour_label_accuracy = 0.1) {
  
  if (is.null(df_panel) || nrow(df_panel) == 0) {
    return(blank_panel())
  }
  
  max_depth <- max(df_panel$depth_kriged_m, na.rm = TRUE)
  bathy_poly <- build_bathy_poly(cruise, transect, max_depth)
  mask_panel <- get_mask_panel(mask_all, cruise, transect)
  
  p <- ggplot(df_panel, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_raster(aes(fill = predicted_kriged_value)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 1.2
    ) +
    scale_x_continuous(breaks = get_rel_dist_breaks(df_panel$rel_dist), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fun +
    labs(
      title = title_text,
      x = if (show_x) "Relative distance" else NULL,
      y = if (show_y) "Depth (m)" else NULL
    ) +
    panel_theme(show_x = show_x, show_y = show_y) +
    theme(legend.position = "none")
  
  if (add_contours) {
    p <- p +
      geom_contour(
        aes(z = predicted_kriged_value),
        color = "gray55",
        alpha = 0.5,
        bins = contour_bins
      ) +
      metR::geom_text_contour(
        aes(
          z = predicted_kriged_value,
          label = scales::number(after_stat(level), accuracy = contour_label_accuracy)
        ),
        stroke = 0.2,
        size = 3,
        check_overlap = TRUE,
        skip = 0.5,
        min.size = 20,
        color = "black"
      )
  }
  
  if (!is.null(mask_panel) && nrow(mask_panel) > 0) {
    p <- p +
      geom_tile(
        data = mask_panel,
        aes(x = rel_dist, y = depth_kriged_m),
        inherit.aes = FALSE,
        fill = "white",
        alpha = 0.9
      )
  }
  
  if (!is.null(bathy_poly) && nrow(bathy_poly) > 0) {
    p <- p +
      geom_polygon(
        data = bathy_poly,
        aes(x = rel_dist, y = B_smooth),
        inherit.aes = FALSE,
        fill = "gray40",
        alpha = 1
      )
  }
  
  p
}

make_env_panel <- function(df_panel, cruise, transect,
                           variable, breaks_vec, fill_name,
                           palette_vals, level_labels,
                           title_text = NULL,
                           show_x = FALSE, show_y = FALSE) {
  
  if (is.null(df_panel) || nrow(df_panel) == 0) {
    return(blank_panel())
  }
  
  max_depth <- max(df_panel$depth_kriged_m, na.rm = TRUE)
  bathy_poly <- build_bathy_poly(cruise, transect, max_depth)
  mask_panel <- get_mask_panel(mask_all, cruise, transect)
  
  p <- ggplot(df_panel, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_contour_filled(
      aes(z = predicted_kriged_value),
      breaks = breaks_vec
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    scale_x_continuous(breaks = get_rel_dist_breaks(df_panel$rel_dist), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_fill_manual(
      values = stats::setNames(palette_vals, level_labels),
      limits = level_labels,
      drop = FALSE,
      name = fill_name
    ) +
    labs(
      title = title_text,
      x = if (show_x) "Relative distance" else NULL,
      y = if (show_y) "Depth (m)" else NULL
    ) +
    panel_theme(show_x = show_x, show_y = show_y) +
    theme(legend.position = "none")
  
  if (!is.null(mask_panel) && nrow(mask_panel) > 0) {
    p <- p +
      geom_tile(
        data = mask_panel,
        aes(x = rel_dist, y = depth_kriged_m),
        inherit.aes = FALSE,
        fill = "white",
        alpha = 0.9
      )
  }
  
  if (!is.null(bathy_poly) && nrow(bathy_poly) > 0) {
    p <- p +
      geom_polygon(
        data = bathy_poly,
        aes(x = rel_dist, y = B_smooth),
        inherit.aes = FALSE,
        fill = "gray40",
        alpha = 1
      )
  }
  
  p
}

## scale variable helpers:
extract_legend <- function(p) {
  gt <- ggplotGrob(p)
  guide_idx <- which(grepl("^guide-box", gt$layout$name))
  
  if (length(guide_idx) == 0) return(NULL)
  
  # keep the first guide grob that actually contains something
  for (i in guide_idx) {
    g <- gt$grobs[[i]]
    if (!is.null(g)) return(g)
  }
  
  NULL
}

make_discrete_legend_plot <- function(labels_vec, palette_vals, fill_name,
                                      direction = "vertical") {
  lev_ids <- paste0("bin_", seq_along(labels_vec))
  
  df_leg <- data.frame(
    x = 1,
    y = seq_along(labels_vec),
    fill_bin = factor(lev_ids, levels = lev_ids)
  )
  
  ggplot(df_leg, aes(x = x, y = y, fill = fill_bin)) +
    geom_tile() +
    scale_fill_manual(
      values = stats::setNames(palette_vals, lev_ids),
      breaks = lev_ids,
      labels = labels_vec,
      drop = FALSE,
      name = fill_name,
      na.value = "white"
    ) +
    theme_void() +
    theme(
      legend.position = if (direction == "horizontal") "top" else "right",
      legend.direction = direction,
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.width = if (direction == "horizontal") unit(1.8, "cm") else unit(1.0, "cm"),
      legend.key.height = unit(0.7, "cm")
    )
}

make_env_legend_plot <- function(level_labels, palette_vals, fill_name,
                                 direction = "vertical") {
  lev_ids <- paste0("bin_", seq_along(level_labels))
  
  df_leg <- data.frame(
    x = 1,
    y = seq_along(level_labels),
    fill_bin = factor(lev_ids, levels = lev_ids)
  )
  
  ggplot(df_leg, aes(x = x, y = y, fill = fill_bin)) +
    geom_tile() +
    scale_fill_manual(
      values = stats::setNames(palette_vals, lev_ids),
      breaks = lev_ids,
      labels = level_labels,
      drop = FALSE,
      name = fill_name
    ) +
    theme_void() +
    theme(
      legend.position = if (direction == "horizontal") "top" else "right",
      legend.direction = direction,
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.width = if (direction == "horizontal") unit(1.8, "cm") else unit(1.0, "cm"),
      legend.key.height = unit(0.7, "cm")
    )
}

make_continuous_legend_plot <- function(scale_obj, direction = "vertical") {
  ggplot(data.frame(x = 1, y = 1, z = 0), aes(x, y, fill = z)) +
    geom_tile() +
    scale_obj +
    theme_void() +
    theme(
      legend.position = if (direction == "horizontal") "top" else "right",
      legend.direction = direction,
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.width = if (direction == "horizontal") unit(1.8, "cm") else unit(1.0, "cm"),
      legend.key.height = unit(0.7, "cm")
    )
}

make_slope_anomaly_legend_plot <- function() {
  ggplot(data.frame(x = 1, y = 1, z = -1), aes(x, y, fill = z)) +
    geom_point(alpha = 0) +
    colorspace::scale_fill_continuous_divergingx(
      palette = "RdBu",
      limits = c(-2, 0),
      mid = -1,
      rev = TRUE,
      oob = scales::squish,
      na.value = "white",
      name = expression(paste("Slope (", m^{-3}, " ", mu, "m"^{-3}, ")")),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(10, "cm"),
        barheight = unit(0.6, "cm")
      )
    ) +
    theme_void() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
}

# ---------------------------
# Fixed continuous scales / transforms
# ---------------------------

biovol_vars <- c("total_biovolume_sum", "integrated_biovol_all", "integrated_biovol_sample")

get_fill_value <- function(x, variable_name) {
  if (variable_name %in% biovol_vars) {
    return(log10(x))
  }
  
  if (variable_name == "intercept") {
    return(x)
  }
  
  if (variable_name == "chl_a_ug_l") {
    x <- ifelse(x > 0, x, NA_real_)
    return(log10(x))
  }
  
  x
}

get_var_scale <- function(variable_name, legend_direction = "vertical") {
  
  guide_obj <- if (legend_direction == "horizontal") {
    guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.6, "cm")
    )
  } else {
    guide_colorbar(
      direction = "vertical",
      title.position = "top",
      barwidth = unit(0.7, "cm"),
      barheight = unit(6.5, "cm")
    )
  }
  
  if (variable_name %in% biovol_vars) {
    return(
      scale_fill_viridis_c(
        option = "D",
        direction = 1,
        limits = c(10, 17),                  # log10(1e5) to log10(1e25)
        oob = scales::squish,
        na.value = "white",
        breaks = seq(10, 17, by = 3),
        labels = function(x) scales::label_scientific(digits = 1)(10^x),
        name = expression(paste("Biovolume (", mu, "m"^3, " ", m^{-3}, ")")),
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "slope") {
    return(
      scale_fill_viridis_c(
        option = "D",
        direction = 1,
        limits = c(-2, 0),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(-2, 0, by = 0.25),
        name = expression(paste("Slope (", m^{-3}, " ", mu, "m"^{-3}, ")")),
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "slope_anomaly") {
    return(
      colorspace::scale_fill_continuous_divergingx(
        palette = "RdBu",
        limits = c(-2, 0),
        mid = -1,
        rev = TRUE,
        oob = scales::squish,
        na.value = "white",
        name = expression(paste("Slope (", m^{-3}, " ", mu, "m"^{-3}, ")")),
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          barwidth = unit(10, "cm"),
          barheight = unit(0.6, "cm")
        )
      )
    )
  }
  
  if (variable_name == "intercept") {
    return(
      scale_fill_viridis_c(
        option = "D",
        limits = c(5, 15),                  # corresponds to 1e2 to 1e30
        oob = scales::squish,
        na.value = "white",
        breaks = seq(5, 15, by = 5),
        labels = function(x) scales::label_scientific(digits = 1)(10^x),
        name = expression(paste("Intercept (", mu, "m"^3, " ", m^{-3}, " ", mu, "m"^{-3}, ")")),
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "Temperature") {
    return(
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        limits = c(8, 18),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(8, 18, by = 2),
        name = expression(paste(" Summer Temperature (", degree*C, ")")),
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "Oxygen_ml_l") {
    return(
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        limits = c(1, 8.5),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(1, 8.5, by = 1.5),
        name = expression(paste("Oxygen (mL ", L^{-1}, ")")),
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "Salinity") {
    return(
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        limits = c(30, 34),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(30, 34, by = 0.5),
        name = "Salinity (PSU)",
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "chl_a_ug_l") {
    chl_breaks_raw <- c( 0.01, 0.1, 1, 5, 10, 24)
    
    return(
      scale_fill_viridis_c(
        option = "D",
        limits = log10(c(0.01, 24)),
        oob = scales::squish,
        na.value = "white",
        breaks = log10(chl_breaks_raw),
        labels = scales::label_number()(chl_breaks_raw),
        name = expression(paste("Chl-a (", mu, "g ", L^{-1}, ")")),
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "pH") {
    return(
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        limits = c(7.5, 8.5),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(7.5, 8.5, by = 0.25),
        name = "pH",
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "PAR") {
    return(
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        limits = c(0, 0.030),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(0, 0.030, by = 0.005),
        labels = scales::label_number(accuracy = 0.001),
        name = "PAR (unitless)",
        guide = guide_obj
      )
    )
  }
  
  if (variable_name == "biovol_diff") {
    return(
      colorspace::scale_fill_continuous_divergingx(
        palette = "RdBu",
        mid = 0,
        limits = c(-200, 200),
        rev = TRUE,
        oob = scales::squish,
        na.value = "white",
        name = "% difference",
        guide = guide_obj
      )
    )
  }
  
  stop("No scale defined for variable_name: ", variable_name)
}

get_var_scale_by_cruise <- function(variable_name, cruise, legend_direction = "vertical") {
  
  if (variable_name == "Temperature" && cruise %in% c("w22", "w23")) {
    
    guide_obj <- if (legend_direction == "horizontal") {
      guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        barwidth = unit(8, "cm"),
        barheight = unit(0.5, "cm")
      )
    } else {
      guide_colorbar(
        direction = "vertical",
        title.position = "top",
        barwidth = unit(0.7, "cm"),
        barheight = unit(6.5, "cm")
      )
    }
    
    return(
      scale_fill_distiller(
        palette = "RdYlBu",
        direction = -1,
        limits = c(7.5, 10),
        oob = scales::squish,
        na.value = "white",
        breaks = seq(7.5, 10, by = 0.5),
        name = expression(paste("Winter temperature (", degree*C, ")")),
        guide = guide_obj
      )
    )
  }
  
  get_var_scale(variable_name, legend_direction = legend_direction)
}

get_rel_dist_breaks <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(0)
  }
  
  rng <- range(x, na.rm = TRUE)
  
  pretty_breaks <- pretty(rng, n = 4)
  breaks <- sort(unique(c(pretty_breaks, 0)))
  
  breaks[breaks >= rng[1] & breaks <= rng[2]]
}

make_fixed_panel_old <- function(df_panel, cruise, transect, variable_name,
                             scale_obj,
                             title_text = NULL,
                             show_x = FALSE, show_y = FALSE,
                             add_contours = FALSE,
                             contour_breaks = NULL,
                             contour_label_accuracy = 0.1) {
  
  if (is.null(df_panel) || nrow(df_panel) == 0) {
    return(blank_panel())
  }
  
  df_plot <- df_panel %>%
    mutate(fill_value = get_fill_value(predicted_kriged_value, variable_name)) %>%
    filter(is.finite(fill_value))
  
  if (nrow(df_plot) == 0) {
    return(blank_panel())
  }
  
  max_depth <- max(df_plot$depth_kriged_m, na.rm = TRUE)
  bathy_poly <- build_bathy_poly(cruise, transect, max_depth)
  mask_panel <- get_mask_panel(mask_all, cruise, transect)
  
  p <- ggplot(df_plot, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_raster(aes(fill = fill_value)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 1.1
    ) +
    scale_x_continuous(breaks = get_rel_dist_breaks(df_plot$rel_dist), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_obj +
    labs(
      title = ifelse(is.null(title_text), paste0(cruise, " / ", transect), title_text),
      x = if (show_x) "Relative distance" else NULL,
      y = if (show_y) "Depth (m)" else NULL
    ) +
    panel_theme(show_x = show_x, show_y = show_y) +
    theme(legend.position = "none")
  
  if (add_contours) {
    p <- p +
      geom_contour(
        aes(z = fill_value),
        color = "gray55",
        alpha = 0.5,
        breaks = contour_breaks
      ) +
      metR::geom_text_contour(
        aes(
          z = fill_value,
          label = scales::number(after_stat(level), accuracy = contour_label_accuracy)
        ),
        stroke = 0.2,
        size = 3,
        check_overlap = TRUE,
        skip = 0.5,
        min.size = 20,
        color = "black"
      )
  }
  
  if (!is.null(mask_panel) && nrow(mask_panel) > 0) {
    p <- p +
      geom_tile(
        data = mask_panel,
        aes(x = rel_dist, y = depth_kriged_m),
        inherit.aes = FALSE,
        fill = "white",
        alpha = 0.9
      )
  }
  
  if (!is.null(bathy_poly) && nrow(bathy_poly) > 0) {
    p <- p +
      geom_polygon(
        data = bathy_poly,
        aes(x = rel_dist, y = B_smooth),
        inherit.aes = FALSE,
        fill = "gray40",
        alpha = 1
      )
  }
  
  p
}



#new ver of make_fixed_panel
make_fixed_panel <- function(df_panel, cruise, transect, variable_name,
                             scale_obj,
                             title_text = NULL,
                             show_x = FALSE, show_y = FALSE,
                             add_contours = FALSE,
                             contour_breaks = NULL,
                             contour_label_accuracy = 0.1) {
  
  if (is.null(df_panel) || nrow(df_panel) == 0) {
    return(blank_panel())
  }
  
  df_plot <- df_panel %>%
    dplyr::mutate(fill_value = get_fill_value(predicted_kriged_value, variable_name)) %>%
    dplyr::filter(is.finite(fill_value))
  
  if (nrow(df_plot) == 0) {
    return(blank_panel())
  }
  
  max_depth <- max(df_plot$depth_kriged_m, na.rm = TRUE)
  bathy_poly <- build_bathy_poly(cruise, transect, max_depth)
  mask_panel <- get_mask_panel(mask_all, cruise, transect)
  
  p <- ggplot(df_plot, aes(x = rel_dist, y = depth_kriged_m)) +
    geom_raster(aes(fill = fill_value), interpolate = TRUE) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 1.1
    ) +
    scale_x_continuous(
      breaks = get_rel_dist_breaks(df_plot$rel_dist),
      expand = c(0, 0)
    ) +
    scale_y_reverse(expand = c(0, 0)) +
    scale_obj +
    labs(
      title = ifelse(is.null(title_text), paste0(cruise, " / ", transect), title_text),
      x = if (show_x) "Relative distance" else NULL,
      y = if (show_y) "Depth (m)" else NULL
    ) +
    panel_theme(show_x = show_x, show_y = show_y) +
    theme(legend.position = "none")
  
  if (add_contours) {
    p <- p +
      geom_contour(
        aes(z = fill_value),
        color = "gray55",
        alpha = 0.5,
        breaks = contour_breaks
      ) +
      metR::geom_text_contour(
        aes(
          z = fill_value,
          label = scales::number(after_stat(level), accuracy = contour_label_accuracy)
        ),
        stroke = 0.2,
        size = 3,
        check_overlap = TRUE,
        skip = 0.5,
        min.size = 20,
        color = "black"
      )
  }
  
  if (!is.null(mask_panel) && nrow(mask_panel) > 0) {
    p <- p +
      geom_tile(
        data = mask_panel,
        aes(x = rel_dist, y = depth_kriged_m),
        inherit.aes = FALSE,
        fill = "white",
        alpha = 0.9
      )
  }
  
  if (!is.null(bathy_poly) && nrow(bathy_poly) > 0) {
    p <- p +
      geom_polygon(
        data = bathy_poly,
        aes(x = rel_dist, y = B_smooth),
        inherit.aes = FALSE,
        fill = "gray40",
        alpha = 1
      )
  }
  
  p
}

make_panel_letters <- function(df, cruise_levels, transect_levels) {
  df %>%
    dplyr::distinct(Cruise, Transect_ID) %>%
    dplyr::mutate(
      Cruise = factor(as.character(Cruise), levels = cruise_levels),
      Transect_ID = factor(as.character(Transect_ID), levels = transect_levels)
    ) %>%
    dplyr::arrange(Transect_ID, Cruise) %>%
    dplyr::mutate(
      panel_letter = paste0(letters[dplyr::row_number()], ".")
    )
}

# ---------------------------
# Build one full 6 x 4 figure by variable
# ---------------------------
build_variable_grid <- function(variable_name) {
  
  legend_direction <- "horizontal"
  
  legend_grob <- NULL
  
  if (variable_name != "Temperature") {
    scale_obj <- get_var_scale(variable_name, legend_direction = legend_direction)
    
    legend_grob <- extract_legend(
      make_continuous_legend_plot(
        scale_obj = scale_obj,
        direction = legend_direction
      )
    )
  }
  
  plots <- list()
  
  df_for_letters <- kriged_all %>%
    dplyr::filter(
      variable == variable_name,
      is.finite(predicted_kriged_value),
      is.finite(rel_dist),
      is.finite(depth_kriged_m)
    )
  
  # Do not plot Summer 2023 oxygen panels
  if (variable_name == "Oxygen_ml_l") {
    df_for_letters <- df_for_letters %>%
      dplyr::filter(Cruise != "s23")
  }
  
  panel_letter_df <- df_for_letters %>%
    make_panel_letters(
      cruise_levels = cruise_levels,
      transect_levels = transect_levels
    )
  
  source_var_for_letters <- variable_name
  
  panel_letter_df <- kriged_all %>%
    dplyr::filter(
      variable == source_var_for_letters,
      is.finite(predicted_kriged_value),
      is.finite(depth_kriged_m),
      is.finite(rel_dist)
    ) %>%
    make_panel_letters(
      cruise_levels = cruise_levels,
      transect_levels = transect_levels
    )
  
  for (tr in transect_levels) {
    for (cr in cruise_levels) {
      
      
      df_panel <- kriged_all %>%
        dplyr::filter(
          variable == variable_name,
          Cruise == cr,
          Transect_ID == tr
        ) %>%
        drop_na(predicted_kriged_value, depth_kriged_m, rel_dist)
      
      # Skip Summer 2023 oxygen panels
      if (variable_name == "Oxygen_ml_l" && cr == "s23") {
        df_panel <- df_panel[0, ]
      }
      
      show_x <- TRUE
      show_y <- cr == head(cruise_levels, 1)
      
      panel_letter <- panel_letter_df %>%
        dplyr::filter(
          as.character(Cruise) == cr,
          as.character(Transect_ID) == tr
        ) %>%
        dplyr::pull(panel_letter)
      
      if (length(panel_letter) == 0) {
        panel_letter <- NULL
      }
      
      scale_obj_panel <- get_var_scale_by_cruise(
        variable_name = variable_name,
        cruise = cr,
        legend_direction = legend_direction
      )
      
      contour_breaks_panel <- if (variable_name == "slope") {
        seq(-2, 0, by = 0.25)
      } else if (variable_name %in% c("integrated_biovol_all", "integrated_biovol_sample")) {
        # fill_value is log10(biovolume), so contour levels are also in log10 space
        seq(10, 17, by = 1)
      } else {
        NULL
      }
      
      p <- make_fixed_panel(
        df_panel = df_panel,
        cruise = cr,
        transect = tr,
        variable_name = variable_name,
        scale_obj = scale_obj_panel,
        title_text = panel_letter,
        show_x = show_x,
        show_y = show_y,
        add_contours = TRUE,
        contour_breaks = contour_breaks_panel,
        contour_label_accuracy = if (variable_name == "slope") 0.1 else 1
      )
      
      plots[[paste(tr, cr, sep = "_")]] <- p
    }
  }
  
  header_row <- c(
    lapply(cruise_levels, function(cr) column_label_plot(cruise_label_map[[cr]])),
    list(ggplot() + theme_void())
  )
  
  all_rows <- list(header_row)
  
  for (tr in transect_levels) {
    row_plots <- lapply(cruise_levels, function(cr) {
      plots[[paste(tr, cr, sep = "_")]]
    })
    row_plots <- c(row_plots, list(row_label_plot(transect_label_map[[tr]])))
    all_rows[[length(all_rows) + 1]] <- row_plots
  }
  
  plot_list <- unlist(all_rows, recursive = FALSE)
  
  main_grid <- wrap_plots(
    plotlist = plot_list,
    ncol = length(cruise_levels) + 1,
    byrow = TRUE,
    widths = c(rep(1, length(cruise_levels)), 0.12),
    heights = c(0.10, rep(1, length(transect_levels)))
  ) &
    theme(panel.spacing.y = unit(0.05, "lines"))
  
  # For-press layout: use the same top-horizontal legend configuration
  # for slope and both integrated-biovolume figures.
  combined <- cowplot::plot_grid(
    legend_grob,
    main_grid,
    ncol = 1,
    rel_heights = c(0.05, 0.95)
  )
  
  final_plot <- cowplot::ggdraw() +
    cowplot::draw_plot(combined, x = 0.03, y = 0.04, width = 0.94, height = 0.925) +
    cowplot::draw_label(
      "Depth (m)",
      x = 0.025, y = 0.5,
      angle = 90,
      fontface = "bold",
      size = 18
    ) +
    cowplot::draw_label(
      "Relative distance",
      x = 0.5, y = 0.035,
      fontface = "bold",
      size = 18
    )
  
  final_plot
}


# ---------------------------
# Build and save ONLY the three for-press figures
# ---------------------------
var_file_stub <- c(
  slope = "multi_panel_plankton_slope_for_press",
  integrated_biovol_all = "multi_panel_plankton_biovol_integratedAll_for_press",
  integrated_biovol_sample = "multi_panel_plankton_biovol_integratedSample_for_press"
)

for (v in target_vars) {
  message("Building for-press figure for variable: ", v)
  
  p_grid <- build_variable_grid(v)
  
  out_png <- file.path(output_dir, paste0(var_file_stub[[v]], ".png"))
  out_pdf <- file.path(output_dir, paste0(var_file_stub[[v]], ".pdf"))
  
  ggsave(
    filename = out_png,
    plot = p_grid,
    width = 16,
    height = 18,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = out_pdf,
    plot = p_grid,
    width = 16,
    height = 18,
    units = "in",
    bg = "white"
  )
}

message("Done. Saved only slope and integrated-biovolume for-press figures.")
