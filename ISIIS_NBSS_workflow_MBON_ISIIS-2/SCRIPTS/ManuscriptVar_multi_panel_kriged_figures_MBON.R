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
# dataSource:
#   plankton     -> ISIIS
#   allParticles -> LISST-ISIIS
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
option_list <- list(
  make_option(
    c("-d", "--dataSource"),
    type = "character",
    default = "plankton",
    help = "Data source: 'plankton' or 'allParticles'"
  ),
  make_option(
    c("-o", "--outputDir"),
    type = "character",
    default = NULL,
    help = "Optional output directory"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!opt$dataSource %in% c("plankton", "allParticles", "LISST_only")) {
  stop("--dataSource must be 'plankton', 'allParticles', or 'LISST_only'")
}

instrument_dir <- if (opt$dataSource == "plankton") {
  "ISIIS"
} else if (opt$dataSource == "allParticles") {
  "LISST-ISIIS"
} else if (opt$dataSource == "LISST_only") {
  "LISST"
} else {
  stop("Unknown dataSource option")
}

input_dir <- if (opt$dataSource == "plankton") {
  file.path(base_dir, "data", "interp_and_krige_dist10", "ISIIS", "plankton", "size_spectra_coeffs")
} else if (opt$dataSource == "allParticles") {
  file.path(base_dir, "data", "interp_and_krige_dist10", "LISST-ISIIS", "size_spectra_coeffs")
} else if (opt$dataSource == "LISST_only") {
  file.path(base_dir, "data", "interp_and_krige_dist10", "LISST", "size_spectra_coeffs")
}

input_dir_abundance <- file.path(base_dir, "data", "interp_and_krige_dist10", "abundance_data")

output_dir <- if (is.null(opt$outputDir)) {
  file.path(
    base_dir, "figures", "interp_and_krige_dist10",
    "multi_panel_kriged_preserved_panels", instrument_dir
  )
} else {
  opt$outputDir
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading kriged files from: ", input_dir)
message("Writing figures to: ", output_dir)

# ---------------------------
# Constants
# ---------------------------
transect_levels <- c("LP", "GH","NH", "HH", "CC")

panel_layout <- tibble(
    Cruise=c(
        "f20","f20",
        "sp21","sp21"
    ),
    time_of_day=c(
        "day","night",
        "day","night"
    ),
    label=c(
        "Fall 2020\nDay",
        "Fall 2020\nNight",
        "Spring 2021\nDay",
        "Spring 2021\nNight"
    )
)

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

target_vars <- c(metric_vars, env_vars)

transect_label_map <- c(
    LP="La Push",
    GH="Grays Harbor",
    NH="Newport Line",
    HH="Heceta Head",
    CC="Crescent City"
)

fmt_sci <- scales::label_scientific(digits = 2)

# ---------------------------
# Read all kriged outputs, size spectra and abundance
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

#abundance
csv_files_abundance <- list.files(
  input_dir_abundance,
  pattern = "plankton-abundance_envs_interp_and_krige_dist10\\.csv$",
  full.names = TRUE,
  recursive = TRUE
) 
  
if (length(csv_files_abundance) == 0) {
  stop(
    "No abundance kriged CSV files found in: ",
    input_dir_abundance
  )
}

kriged_all_abundance <- purrr::map_dfr(csv_files_abundance, function(f) {
  readr::read_csv(f, show_col_types = FALSE)
})
  
#assing day night labels:

tow_lookup <- readr::read_csv(file.path(base_dir, "ancillary", "towId_timeOfDay.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(Cruise = as.character(Cruise), tow_id = as.character(tow_id), time_of_day = tolower(trimws(as.character(time_of_day))) ) %>%
  dplyr::distinct(Cruise, tow_id, .keep_all = TRUE)

kriged_all <- kriged_all %>% 
    dplyr::mutate(Cruise = as.character(Cruise), tow_id = as.character(tow_id)) %>%
    dplyr::left_join(tow_lookup, by = c("Cruise", "tow_id"))

kriged_all_abundance <- kriged_all_abundance %>%
    dplyr::mutate(Cruise = as.character(Cruise), tow_id = as.character(tow_id)) %>%
    dplyr::left_join(tow_lookup, by = c("Cruise", "tow_id"))


#Missing panel check
missing_abundance_panels <- tidyr::crossing(
  Transect_ID = transect_levels,
  panel_layout %>%
    dplyr::select(Cruise, time_of_day)
) %>%
  dplyr::anti_join(
    kriged_all_abundance %>%
      dplyr::distinct(
        Cruise,
        Transect_ID,
        time_of_day
      ),
    by = c(
      "Cruise",
      "Transect_ID",
      "time_of_day"
    )
  )

if (nrow(missing_abundance_panels) > 0) {
  message("Abundance panels with no data:")
  print(missing_abundance_panels, n = Inf)
}


#read the no data masks
mask_files <- list.files(
  input_dir,
  pattern = "no_data_mask_interp_and_krige_dist10\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

message(
  "Size-spectrum rows missing time_of_day: ",
  sum(is.na(kriged_all$time_of_day))
)

message(
  "Abundance rows missing time_of_day: ",
  sum(is.na(kriged_all_abundance$time_of_day))
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
      Cruise %in% unique(panel_layout$Cruise),
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

#clean up data
missing_cols <- setdiff(required_cols, names(kriged_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns in kriged data: ", paste(missing_cols, collapse = ", "))
}

kriged_all <- kriged_all %>%
  mutate(
    Cruise = as.character(Cruise),
    Transect_ID = as.character(Transect_ID),
    tow_id = as.character(tow_id),
    time_of_day = tolower(trimws(as.character(time_of_day))),
    rel_dist = as.numeric(rel_dist),
    depth_kriged_m = as.numeric(depth_kriged_m),
    predicted_kriged_value = as.numeric(predicted_kriged_value)
  ) %>%
  filter(
    variable %in% target_vars,
    Cruise %in% unique(panel_layout$Cruise),
    Transect_ID %in% transect_levels,
    time_of_day %in% c("day", "night")
  )


missing_cols_abundance <- setdiff(
  c(required_cols, "tow_id", "time_of_day"),
  names(kriged_all_abundance)
)

if (length(missing_cols_abundance) > 0) {
  stop(
    "Missing required columns in abundance kriged data: ",
    paste(missing_cols_abundance, collapse = ", ")
  )
}

kriged_all_abundance <- kriged_all_abundance %>%
  dplyr::mutate(
    Cruise = as.character(Cruise),
    Transect_ID = as.character(Transect_ID),
    tow_id = as.character(tow_id),
    time_of_day = tolower(trimws(as.character(time_of_day))),
    rel_dist = as.numeric(rel_dist),
    depth_kriged_m = as.numeric(depth_kriged_m),
    predicted_kriged_value = as.numeric(predicted_kriged_value),
    variable = as.character(variable)
  ) %>%
  dplyr::filter(
    variable == "ind_m3_total",
    Cruise %in% unique(panel_layout$Cruise),
    Transect_ID %in% transect_levels,
    time_of_day %in% c("day", "night")
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
      plot.title = element_text(size = 18, hjust = 0, face = "bold"),
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
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 5.5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(12, 2, 12, 2)
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
 
  if (variable_name == "ind_m3_total") {
    x <- ifelse(x > 0, x, NA_real_)
    return(log10(x))
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
  
  if (variable_name %in% biovol_vars) {
    return(
      scale_fill_viridis_c(
        option = "D",
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
        name = expression(paste("Temperature (", degree*C, ")")),
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

  if (variable_name == "ind_m3_total") {

  abundance_breaks <- c(1,10,100,1000,10000,100000)

  return(
    scale_fill_viridis_c(
      option = "D",
      limits = log10(range(abundance_breaks)),
      oob = scales::squish,
      na.value = "white",

      breaks = log10(abundance_breaks),
      labels = scales::label_number()(abundance_breaks),

      name = expression(paste("Abundance (ind ",m^{-3},")")),

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
build_variable_grid <- function(variable_name, kriged_dataset) {

  legend_direction <- if (
    variable_name == "slope_anomaly"
  ) {
    "horizontal"
  } else {
    "vertical"
  }

  scale_obj <- get_var_scale(
    variable_name,
    legend_direction = legend_direction
  )

  legend_grob <- extract_legend(
    make_continuous_legend_plot(
      scale_obj = scale_obj,
      direction = legend_direction
    )
  )

  plots <- list()

  source_variable <- if (
    variable_name == "slope_anomaly"
  ) {
    "slope"
  } else {
    variable_name
  }

  for (tr in transect_levels) {
    for (i in seq_len(nrow(panel_layout))) {

      cr  <- panel_layout$Cruise[i]
      tod <- panel_layout$time_of_day[i]

      df_panel <- kriged_dataset %>%
        dplyr::filter(
          variable == source_variable,
          Cruise == cr,
          Transect_ID == tr,
          time_of_day == tod
        ) %>%
        tidyr::drop_na(
          predicted_kriged_value,
          depth_kriged_m,
          rel_dist
        )
      if (cr == "sp21" && tr == "CC") {
         df_panel <- df_panel[0, ]
         }

      show_x <- TRUE
      show_y <- i == 1

      scale_obj_panel <- get_var_scale(
        variable_name,
        legend_direction = legend_direction
      )

      if (nrow(df_panel) > 0) {
            cat(
              variable_name,
              cr,
              tr,
              tod,
              range(df_panel$predicted_kriged_value, na.rm = TRUE),
              "\n"
            )
          } else {
            cat(
              variable_name,
              cr,
              tr,
              tod,
              "NO DATA\n"
            )
          }

      p <- make_fixed_panel(
        df_panel = df_panel,
        cruise = cr,
        transect = tr,
        variable_name = variable_name,
        scale_obj = scale_obj_panel,
        title_text = NULL,
        show_x = show_x,
        show_y = show_y,
        add_contours = variable_name == "slope_anomaly",
        contour_breaks = if (
          variable_name == "slope_anomaly"
        ) {
          seq(-2.5, 0, by = 0.25)
        } else {
          NULL
        },
        contour_label_accuracy = 0.1
      )

      plots[[paste(tr, cr, tod, sep = "_")]] <- p
    }
  }

  header_row <- c(
    lapply(panel_layout$label, column_label_plot),
    list(ggplot() + theme_void())
  )

  all_rows <- list(header_row)

  for (tr in transect_levels) {
    row_plots <- lapply(
      seq_len(nrow(panel_layout)),
      function(i) {
        plots[[
          paste(
            tr,
            panel_layout$Cruise[i],
            panel_layout$time_of_day[i],
            sep = "_"
          )
        ]]
      }
    )

    row_plots <- c(
      row_plots,
      list(row_label_plot(transect_label_map[[tr]]))
    )

    all_rows[[length(all_rows) + 1]] <- row_plots
  }

  plot_list <- unlist(all_rows, recursive = FALSE)

  n_panel_cols <- nrow(panel_layout)

  main_grid <- patchwork::wrap_plots(
    plotlist = plot_list,
    ncol = n_panel_cols + 1,
    byrow = TRUE,
    widths = c(rep(1, n_panel_cols), 0.12),
    heights = c(
      0.25,
      rep(1, length(transect_levels))
    )
  )

  if (variable_name == "slope_anomaly") {

    legend_plot <- make_slope_anomaly_legend_plot()

    combined <- cowplot::plot_grid(
      legend_plot,
      main_grid,
      ncol = 1,
      rel_heights = c(0.05, 0.95)
    )

  } else {

    combined <- cowplot::plot_grid(
      main_grid,
      legend_grob,
      nrow = 1,
      rel_widths = c(1, 0.18),
      align = "h"
    )
  }

  cowplot::ggdraw() +
    cowplot::draw_plot(
      combined,
      x = 0.03,
      y = 0.04,
      width = 0.94,
      height = 0.94
    ) +
    cowplot::draw_label(
      "Depth (m)",
      x = 0.018,
      y = 0.5,
      angle = 90,
      fontface = "bold",
      size = 18
    ) +
    cowplot::draw_label(
      "Relative distance",
      x = 0.5,
      y = 0.045,
      fontface = "bold",
      size = 18
    )
}

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
#  # Do not plot Summer 2023 oxygen panels
#  
#  
#  #
#  
#        
#      
#      
#     
#                   variable_name
#                 rel_dist
#
#
# 
#      # Skip Summer 2023 oxygen panels
#      
#
#      #
#      
#
#      
#      ##DEBUG slope anomaly
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
#
# 
#    
#    
#    
#    
#    
#      
# 
#      #
#      #
#      #
#      
#      
#    
#  
#  final_plot


# ---------------------------
# Build and save all figures
# ---------------------------
var_file_stub <- c(
  slope = "multi_panel_plankton_slope",
  slope_anomaly = "multi_panel_plankton_slope_anomaly",
  intercept = "multi_panel_plankton_intercept",
  total_biovolume_sum = "multi_panel_plankton_biovol_totalSum",
  integrated_biovol_all = "multi_panel_plankton_biovol_integratedAll",
  integrated_biovol_sample = "multi_panel_plankton_biovol_integratedSample",
  biovol_diff = "multi_panel_plankton_biovol_diff",
  Temperature = "multi_panel_plankton_temperature",
  Salinity = "multi_panel_plankton_salinity",
  Oxygen_ml_l = "multi_panel_plankton_oxygen",
  PAR = "multi_panel_plankton_PAR",
  pH = "multi_panel_plankton_pH",
  chl_a_ug_l = "multi_panel_plankton_chla"
)

if (opt$dataSource == "allParticles") {
  var_file_stub <- sub("plankton", "allParticles", var_file_stub)
}

for (v in target_vars) {
  message("Building figure for variable: ", v)
  
  p_grid <- build_variable_grid(variable_name=v, kriged_dataset=kriged_all)
  
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


message("Building abundance figure")

p_grid <- build_variable_grid(
    variable_name = "ind_m3_total",
    kriged_dataset = kriged_all_abundance
)

ggsave(
    file.path(output_dir,
              "multi_panel_plankton_abundance.png"),
    p_grid,
    width = 16,
    height = 18,
    dpi = 300
)

ggsave(
    file.path(output_dir,
              "multi_panel_plankton_abundance.pdf"),
    p_grid,
    width = 16,
    height = 18
)


message("Done.")




