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
library(emmeans)
library(multcomp)
library(rlang)
base_dir<- getwd()

source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


option_list <- list(
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)


size_var<- opt$sizeVar # size_var<-'contour'
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "ancillary", "transect_shelf_land.csv")
sizeSpectraDirPlankton <- file.path(base_dir,  'data', 'size_spectra_plankton_5-10m-bins_thresholded')
sizeSpectraDirAll <- file.path(base_dir,  'data', 'ISIIS-LISST_mergedNBSS_5m_allParticles')
sizeSpectraDirAllISIIS <- file.path(base_dir,  'data', 'size_spectra_allParticles_5-10m-bins_thresholded')
sizeSpectraDirDetritus <- file.path(base_dir,  'data', 'size_spectra_detritus_5-10m-bins_thresholded')

LISSTdir <-  file.path(base_dir,  'data', "LISST_data", "LISST_size_spectra_thresholded", "5m_depthBin") 


output_dir_figs_NB <- file.path(base_dir, "figures", 'size_spectra')
if (!dir.exists(output_dir_figs_NB)) {
  dir.create(output_dir_figs_NB, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_NB)
}


df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = "by_taxa\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_full <- list.files(sizeSpectraDirAll, pattern = "LISST-ISIIS_NBSS\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_full_ISIIS <- list.files(sizeSpectraDirAllISIIS, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_detritus <- list.files(sizeSpectraDirDetritus, pattern = "by_taxa\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_LISST <- list.files(LISSTdir, pattern = "5m-depthBin_LISST_size_spectra\\.csv$", full.names = TRUE, recursive = TRUE)


load_ISIIS_data <- function(df_list,
                            env_vars = c("Temperature","Salinity","Oxygen_ml_l","pH","PAR","chl_a_ug_l"),
                            x1_full = 12.78849,
                            x2_full = 4.547442e13) {
  
  dat_ISIIS <- map_dfr(
    df_list,
    ~ read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        "group_id",
        dplyr::any_of(c("total_biovolume", "biovolume_um3")),
        "biovol_um3_mid", "biovol_um3_width",
        "time_Pacific", "rel_dist", "depth_mid", "Latitude", "Longitude",
        "slope", "intercept", "Cruise", "Transect_ID",
        "logNB", "ECD_um_mid",
        dplyr::any_of("grouped_taxon"),
        dplyr::any_of(env_vars)
      )
    )
  )
  
  # make optional columns exist even if absent in every file
  if (!"total_biovolume" %in% names(dat_ISIIS)) {
    dat_ISIIS$total_biovolume <- NA_real_
  }
  if (!"biovolume_um3" %in% names(dat_ISIIS)) {
    dat_ISIIS$biovolume_um3 <- NA_real_
  }
  if (!"grouped_taxon" %in% names(dat_ISIIS)) {
    dat_ISIIS$grouped_taxon <- NA_character_
  }
  
  # standardize biovolume name
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(
      total_biovolume = dplyr::coalesce(total_biovolume, biovolume_um3)
    ) %>%
    dplyr::select(-dplyr::any_of("biovolume_um3")) %>%
    dplyr::mutate(
      time_Pacific_orig = time_Pacific,
      mean_time_Pacific = time_Pacific - lubridate::hours(1)
    ) %>%
    dplyr::arrange(time_Pacific)
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::rename(group_id_old = group_id) %>%
    dplyr::mutate(
      group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
      total_biovolume  = as.numeric(total_biovolume),
      biovol_um3_mid   = as.numeric(biovol_um3_mid),
      biovol_um3_width = as.numeric(biovol_um3_width),
      slope            = as.numeric(slope),
      intercept        = as.numeric(intercept),
      lower_size_edge  = biovol_um3_mid - biovol_um3_width / 2,
      upper_size_edge  = biovol_um3_mid + biovol_um3_width / 2,
      b_lin            = 10^intercept
    )
  
  group_summ <- dat_ISIIS %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(
      total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
      min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
      max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
      slope_g = dplyr::first(slope),
      intercept_g = dplyr::first(intercept),
      b_lin_g = dplyr::first(b_lin),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      integrated_biovol_all = dplyr::if_else(
        is.finite(slope_g) & is.finite(b_lin_g),
        ifelse(
          abs(slope_g + 1) < 1e-12,
          b_lin_g * log(x2_full / x1_full),
          (b_lin_g / (slope_g + 1)) * (x2_full^(slope_g + 1) - x1_full^(slope_g + 1))
        ),
        NA_real_
      ),
      integrated_biovol_sample = dplyr::if_else(
        is.finite(slope_g) & is.finite(b_lin_g) &
          is.finite(min_size) & is.finite(max_size) &
          min_size > 0 & max_size > 0 & max_size >= min_size,
        ifelse(
          abs(slope_g + 1) < 1e-12,
          b_lin_g * log(max_size / min_size),
          (b_lin_g / (slope_g + 1)) * (max_size^(slope_g + 1) - min_size^(slope_g + 1))
        ),
        NA_real_
      ),
      biovol_diff = dplyr::if_else(
        is.finite(integrated_biovol_sample) &
          is.finite(total_biovolume_sum) &
          total_biovolume_sum > 0,
        ((integrated_biovol_sample - total_biovolume_sum) / total_biovolume_sum) * 100,
        NA_real_
      )
    ) %>%
    dplyr::select(
      group_id,
      total_biovolume_sum,
      min_size,
      max_size,
      integrated_biovol_all,
      integrated_biovol_sample,
      biovol_diff
    )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::left_join(group_summ, by = "group_id") %>%
    dplyr::mutate(
      slope_anomaly = slope + 1,
      NB = 10^(logNB)
    )
  
  return(dat_ISIIS)
}

dat_ISIIS_full     <- load_ISIIS_data(df_list_full, env_vars = env_vars)
dat_ISIIS_plankton <- load_ISIIS_data(df_list_plankton, env_vars = env_vars)
dat_ISIIS_detritus <- load_ISIIS_data(df_list_detritus, env_vars = env_vars)


#load and format LISST data
dat_LISST <- map_dfr(df_list_LISST, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  # Extract filename only
  fname <- basename(f)
  # Split by "_"  
  parts <- strsplit(fname, "_")[[1]]
  # parts should be:
  # [1] "SPECTRA"  "s22"  "HH"  "5m-depthBin"  "LISST"  "size"  "spectra.csv"
  Cruise      <- parts[2]
  Transect_ID <- parts[3]
  df %>%
    dplyr::mutate(
      Cruise      = Cruise,
      Transect_ID = Transect_ID
    ) %>% assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude') %>%
    dplyr::mutate(
      shelf_loc_200 = case_when(
        rel_dist >= 0 ~ "on_shelf",
        rel_dist <  0 ~ "off_shelf",
        TRUE          ~ NA_character_
      ),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf"))
    )
})


#     # keep the original for reference
#     # new globally consecutive id: increments whenever group_id_old changes
#         TRUE          ~ NA_character_

clean_ISIIS_dat <- function(dat_ISIIS){
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(
      shelf_loc_200 = case_when(
        rel_dist >= 0 ~ "on_shelf",
        rel_dist <  0 ~ "off_shelf",
        TRUE          ~ NA_character_
      ),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf")),
      Cruise        = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      Transect_ID   = factor(Transect_ID,
                             levels = c("GH", "OB", "CR", "CM", "NH", "HH", "RR"))
    ) %>%
    dplyr::filter(!(Cruise == "w22" & Transect_ID == "CM" & Longitude < -124.75))
  
  return(dat_ISIIS)
}


dat_ISIIS_full<- clean_ISIIS_dat(dat_ISIIS_full)
dat_ISIIS_plankton<- clean_ISIIS_dat(dat_ISIIS_plankton)
dat_ISIIS_detritus<- clean_ISIIS_dat(dat_ISIIS_detritus)




PSS_summary_plot <- function(dat,
                             string_ID,
                             PSS_var = "slope",
                             group_col = "season",
                             output_dir = output_dir_figs_NB,
                             y_limits = NULL) {
  

  
  # ----------------------------
  # Validate inputs
  # ----------------------------
  if (!PSS_var %in% names(dat)) {
    stop("PSS_var not found in dat: ", PSS_var)
  }
  if (!group_col %in% names(dat)) {
    stop("group_col not found in dat: ", group_col)
  }
  
  yvar <- sym(PSS_var)
  gvar <- sym(group_col)
  
  df <- dat %>%
    dplyr::mutate(
      .group = factor(!!gvar),
      .y     = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    dplyr::filter(!is.na(.group), !is.na(.y))
  
  # ----------------------------
  # Mean per group (for line)
  # ----------------------------
  df_mean <- df %>%
    group_by(.group) %>%
    summarise(mean_val = mean(.y, na.rm = TRUE), .groups = "drop")
  
  # ----------------------------
  # Plot
  # ----------------------------
  p <- ggplot(df, aes(x = .group, y = .y)) +
    
    # violin
    geom_violin(
      aes(fill = .group),
      trim = FALSE,
      alpha = 0.6,
      linewidth = 0.3
    ) +
    
    # boxplot (median + IQR)
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    
    # mean points
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 21,
      size = 3,
      fill = "white",
      color = "black"
    ) +
    
    # dashed line connecting means
    geom_line(
      data = df_mean,
      aes(x = .group, y = mean_val, group = 1),
      color = "black",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    
    scale_fill_brewer(palette = "Set2") +
    
    labs(
      title = paste0("NBSS ", PSS_var, " summary — ", string_ID),
      x = group_col,
      y = PSS_var
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  print(p)
  
  # ----------------------------
  # Save
  # ----------------------------
  out_path <- file.path(
    output_dir,
    paste0("SPECTRA_", PSS_var, "_summary_by_", group_col, "_", string_ID, ".png")
  )
  
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = out_path,
    plot = p,
    width = 9,
    height = 6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  message("Saved plot to: ", out_path)
  
  return(p)
}

# Compare by season
PSS_summary_plot(dat_ISIIS_plankton, "plankton", "slope", "season", y_limits = c(-2, 0))

# Compare by shelf location
PSS_summary_plot(dat_ISIIS_plankton, "plankton", "slope", "shelf_loc_200", y_limits = c(-2, 0))


# super summary figure for seasons, top 50 meters, by shelf location:
PSS_summary_shelf_season <- function(dat,
                                     string_ID,
                                     PSS_var = "slope",
                                     output_dir = output_dir_figs_NB,
                                     y_limits = NULL,
                                     depth_cut = 50,
                                     log10_transform = FALSE) {

  
  req <- c(PSS_var, "season", "shelf_loc_200", "depth_mid")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))
  
  yvar <- sym(PSS_var)
  
  df <- dat %>%
    dplyr::filter(depth_mid < depth_cut) %>%
    dplyr::mutate(
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      .y_raw = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    dplyr::filter(!is.na(season), !is.na(shelf_loc_200), !is.na(.y_raw))
  
  # Optional log10 transform (safe)
  if (log10_transform) {
    df <- df %>%
      dplyr::mutate(.y = ifelse(.y_raw > 0, log10(.y_raw), NA_real_)) %>%
      dplyr::filter(!is.na(.y))
    y_label <- paste0("log10(", PSS_var, ")")
  } else {
    df <- df %>% dplyr::mutate(.y = .y_raw)
    y_label <- PSS_var
  }
  
  df_med <- df %>%
    dplyr::group_by(.data$shelf_loc_200, .data$season) %>%
    dplyr::summarise(
      med_val = stats::median(.data$.y, na.rm = TRUE),
      .groups = "drop"
    )
  
  season_colors <- c(
    "Winter" = "#0072B2",
    "Summer" = "#E69F00"
  )
  
  dodge_w <- 0.75
  
  shelf_levels <- levels(df$shelf_loc_200)
  x_base <- setNames(seq_along(shelf_levels), shelf_levels)
  x_offsets <- c("Winter" = -dodge_w/4, "Summer" = dodge_w/4)
  
  df_med2 <- df_med %>%
    dplyr::mutate(x = x_base[as.character(shelf_loc_200)] + x_offsets[as.character(season)])
  
  p <- ggplot(df, aes(x = shelf_loc_200, y = .y)) +
    
    geom_violin(
      aes(fill = season),
      position = position_dodge(width = dodge_w),
      trim = FALSE,
      alpha = 0.65,
      linewidth = 0.3
    ) +
    
    geom_boxplot(
      aes(group = interaction(shelf_loc_200, season)),
      position = position_dodge(width = dodge_w),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    
    # median points
    geom_point(
      data = df_med2,
      aes(x = x, y = med_val),
      shape = 21, size = 3, fill = "white", color = "black",
      inherit.aes = FALSE
    ) +
    
    # dashed line connecting medians (Winter -> Summer within each shelf_loc_200)
    geom_line(
      data = df_med2,
      aes(x = x, y = med_val, group = shelf_loc_200),
      linetype = "dashed", linewidth = 0.9, color = "black",
      inherit.aes = FALSE
    ) +
    
    scale_fill_manual(values = season_colors) +
    
    labs(
      title = paste0("NBSS ", y_label, " — depth < ", depth_cut, " m — ", string_ID),
      x = "Shelf location (200 m)",
      y = y_label_expr,
      fill = "Season"
    ) +
    
    theme(
      # Axis titles
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      
      # Tick labels (numbers / categories)
      axis.text.x  = element_text(size = 16),
      axis.text.y  = element_text(size = 16),
      
      # Tick marks themselves
      axis.ticks.length = unit(6, "pt"),
      axis.ticks = element_line(linewidth = 1.2),
      # Legend title
      legend.title = element_text(size = 18, face = "bold"),
      
      # Legend labels
      legend.text  = element_text(size = 16),
      
      # Legend keys (color boxes / lines)
      legend.key.size = unit(1.4, "lines"),
      
      # Legend spacing
      legend.spacing.x = unit(0.6, "cm"),
      legend.spacing.y = unit(0.4, "cm"),
      
      # Optional: keep legend on top and horizontal
      legend.position = "top",
      legend.direction = "horizontal"
    )
  
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  
  print(p)
  
  out_path <- file.path(
    output_dir,
    paste0("SPECTRA_", PSS_var,
           if (log10_transform) "_log10" else "",
           "_violin_shelf_season_depthlt", depth_cut, "_", string_ID, ".png")
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(out_path, p, width = 9, height = 6, dpi = 300, bg = "white")
  message("Saved plot to: ", out_path)
  
  return(p)
}


# Example
PSS_summary_shelf_season(dat_ISIIS_plankton, "plankton", PSS_var = "slope", y_limits = c(-2, 0))

PSS_summary_shelf_season(
  dat_ISIIS_plankton,
  "plankton",
  PSS_var = "integrated_biovol",
  log10_transform = TRUE, y_limits = c(log10(min(dat_ISIIS_plankton$integrated_biovol)), 17))


PSS_summary_shelf_season(
  dat_ISIIS_plankton,
  "plankton",
  PSS_var = "intercept",
   y_limits = c(min(dat_ISIIS_plankton$intercept), max(dat_ISIIS_plankton$intercept)))

#accompanying mixed model to control for spatiotemporal correlation and heteroskedasticity

#ambitious function to get graph and do mixed model at the same time:
PSS_violin_plus_mixedModel <- function(dat,
                                       string_ID,
                                       PSS_var = "slope",
                                       depth_cut = 50,
                                       log10_transform = FALSE,
                                       output_dir = output_dir_figs_NB,
                                       y_limits = NULL,
                                       # mixed model controls
                                       random_effect = "(1 | Cruise/Transect_ID)",
                                       include_covariates = TRUE,
                                       # plotting controls
                                       season_colors = c("Winter" = "#0072B2",
                                                         "Summer" = "#E69F00"),
                                       gap = 0.35,
                                       ci_cap_width = 0.03,
                                       width = 12,
                                       height = 6,
                                       dpi = 300,
                                       save_plot = TRUE,
                                       # save Mixed models results arguments
                                       save_model   = TRUE,
                                       model_dir    = output_dir,   # or output_dir_figs_NB
                                       model_prefix = NULL          # default builds from string_ID + PSS_var
                                       ) {
  
  # packages (explicit namespaces everywhere)
  req_pkgs <- c("dplyr","ggplot2","rlang","lme4","lmerTest","emmeans","patchwork","tidyr","grid")
  miss_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_pkgs) > 0) stop("Missing packages: ", paste(miss_pkgs, collapse = ", "))
  
  # ----------------------------
  # Validate required columns
  # ----------------------------
  needed <- c(PSS_var, "season", "shelf_loc_200", "depth_mid",
              "Cruise","Transect_ID","group_id","distance_sample_km")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))
  
  yvar <- rlang::sym(PSS_var)
  
  # ----------------------------
  # Filter + coerce outcome
  # ----------------------------
  df <- dat %>%
    dplyr::filter(.data$depth_mid < depth_cut) %>%
    dplyr::mutate(
      season = factor(.data$season, levels = c("Winter","Summer")),
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf","on_shelf")),
      Cruise = factor(.data$Cruise),
      Transect_ID = factor(.data$Transect_ID),
      .y_raw = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    dplyr::filter(!is.na(.data$season), !is.na(.data$shelf_loc_200), is.finite(.data$.y_raw))
  
  if (nrow(df) == 0) stop("After filtering depth/NA/non-finite, no rows remain.")
  
  # optional log10 transform
  if (log10_transform) {
    df <- df %>%
      dplyr::mutate(.y = dplyr::if_else(.data$.y_raw > 0, log10(.data$.y_raw), NA_real_)) %>%
      dplyr::filter(is.finite(.data$.y))
    y_label <- paste0("log10(", PSS_var, ")")
  } else {
    df <- df %>% dplyr::mutate(.y = .data$.y_raw)
    y_label <- PSS_var
  }
  
  # -------------------------------------------------
  # Custom y-axis label expressions by PSS_var
  # -------------------------------------------------
  if (PSS_var == "integrated_biovol") {
    
    y_label_expr <- expression(
      Integrated~biovolume~(mu*m^3~m^-3)
    )
    
  } else if (PSS_var == "intercept") {
    
    y_label_expr <- expression(
      Intercept~(mu*m^3~m^-3~mu*m^-3)
    )
    
  } else if (PSS_var == "slope") {
    
    y_label_expr <- expression(
      slope~(m^-3~mu*m^-3)
    )
    
  } else {
    
    # fallback for any other variable
    y_label_expr <- PSS_var
    
  }
  
  # -------------------------------------------------
  # Conditional y-axis formatting (display as 10^x)
  # -------------------------------------------------
  if (PSS_var %in% c("intercept", "integrated_biovol")) {
    
    y_scale_obj <- ggplot2::scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", round(x,2)))
    )
    
  } else {
    
    y_scale_obj <- NULL
    
  }
  
  
  # ----------------------------
  # Build df_mod at group_id level (your approach)
  # ----------------------------
  df_mod <- df %>%
    dplyr::group_by(.data$Cruise, .data$Transect_ID, .data$group_id, .data$season, .data$shelf_loc_200) %>%
    dplyr::summarise(
      y = dplyr::first(.data$.y),
      depth_mid = dplyr::first(.data$depth_mid),
      distance_sample_km = dplyr::first(.data$distance_sample_km),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(.data$y))
  
  # ----------------------------
  # Fit mixed model
  # ----------------------------
  fixed_part <- if (include_covariates) {
    "y ~ season * shelf_loc_200 + scale(depth_mid) + scale(distance_sample_km)"
  } else {
    "y ~ season * shelf_loc_200"
  }
  fml <- stats::as.formula(paste(fixed_part, "+", random_effect))
  
  m1 <- lmerTest::lmer(fml, data = df_mod, REML = TRUE)
  
  # ----------------------------
  # emmeans for season × shelf at typical covariates
  # ----------------------------
  at_list <- if (include_covariates) {
    list(`scale(depth_mid)` = 0, `scale(distance_sample_km)` = 0)
  } else {
    list()
  }
  
  emm <- emmeans::emmeans(m1, ~ season * shelf_loc_200, at = at_list)
  pred_df <- as.data.frame(emm)
  
  # handle CI column names from emmeans (large-n => asymp.LCL/UCL; otherwise lower.CL/upper.CL)
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      dplyr::rename(lower.CL = .data$asymp.LCL, upper.CL = .data$asymp.UCL)
  } else if (all(c("lower.CL", "upper.CL") %in% names(pred_df))) {
    # already fine; do nothing
  } else if (all(c("LCL", "UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      dplyr::rename(lower.CL = .data$LCL, upper.CL = .data$UCL)
  } else {
    stop(
      "emmeans output did not include recognizable CI columns. Found: ",
      paste(names(pred_df), collapse = ", ")
    )
  }
  
  # ----------------------------
  # Panel A: violins (raw distribution)
  # ----------------------------
  dodge_w <- 0.75
  
  p_violin <- ggplot2::ggplot(df, ggplot2::aes(x = .data$shelf_loc_200, y = .data$.y)) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$season),
      position = ggplot2::position_dodge(width = dodge_w),
      trim = FALSE,
      alpha = 0.65,
      linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(.data$shelf_loc_200, .data$season)),
      position = ggplot2::position_dodge(width = dodge_w),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(values = season_colors, name = "Season") +
    ggplot2::labs(
      title = paste0("Raw data (depth < ", depth_cut, " m) — ", string_ID),
      x = "Shelf location (200 m)",
      y = y_label_expr
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text  = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text  = ggplot2::element_text(size = 13)
    )
  
  # ----------------------------
  # Panel B: model predictions (compressed x + connectors)
  # ----------------------------
  pred_df2 <- pred_df %>%
    dplyr::mutate(
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf","on_shelf")),
      season        = factor(.data$season,        levels = c("Winter","Summer")),
      x_base_comp   = 1 + (as.numeric(.data$shelf_loc_200) - 1) * gap,
      x             = .data$x_base_comp + ifelse(.data$season == "Winter", -0.12, 0.12)
    )
  
  seg_df <- pred_df2 %>%
    dplyr::select(.data$shelf_loc_200, .data$season, .data$x, .data$emmean) %>%
    tidyr::pivot_wider(names_from = .data$season, values_from = c(.data$x, .data$emmean)) %>%
    dplyr::transmute(
      shelf_loc_200,
      x    = .data$x_Winter,
      xend = .data$x_Summer,
      y    = .data$emmean_Winter,
      yend = .data$emmean_Summer
    )
  
  p_model <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg_df,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      color = "black", linetype = "dashed", linewidth = 0.9
    ) +
    ggplot2::geom_errorbar(
      data = pred_df2,
      ggplot2::aes(x = .data$x, y = .data$emmean, ymin = .data$lower.CL, ymax = .data$upper.CL, color = .data$season),
      width = ci_cap_width, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = pred_df2,
      ggplot2::aes(x = .data$x, y = .data$emmean, color = .data$season),
      size = 3
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("off_shelf","on_shelf")
    ) +
    ggplot2::scale_color_manual(values = season_colors, name = "Season") +
    ggplot2::labs(
      title = "Mixed model (emmeans ± CI)",
      x = "Shelf location (200 m)",
      y = y_label_expr
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text  = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text  = ggplot2::element_text(size = 13)
    )
  
  # ----------------------------
  # Align y-limits across panels
  # ----------------------------
  if (is.null(y_limits)) {
    y_limits <- range(df$.y, pred_df2$lower.CL, pred_df2$upper.CL, na.rm = TRUE)
  }
  
  p_violin <- p_violin + ggplot2::coord_cartesian(ylim = y_limits)
  p_model  <- p_model  + ggplot2::coord_cartesian(ylim = y_limits)
  
  # Apply conditional y-scale if needed
  if (!is.null(y_scale_obj)) {
    p_violin <- p_violin + y_scale_obj
    p_model  <- p_model  + y_scale_obj
  }
  
  # ----------------------------
  # Combine (side-by-side)
  # ----------------------------
  p_combined <- p_violin + p_model + patchwork::plot_layout(widths = c(1.3, 1))
  
  print(p_combined)
  
  # ----------------------------
  # Save
  # ----------------------------
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_", PSS_var,
             if (log10_transform) "_log10" else "",
             "_violin_plusMixedModel_depthlt", depth_cut, "_", string_ID, ".png")
    )
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_path, p_combined, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved plot to: ", out_path)
  }
  
  # ----------------------------
  # Save model outputs (RDS + text)
  # ----------------------------
  if (save_model) {
    
    if (is.null(model_prefix)) {
      model_prefix <- paste0(
        "SPECTRA_", PSS_var,
        if (log10_transform) "_log10" else "",
        "_mixedModel_depthlt", depth_cut, "_", string_ID
      )
    }
    
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1) Full objects for exact reproducibility
    rds_path <- file.path(model_dir, paste0(model_prefix, ".rds"))
    save_obj <- list(
      model        = m1,
      formula      = fml,
      df_mod       = df_mod,
      emmeans_raw  = pred_df,   # as.data.frame(emm) before x compression
      emmeans_plot = pred_df2,  # with x positions
      seg_df       = seg_df,
      settings     = list(
        PSS_var = PSS_var,
        depth_cut = depth_cut,
        log10_transform = log10_transform,
        random_effect = random_effect,
        include_covariates = include_covariates,
        at_list = at_list
      )
    )
    saveRDS(save_obj, rds_path)
    
    # 2) Human-readable summaries
    sum_path   <- file.path(model_dir, paste0(model_prefix, "_summary.txt"))
    anova_path <- file.path(model_dir, paste0(model_prefix, "_anova.txt"))
    emm_path   <- file.path(model_dir, paste0(model_prefix, "_emmeans.csv"))
    
    # write summary + anova as text files
    capture.output(summary(m1), file = sum_path)
    capture.output(anova(m1),   file = anova_path)
    
    # emmeans table as csv
    utils::write.csv(pred_df, emm_path, row.names = FALSE)
    
    message("Saved model bundle to: ", rds_path)
    message("Saved summary to:      ", sum_path)
    message("Saved anova to:        ", anova_path)
    message("Saved emmeans to:      ", emm_path)
  }
  
  invisible(list(
    plot = p_combined,
    plot_violin = p_violin,
    plot_model = p_model,
    model = m1,
    emmeans = pred_df2,
    df_mod = df_mod
  ))
}
# end of Mixed model function

res_slope <- PSS_violin_plus_mixedModel(
  dat       = dat_ISIIS_plankton,
  string_ID = "ISIIS_plankton",
  PSS_var   = "slope",
  depth_cut = 50,
  output_dir = output_dir_figs_NB,
  y_limits       = c(-0.5, -1.25), 
  save_plot  = TRUE,
  save_model = TRUE
)

res_biovol <- PSS_violin_plus_mixedModel(
  dat       = dat_ISIIS_plankton,
  string_ID = "ISIIS_plankton",
  PSS_var   = "integrated_biovol",
  depth_cut = 50,
  log10_transform = TRUE,
  y_limits       = c(10, 14.5),  
  output_dir = output_dir_figs_NB,
  save_plot  = TRUE,
  save_model = TRUE
)

res_biovol_noCruise <- PSS_violin_plus_mixedModel(
  dat       = dat_ISIIS_plankton,
  string_ID = "ISIIS_plankton",
  PSS_var   = "integrated_biovol",
  random_effect = "(1 | Transect_ID)",
  depth_cut = 50,
  log10_transform = TRUE,
  y_limits       = c(8, 17),  
  output_dir = output_dir_figs_NB,
  save_plot  = TRUE,
  save_model = TRUE
)


res_intercept <- PSS_violin_plus_mixedModel(
  dat       = dat_ISIIS_plankton,
  string_ID = "ISIIS_plankton",
  PSS_var   = "intercept",
  depth_cut = 50,
  y_limits       = c(5,12.5),
  output_dir = output_dir_figs_NB,
  save_plot  = TRUE,
  save_model = TRUE
)



PSS_twoPanel_medianFits_bySeasonShelf <- function(dat,
                                                  string_ID = "dataset",
                                                  x_var = "logBiovol",
                                                  slope_var = "slope",
                                                  intercept_var = "intercept",
                                                  shelf_col = "shelf_loc_200",
                                                  season_col = "season",
                                                  season_levels = c("Winter", "Summer"),
                                                  shelf_levels = c("off_shelf", "on_shelf"),
                                                  depth_cut = 50,
                                                  season_colors = c("Winter" = "#0072B2",
                                                                    "Summer" = "#E69F00"),
                                                  y_pad_frac = 0.05,   # 5% padding on y-range
                                                  output_dir = NULL,
                                                  save_plot = FALSE,
                                                  out_file = NULL,
                                                  width = 12,
                                                  height = 6,
                                                  dpi = 300) {
  
  # ----------------------------
  # Validate required columns
  # ----------------------------
  req <- c(x_var, slope_var, intercept_var, shelf_col, season_col)
  missing_cols <- setdiff(req, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in `dat`: ", paste(missing_cols, collapse = ", "))
  }
  
  # ----------------------------
  # Clean + factorize
  # ----------------------------
  df0 <- dat %>%
    dplyr::filter(depth_mid < depth_cut) %>%
    dplyr::mutate(
      .x      = suppressWarnings(as.numeric(.data[[x_var]])),
      .slope  = suppressWarnings(as.numeric(.data[[slope_var]])),
      .int    = suppressWarnings(as.numeric(.data[[intercept_var]])),
      .shelf  = factor(.data[[shelf_col]], levels = shelf_levels),
      .season = factor(.data[[season_col]], levels = season_levels)
    ) %>%
    dplyr::filter(
      is.finite(.x), is.finite(.slope), is.finite(.int),
      !is.na(.shelf), !is.na(.season)
    )
  
  if (nrow(df0) == 0) {
    stop("After filtering to finite x/slope/intercept and non-NA shelf/season, no rows remain.")
  }
  
  # ----------------------------
  # x-range from THIS dataset
  # ----------------------------
  xr <- range(df0$.x, na.rm = TRUE)
  
  # ----------------------------
  # MEDIAN slope/intercept per shelf × season
  # ----------------------------
  df_stats <- df0 %>%
    dplyr::group_by(.shelf, .season) %>%
    dplyr::summarise(
      slope_med = stats::median(.slope, na.rm = TRUE),
      int_med   = stats::median(.int,   na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    )
  
  # Expand to endpoints to draw lines across xr
  df_lines <- df_stats %>%
    tidyr::crossing(.x = xr) %>%
    dplyr::mutate(logNB_fit = int_med + slope_med * .x)
  
  # ----------------------------
  # Tight y-limits from fitted lines only
  # ----------------------------
  y_rng <- range(df_lines$logNB_fit, na.rm = TRUE)
  y_span <- diff(y_rng)
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  y_pad <- y_pad_frac * y_span
  y_limits <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  
  # ----------------------------
  # Plot
  # ----------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df_lines,
      ggplot2::aes(x = .x, y = logNB_fit, color = .season, group = .season),
      linewidth = 1.5, linetype = "dashed", alpha = 0.95
    ) +
    ggplot2::facet_grid(. ~ .shelf, scales = "fixed") +
    ggplot2::scale_color_manual(values = season_colors, name = "Season") +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(
      title = paste0("Median NBSS fits by season × shelf — ", string_ID),
      x = expression(log[10]*"(Biovolume mid, " * mu * "m"^3 * ")"),
      y = expression(log[10]*"(Normalized biovolume)")
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = 14),
      axis.text.y  = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      
      axis.ticks.length = grid::unit(6, "pt"),
      axis.ticks = ggplot2::element_line(linewidth = 0.9),
      
      legend.title = ggplot2::element_text(size = 16, face = "bold"),
      legend.text  = ggplot2::element_text(size = 14),
      legend.key.size = grid::unit(1.1, "lines"),
      legend.position = "top",
      legend.direction = "horizontal",
      
      strip.text = ggplot2::element_text(size = 14, face = "bold")
    )
  
  print(p)
  
  # ----------------------------
  # Optionally save
  # ----------------------------
  if (save_plot) {
    if (is.null(out_file)) {
      if (is.null(output_dir)) stop("If save_plot=TRUE, provide output_dir or out_file.")
      out_file <- file.path(output_dir, paste0("NBSS_medianFits_season_shelf_", string_ID, ".png"))
    }
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved plot to: ", out_file)
  }
  
  return(list(
    plot = p,
    stats = df_stats,
    y_limits = y_limits
  ))
}

# Example
res <- PSS_twoPanel_medianFits_bySeasonShelf(
  dat = dat_ISIIS_plankton,
  string_ID = "ISIIS_plankton",
  save_plot = TRUE,
  output_dir = output_dir_figs_NB
)

res$stats
res$y_limits





#1- mean slope per transect, but distinguishing inshore vs offshore, one panel per cruise
# Ensure factors are in the right order
# Violin plot


#slope summary function
slope_summary_plot <- function(dat_ISIIS, string_ID){
    p_slope<- ggplot(dat_ISIIS, aes(x = slope, y = Transect_ID, fill = shelf_loc_200)) +
      geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.3) +
      facet_wrap(~ Cruise, ncol = 4) +
      coord_cartesian(xlim = c(-2, 0)) +
      scale_fill_brewer(palette = "Set2") +
  
      # GH → RR order (top → bottom)
      scale_y_discrete(
        limits = rev(c("GH","OB","CR","CM","NH","HH","RR")),
        expand = expansion(mult = c(0.01, 0.01))
       ) +
  
      labs(
        title = paste0('NBSS slope for ', string_ID),
        x = "Size-spectrum slope",
        y = "",
        fill = "Shelf location (200 m)"
        ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "top",                   # ⬅⬅ legend on top
        legend.direction = "horizontal",           # optional: lay out horizontally
        strip.background = element_rect(fill = "gray90"),
        axis.text.y = element_text(size = 18),
        strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(1.2, "lines")
      )

    print(p_slope)

    out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_slope_summary_", string_ID, ".png"))
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

    ggsave(filename = out_path, plot = p_slope, width = 9, height = 15, units = "in", dpi = 300, bg = "white")
}


PSS_panel_plot <- function(dat_ISIIS, string_ID, PSS_var = "slope") {
  
  library(tidyverse)
  
  # ------------------------------------------
  # Validate PSS_var input
  # ------------------------------------------
  if (!PSS_var %in% c("slope", "intercept")) {
    stop("PSS_var must be 'slope' or 'intercept'")
  }
  
  # Name of variable being plotted
  yvar <- sym(PSS_var)
  
  # Nice label
  ylab <- if (PSS_var == "slope") {
    "Mean slope (± SD)"
  } else {
    "Mean intercept (± SD)"
  }
  
  # ------------------------------------------
  # Ordered levels
  # ------------------------------------------
  transect_levels <- c("GH","OB","CR","CM","NH","HH","RR")
  cruise_levels   <- c("w22", "s22", "w23", "s23")
  cruise_labels   <- c(
    "w22" = "Winter 2022",
    "s22" = "Summer 2022",
    "w23" = "Winter 2023",
    "s23" = "Summer 2023"
  )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(
      Transect_ID   = factor(Transect_ID, levels = transect_levels),
      Cruise        = factor(Cruise,      levels = cruise_levels),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("on_shelf","off_shelf"))
    )
  
  # ------------------------------------------
  # Summary stats for the chosen variable
  # ------------------------------------------
  df_sum <- dat_ISIIS %>%
    group_by(Transect_ID, Cruise, shelf_loc_200) %>%
    summarise(
      mean_val = mean(!!yvar, na.rm = TRUE),
      sd_val   = sd(!!yvar, na.rm = TRUE),
      n        = n(),
      se_val   = sd_val / sqrt(n),
      .groups  = "drop"
    )
  
  # ------------------------------------------
  # Plot
  # ------------------------------------------
  p <- ggplot(df_sum, aes(x = Cruise, y = mean_val, color = shelf_loc_200)) +
    
    # connecting lines across seasons
    geom_line(
      aes(group = shelf_loc_200, linetype = shelf_loc_200),
      position = position_dodge(width = 0.4),
      linewidth = 0.9
    ) +
    
    scale_linetype_manual(values = c("on_shelf" = "dotted",
                                     "off_shelf" = "dashed")) +
    
    # error bars (± SD)
    geom_errorbar(
      aes(ymin = mean_val - sd_val,
          ymax = mean_val + sd_val),
      width = 0.15,
      position = position_dodge(width = 0.4),
      linewidth = 0.5
    ) +
    
    # mean points
    geom_point(
      position = position_dodge(width = 0.4),
      size = 2.8,
      alpha = 0.9
    ) +
    
    facet_wrap(~ Transect_ID, ncol = 1, strip.position = "left") +
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(labels = cruise_labels) +
    
    labs(
      title = paste0("NBSS ", PSS_var, " across seasons & shelf location: ", string_ID),
      x = "Season",
      y = ylab,
      color = "Shelf location (200 m)",
      linetype = "Shelf location (200 m)"
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      strip.background = element_rect(fill = "gray95"),
      strip.placement = "outside",
      strip.text.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.2, "lines")
    )
  
  print(p)
  
  # ------------------------------------------
  # Save figure with correct variable in name
  # ------------------------------------------
  outfile <- file.path(
    output_dir_figs_NB,
    paste0("SPECTRA_", PSS_var, "_panels_", string_ID, ".png")
  )
  
  ggsave(outfile, p, width = 10, height = 22, dpi = 300, bg = "white")
  
  message("Saved plot to: ", outfile)
  return(p)
}



PSS_panel_plot(dat_ISIIS_plankton, 'plankton', 'slope')
PSS_panel_plot(dat_ISIIS_detritus, 'detritus', 'slope')
PSS_panel_plot(dat_ISIIS_full, 'all', 'slope')

PSS_panel_plot(dat_ISIIS_plankton, 'plankton', 'intercept')
PSS_panel_plot(dat_ISIIS_detritus, 'detritus', 'intercept')
PSS_panel_plot(dat_ISIIS_full, 'all', 'intercept')

PSS_panel_plot_depth <- function(dat_ISIIS, string_ID, PSS_var = "slope") {
  
  library(tidyverse)
  
  # ------------------------------------------
  # Validate PSS_var input
  # ------------------------------------------
  if (!PSS_var %in% c("slope", "intercept")) {
    stop("PSS_var must be 'slope' or 'intercept'")
  }
  
  yvar <- sym(PSS_var)
  
  # ------------------------------------------
  # Labels and fixed axis limits
  # ------------------------------------------
  if (PSS_var == "slope") {
    ylab <- "Mean slope (± SD)"
    y_limits <- c(-2.1, -0.5)
  } else {
    ylab <- "log10 intercept"
    y_limits <- c(5, 20)
  }
  
  # ------------------------------------------
  # Ordered levels
  # ------------------------------------------
  transect_levels <- c("GH","OB","CR","CM","NH","HH","RR")
  cruise_levels   <- c("w22", "s22", "w23", "s23")
  
  cruise_labels <- c(
    "w22" = "Winter 2022",
    "s22" = "Summer 2022",
    "w23" = "Winter 2023",
    "s23" = "Summer 2023"
  )
  
  # ------------------------------------------
  # Add Depth Interval Groups
  # ------------------------------------------
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(
      Transect_ID   = factor(Transect_ID, levels = transect_levels),
      Cruise        = factor(Cruise,      levels = cruise_levels),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf","on_shelf")),
      
      depth_group = case_when(
        depth_mid < 50  ~ "3–50 m",
        depth_mid >= 50 ~ "50–100 m",
        TRUE ~ NA_character_
      ),
      depth_group = factor(depth_group, levels = c("3–50 m", "50–100 m"))
    ) %>%
    drop_na(depth_group)
  
  # ------------------------------------------
  # Summary statistics
  # ------------------------------------------
  df_sum <- dat_ISIIS %>%
    group_by(Transect_ID, shelf_loc_200, Cruise, depth_group) %>%
    summarise(
      mean_val = mean(!!yvar, na.rm = TRUE),
      sd_val   = sd(!!yvar, na.rm = TRUE),
      n        = n(),
      se_val   = sd_val / sqrt(n),
      .groups  = "drop"
    )
  
  # ------------------------------------------
  # Build plot
  # ------------------------------------------
  dodge_amt <- 0.4
  
  p <- ggplot(df_sum, aes(x = Cruise, y = mean_val, color = depth_group)) +
    
    # connecting lines
    geom_line(
      aes(group = depth_group, linetype = depth_group),
      position = position_dodge(width = dodge_amt),
      linewidth = 1.2
    ) +
    
    scale_linetype_manual(values = c(
      "3–50 m"   = "dotted",
      "50–100 m" = "dashed"
    )) +
    
    # whiskers
    geom_errorbar(
      aes(ymin = mean_val - sd_val,
          ymax = mean_val + sd_val),
      width = 0.15,
      position = position_dodge(width = dodge_amt),
      linewidth = 0.8
    ) +
    
    # points
    geom_point(
      position = position_dodge(width = dodge_amt),
      size = 3.5,
      alpha = 0.9
    ) +
    
    facet_grid(
      rows = vars(Transect_ID),
      cols = vars(shelf_loc_200),
      switch = "y"
    ) +
    
    scale_color_brewer(palette = "Set2") +
    scale_x_discrete(labels = cruise_labels) +
    
    # FIXED AXIS LIMITS
    coord_cartesian(ylim = y_limits) +
    
    labs(
      title = paste0("NBSS ", PSS_var," across seasons, shelf location, and depth: ", string_ID),
      x = "Season",
      y = ylab,
      color = "Depth interval (m)",
      linetype = "Depth interval (m)"
    ) +
    
    # ------------------------------------------
  # Global larger fonts + cleaner aesthetic
  # ------------------------------------------
  theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      legend.box      = "horizontal",
      
      strip.background = element_rect(fill = "gray90"),
      strip.text       = element_text(size = 18, face = "bold"),
      
      axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 16),
      
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      
      plot.title  = element_text(size = 22, face = "bold"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.4, "lines")
    )
  
  print(p)
  
  # ------------------------------------------
  # Save
  # ------------------------------------------
  outfile <- file.path(
    output_dir_figs_NB,
    paste0("SPECTRA_", PSS_var, "_DepthSeasonShelf_panels_", string_ID, ".png")
  )
  
  ggsave(outfile, p, width = 14, height = 22, dpi = 300, bg = "white")
  message("Saved plot to: ", outfile)
  
  return(p)
}


PSS_panel_plot_depth(dat_ISIIS_plankton, 'plankton', 'slope')
PSS_panel_plot_depth(dat_ISIIS_detritus, 'detritus', 'slope')
PSS_panel_plot_depth(dat_ISIIS_full, 'all', 'slope')

PSS_panel_plot_depth(dat_ISIIS_plankton, 'plankton', 'intercept')
PSS_panel_plot_depth(dat_ISIIS_detritus, 'detritus', 'intercept')
PSS_panel_plot_depth(dat_ISIIS_full, 'all', 'intercept')




# plotting variance
variance_panel_plot <- function(dat_ISIIS, string_ID, PSS_var = "slope") {
  
  library(tidyverse)
  
  # ------------------------------------------
  # Validate PSS_var
  # ------------------------------------------
  if (!PSS_var %in% c("slope", "intercept")) {
    stop("PSS_var must be 'slope' or 'intercept'")
  }
  
  yvar <- sym(PSS_var)
  
  ylab <- if (PSS_var == "slope") {
    "Variance of NBSS slope"
  } else {
    "Variance of NBSS intercept"
  }
  
  # ------------------------------------------
  # Ordered factor levels
  # ------------------------------------------
  transect_levels <- c("GH","OB","CR","CM","NH","HH","RR")
  cruise_levels   <- c("w22", "s22", "w23", "s23")
  cruise_labels   <- c(
    "w22" = "Winter 2022",
    "s22" = "Summer 2022",
    "w23" = "Winter 2023",
    "s23" = "Summer 2023"
  )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(
      Transect_ID   = factor(Transect_ID, levels = transect_levels),
      Cruise        = factor(Cruise,      levels = cruise_levels),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf","off_shelf"))
    )
  
  # ------------------------------------------
  # Compute variance per group
  # ------------------------------------------
  df_var <- dat_ISIIS %>%
    group_by(Transect_ID, Cruise, shelf_loc_200) %>%
    summarise(
      var_val = var(!!yvar, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ------------------------------------------
  # Plot variance
  # ------------------------------------------
  p <- ggplot(df_var, aes(x = Cruise, y = var_val, color = shelf_loc_200)) +
    
    # connecting lines
    geom_line(
      aes(group = shelf_loc_200, linetype = shelf_loc_200),
      position = position_dodge(width = 0.4),
      linewidth = 0.9
    ) +
    
    scale_linetype_manual(values = c("on_shelf" = "dotted",
                                     "off_shelf" = "dashed")) +
    
    # points
    geom_point(
      position = position_dodge(width = 0.4),
      size = 2.8,
      alpha = 0.9
    ) +
    
    facet_wrap(~ Transect_ID, ncol = 1, strip.position = "left") +
    
    scale_color_brewer(palette = "Set2") +
    
    scale_x_discrete(labels = cruise_labels) +
    
    labs(
      title = paste0("NBSS ", PSS_var, " variance across seasons & shelf location: ", string_ID),
      x = "Season",
      y = ylab,
      color = "Shelf location (200 m)",
      linetype = "Shelf location (200 m)"
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      legend.position  = "top",
      legend.box       = "horizontal",
      strip.background = element_rect(fill = "gray95"),
      strip.placement  = "outside",
      strip.text.y     = element_text(size = 14, face = "bold"),
      axis.text.x      = element_text(size = 11, angle = 0, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.2, "lines")
    )
  
  print(p)
  
  out_path <- file.path(output_dir_figs_NB, paste0("SPECTRA_", PSS_var, "_variance_panels_", string_ID, ".png"))
  ggsave(out_path, p, width = 10, height = 22, dpi = 300, bg = "white")
  
  message("Saved variance plot to: ", out_path)
  return(p)
}


variance_panel_plot(dat_ISIIS_plankton, "plankton", PSS_var = "slope")



#intercept summary
p_intercept <- ggplot(dat_ISIIS, aes(x = intercept, y = Transect_ID, fill = shelf_loc_200)) +
  geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.3) +
  facet_wrap(~ Cruise, ncol = 4) +
  coord_cartesian(xlim = c(0, 20)) +
  scale_fill_brewer(palette = "Set2") +
  
  # GH → RR order (top → bottom)
  scale_y_discrete(
    limits = rev(c("GH","OB","CR","CM","NH","HH","RR")),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    x = "Size-spectrum intercept",
    y = "",
    fill = "Shelf location (200 m)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",                   # ⬅⬅ legend on top
    legend.direction = "horizontal",           # optional: lay out horizontally
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(size = 18),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1.2, "lines")
  )

p_intercept

out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_intercept_summary", ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = out_path, plot = p_intercept, width = 9, height = 15, units = "in", dpi = 300, bg = "white")

#GLOBAL SLOPE RELATIONSHIPS:
df_list_combined <- list.files(sizeSpectraDir, pattern = "ISIIS_LISST_paired\\.csv$", full.names = TRUE, recursive = TRUE)
pairs <- map_dfr( df_list_combined, ~ read_csv( .x, show_col_types = FALSE))

mod <- lm(slope_ISIIS ~ slope_LISST, data = pairs)

# Extract coefficients
a  <- coef(mod)[2]   # slope
b  <- coef(mod)[1]   # intercept
r2 <- summary(mod)$r.squared

# Create formatted equation label
eq_label <- sprintf("y = %.3f x + %.3f\nR² = %.3f", a, b, r2)

# Build the plot
p <- ggplot(pairs, aes(x = slope_LISST, y = slope_ISIIS, color = chl)) +
  geom_point(alpha = 0.85) +
  
  # Linear regression line
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Equation text added to bottom-right corner
  geom_text(
    aes(
      x = max(slope_LISST, na.rm = TRUE),
      y = min(slope_ISIIS, na.rm = TRUE),
      label = eq_label
    ),
    hjust = 1, vjust = 0,
    size = 5,
    inherit.aes = FALSE
  ) +
  
  scale_color_viridis_c(name = expression("Chl-a ("*mu*"g L"^{-1}*")")) +
  labs(x = "LISST slope", y = "ISIIS slope") +
  theme_classic() +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    aspect.ratio = 1
  ) +
  guides(color = guide_colorbar(
    barwidth  = unit(6, "cm"),
    barheight = unit(0.4, "cm")
  ))

p

out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_LISST_ISIIS_slopes.py", ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = out_path, plot = p, width = 9, height = 10, units = "in", dpi = 300, bg = "white")


####. FUNCTION and code to explore relationships between environmental variables and coefficients
plot_two_panel_regression <- function(df, predictor_col) {

  
  #---------------------------------------
  # Prepare long-format dataframe
  #---------------------------------------
  df_plot <- df %>%
    dplyr::select(all_of(predictor_col), slope, intercept, season, shelf_loc_200) %>%
    pivot_longer(
      cols = c(slope, intercept),
      names_to = "metric",
      values_to = "value"
    )
  
  # Nicely formatted panel titles
  df_plot$metric <- factor(df_plot$metric,
                           levels = c("slope", "intercept"),
                           labels = c("Slope", "Intercept"))
  
  # Regression equation helper
  lm_eqn <- function(x, y) {
    model <- lm(y ~ x)
    a <- coef(model)[2]
    b <- coef(model)[1]
    r2 <- summary(model)$r.squared
    sprintf("y = %.3f x + %.3f\nR² = %.3f", a, b, r2)
  }
  
  #---------------------------------------
  # Core plot with GLOBAL regression only
  #---------------------------------------
  p_two <- ggplot(
    df_plot,
    aes_string(
      x = predictor_col,     # dynamic aesthetic
      y = "value",
      color = "season",
      shape = "shelf_loc_200"
    )
  ) +
    
    # points
    geom_point(alpha = 0.3, size = 2) +
    
    # ONE global regression per panel (black dashed)
    geom_smooth(
      data = df_plot,
      aes_string(
        x = predictor_col,
        y = "value"
      ),
      method = "lm",
      se = FALSE,
      linewidth = 1,
      color = "black",
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~ metric, ncol = 2, scales = "free_y") +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = predictor_col,
      y = "Statistic",
      color = "Season",
      shape = "Shelf location (200 m)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold")
    )
  
  #---------------------------------------
  # Regression equation labels
  #---------------------------------------
  eq_labels <- df_plot %>%
    group_by(metric) %>%
    summarise(
      eq = lm_eqn(
        x = .data[[predictor_col]],   # values inside this metric group
        y = value
      ),
      x = max(.data[[predictor_col]], na.rm = TRUE),
      y = min(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_two_final <- p_two +
    geom_text(
      data = eq_labels,
      aes(x = x, y = y, label = eq),
      hjust = 1, vjust = 0, size = 3.5,
      inherit.aes = FALSE
    )
  
  #save fig
  out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_NBSS_coefficients_", predictor_col, ".png"))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = out_path, plot = p_two_final, width = 12, height = 9, units = "in", dpi = 300, bg = "white")
  
  print(p_two_final)
  return(p_two_final)
}

plot_regression_2x4 <- function(df, predictor_col,
                                output_dir = output_dir_figs_NB) {
  
  library(tidyverse)
  
  #-----------------------------------------
  # Sanity checks
  #-----------------------------------------
  if (!predictor_col %in% names(df)) {
    stop(paste("Column", predictor_col, "not found in dataframe"))
  }
  if (!all(c("season","shelf_loc_200","slope","intercept") %in% names(df))) {
    stop("df must contain: season, shelf_loc_200, slope, intercept")
  }
  
  #-----------------------------------------
  # Prep dataframe
  #-----------------------------------------
  df_plot <- df %>%
    dplyr::select(all_of(predictor_col), slope, intercept, season, shelf_loc_200) %>%
    pivot_longer(
      cols = c(slope, intercept),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(metric,
                      levels = c("slope", "intercept"),
                      labels = c("Slope", "Intercept")),
      season = factor(season),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("on_shelf","off_shelf"))
    )
  
  #-----------------------------------------
  # Regression helper: equation + R2 + p-value
  #-----------------------------------------
  lm_eqn <- function(x, y) {
    mod  <- lm(y ~ x)
    summ <- summary(mod)
    
    a    <- coef(mod)[2]            # slope
    b    <- coef(mod)[1]            # intercept
    r2   <- summ$r.squared          # R2
    pval <- summ$coefficients[2, 4] # p-value for slope term
    
    sprintf(
      "y = %.3f x + %.3f\nR² = %.3f\np = %.3f",
      a, b, r2, pval
    )
  }
  
  # Compute text locations and regression per facet
  stat_labels <- df_plot %>%
    group_by(metric, season, shelf_loc_200) %>%
    summarise(
      eq  = lm_eqn(.data[[predictor_col]], value),
      x   = max(.data[[predictor_col]], na.rm = TRUE),
      y   = min(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  #-----------------------------------------
  # Build plot
  #-----------------------------------------
  p_2x4 <- ggplot(df_plot,
                  aes_string(x = predictor_col, y = "value")) +
    
    geom_point(alpha = 0.35, size = 1.8,
               aes(color = season)) +
    
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.9,
      color = "black",
      linetype = "dashed"
    ) +
    
    facet_grid(
      rows = vars(metric),
      cols = vars(season, shelf_loc_200),
      scales = "free_y"
    ) +
    
    labs(
      x = predictor_col,
      y = "Value",
      title = paste("Slope and Intercept vs", predictor_col)
    ) +
    
    scale_color_brewer(palette = "Dark2") +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines")
    ) +
    
    geom_text(
      data = stat_labels,
      aes(label = eq),
      x     = Inf,
      y     = -Inf,
      hjust = 1.1,
      vjust = -0.4,
      size  = 3.5,
      inherit.aes = FALSE
    )
  
  #-----------------------------------------
  # Save figure
  #-----------------------------------------
  outfile <- file.path(
    output_dir,
    paste0("Regression_SeasonXShelf_", predictor_col, ".png")
  )
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = outfile,
    plot = p_2x4,
    width = 16, height = 8,
    dpi = 300, bg = "white"
  )
  
  message("Saved 2×4 regression plot to: ", outfile)
  print(p_2x4)
  return(p_2x4)
}


#filtering based on CTD comparison
dat_ISIIS <- dat_ISIIS %>%
  dplyr::filter(
    Salinity <= 40 | is.na(Salinity),
    between(Temperature, -60, 60) | is.na(Temperature),
    between(Oxygen_ml_l, 0, 10) | is.na(Oxygen_ml_l),
    pH <= 8.25 | is.na(pH)
    
    
  )


plot_regression_2x4(dat_ISIIS, "Temperature")
plot_regression_2x4(dat_ISIIS, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS, "pH")

#with subsetted data

dat_ISIIS_CR <- dat_ISIIS[dat_ISIIS$Transect_ID=='CR',]
dat_ISIIS_NH <- dat_ISIIS[dat_ISIIS$Transect_ID=='NH',]
dat_ISIIS_RR <- dat_ISIIS[dat_ISIIS$Transect_ID=='RR',]

plot_regression_2x4(dat_ISIIS_CR, "Temperature")
plot_regression_2x4(dat_ISIIS_CR, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS_CR, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS_CR, "pH")

plot_regression_2x4(dat_ISIIS_NH, "Temperature")
plot_regression_2x4(dat_ISIIS_NH, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS_NH, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS_NH, "pH")

plot_regression_2x4(dat_ISIIS_RR, "Temperature")
plot_regression_2x4(dat_ISIIS_RR, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS_RR, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS_RR, "pH")



plot_slope_vs_intercept_2x4 <- function(df, output_dir = output_dir_figs_NB) {
  
  library(tidyverse)
  
  req_cols <- c("season", "shelf_loc_200", "slope", "intercept")
  if (!all(req_cols %in% names(df))) {
    stop(paste("df must contain:", paste(req_cols, collapse=", ")))
  }
  
  df_plot <- df %>%
    dplyr::filter(
      !is.na(slope),
      !is.na(intercept),
      is.finite(slope),
      is.finite(intercept)
    ) %>%
    dplyr::select(season, shelf_loc_200, slope, intercept) %>%
    rename(
      intercept_val = intercept,
      slope_val     = slope
    ) %>%
    dplyr::mutate(
      season = factor(season),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("on_shelf","off_shelf"))
    )
  
  lm_eqn <- function(x, y) {
    mod <- lm(y ~ x)
    sprintf(
      "y = %.3f x + %.3f\nR² = %.3f",
      coef(mod)[2], coef(mod)[1], summary(mod)$r.squared
    )
  }
  
  df_labels <- df_plot %>%
    group_by(season, shelf_loc_200) %>%
    summarise(
      eq = lm_eqn(intercept_val, slope_val),
      x = Inf,
      y = -Inf,
      .groups = "drop"
    )
  
  p <- ggplot(df_plot,
              aes(x = intercept_val, y = slope_val)) +
    geom_point(alpha = 0.35, size = 1.8,
               aes(color = season)) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.9,
      color = "black",
      linetype = "dashed"
    ) +
    facet_grid(
      rows = vars(season),
      cols = vars(shelf_loc_200),
      scales = "free"
    ) +
    labs(
      title = "Slope as a function of Intercept",
      x = "Intercept",
      y = "Slope"
    ) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines")
    ) +
    geom_text(
      data = df_labels,
      aes(label = eq, x = x, y = y),
      hjust = 1.1,
      vjust = -0.4,
      size = 3.5,
      inherit.aes = FALSE
    )
  
  out <- file.path(output_dir, "Slope_vs_Intercept_bySeasonShelf.png")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(out, plot = p, width = 12, height = 8,
         bg = "white", dpi = 300)
  
  message("Saved slope-vs-intercept plot to: ", out)
  print(p)
  return(p)
}




plot_slope_vs_intercept_2x4(dat_ISIIS)


### NBSS graph 
compare_LISST_ISIIS_panels <- function(dat_ISIIS, string_ID,
                                       season_filter = NULL,
                                       output_dir = output_dir_figs_NB) {

  #---------------------------------------------------
  # CONSTANTS & FACTORS
  #---------------------------------------------------
  transect_levels <- c("GH","OB","CR","CM","NH","HH","RR")
  shelf_levels    <- c("off_shelf","on_shelf")
  depth_levels    <- c("3–50 m","50–100 m")
  
  depth_source_colors <- c(
    "LISST.3–50 m"    = "#009E73",
    "LISST.50–100 m"  = "#006400",
    "ISIIS.3–50 m"    = "#56B4E9",
    "ISIIS.50–100 m"  = "#CC79A7"
  )
  
  depth_linetypes <- c(
    "3–50 m"   = "dotted",
    "50–100 m" = "dashed"
  )
  
  #---------------------------------------------------
  # SEASON SUBSETTING
  #---------------------------------------------------
  if (!is.null(season_filter)) {
    if (season_filter == "Winter") {
      keep <- c("w22","w23")
    } else if (season_filter == "Summer") {
      keep <- c("s22","s23")
    } else {
      stop("season_filter must be NULL, 'Winter', or 'Summer'")
    }
    
    dat_ISIIS <- dat_ISIIS %>% dplyr::filter(Cruise %in% keep)
    dat_LISST <- dat_LISST %>% dplyr::filter(Cruise %in% keep)
    
    title_suffix <- paste(string_ID, "-", season_filter)
    file_suffix  <- paste0(string_ID, "_", season_filter)
  } else {
    title_suffix <- string_ID
    file_suffix  <- string_ID
  }
  
  #---------------------------------------------------
  # PREP ISIIS DATA
  #---------------------------------------------------
  isiis <- dat_ISIIS %>%
    dplyr::mutate(
      Transect_ID   = factor(Transect_ID, levels = transect_levels),
      shelf_loc_200 = factor(shelf_loc_200, levels = shelf_levels),
      depth_group = case_when(
        depth_mid < 50  ~ "3–50 m",
        depth_mid >= 50 ~ "50–100 m",
        TRUE ~ NA_character_
      ),
      depth_group = factor(depth_group, levels = depth_levels)
    ) %>%
    drop_na(depth_group) %>%
    transmute(
      source     = "ISIIS",
      Transect_ID,
      shelf_loc_200,
      depth_group,
      group_id,
      size_num = as.numeric(ECD_um_mid),
      NB       = as.numeric(NB),
      slope,
      intercept
    ) %>%
    dplyr::filter(size_num > 0, NB > 0)
  
  #---------------------------------------------------
  # PREP LISST DATA
  #---------------------------------------------------
  lisst <- dat_LISST %>%
    dplyr::mutate(
      Transect_ID   = factor(Transect_ID, levels = transect_levels),
      shelf_loc_200 = factor(shelf_loc_200, levels = shelf_levels),
      depth_group = case_when(
        depth_mid < 50  ~ "3–50 m",
        depth_mid >= 50 ~ "50–100 m",
        TRUE ~ NA_character_
      ),
      depth_group = factor(depth_group, levels = depth_levels)
    ) %>%
    drop_na(depth_group) %>%
    transmute(
      source     = "LISST",
      Transect_ID,
      shelf_loc_200,
      depth_group,
      group_id,
      size_num = as.numeric(size_class_mid_um),
      NB       = as.numeric(NB),
      slope,
      intercept
    ) %>%
    dplyr::filter(size_num > 0, NB > 0)
  
  #---------------------------------------------------
  # MERGE
  #---------------------------------------------------
  all_df <- bind_rows(lisst, isiis)
  
  levels_x <- sort(unique(all_df$size_num))
  all_df <- all_df %>% dplyr::mutate(size_fac = factor(size_num, levels = levels_x))
  
  #---------------------------------------------------
  # MEDIAN LINES (NOT mean)
  #---------------------------------------------------
  med_df <- all_df %>%
    dplyr::group_by(source, Transect_ID, shelf_loc_200, depth_group, size_fac) %>%
    dplyr::summarise(NB_med = median(NB, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(size_num = as.numeric(as.character(size_fac)))
  
  #---------------------------------------------------
  # POWER LAW LABEL
  #---------------------------------------------------
  compute_powerlaw_label <- function(df) {
    b <- median(df$intercept, na.rm = TRUE)
    s <- median(df$slope,     na.rm = TRUE)
    if (!is.finite(b) || !is.finite(s)) return("")
    paste0("y==", formatC(b, format="fg", digits=3),
           " * x^{", formatC(s, format="fg", digits=3), "}")
  }
  
  #---------------------------------------------------
  # PANEL RANGES FOR FIXED ANNOTATION POSITIONS
  #---------------------------------------------------
  panel_ranges <- all_df %>%
    group_by(Transect_ID, shelf_loc_200) %>%
    summarise(
      xmax = max(size_num, na.rm = TRUE),
      ymax = max(NB,       na.rm = TRUE),
      .groups = "drop"
    )
  
  offsets <- tibble(
    source = rep(c("LISST","ISIIS"), each = 2),
    depth_group = rep(depth_levels, times = 2),
    rank = 0:3
  )
  
  annot_df <- all_df %>%
    group_by(source, Transect_ID, shelf_loc_200, depth_group) %>%
    summarise(label = compute_powerlaw_label(cur_data()), .groups = "drop") %>%
    left_join(panel_ranges, by=c("Transect_ID","shelf_loc_200")) %>%
    left_join(offsets,       by=c("source","depth_group")) %>%
    dplyr::mutate(
      x_pos = xmax * 0.20,                 # moved left to avoid clipping
      y_pos = ymax * 1. * (10^(-rank))   # vertical stacking
    )
  
  #---------------------------------------------------
  # PLOT
  #---------------------------------------------------
  p <- ggplot() +
    
    geom_line(
      data = all_df,
      aes(x = size_num, y = NB,
          group = interaction(source, depth_group, group_id),
          color = interaction(source, depth_group),
          linetype = depth_group),
      linewidth = 0.25,
      alpha = 0.05
    ) +
    
    geom_line(
      data = med_df,
      aes(x = size_num, y = NB_med,
          color = interaction(source, depth_group),
          linetype = depth_group),
      linewidth = 1.3
    ) +
    
    geom_text(
      data = annot_df,
      aes(x = x_pos, y = y_pos, label = label,
          color = interaction(source, depth_group)),
      parse = TRUE,
      hjust = 0.1,       # small left shift
      vjust = 1,
      size = 6,
      fontface = "bold"
    ) +
    
    facet_grid(
      rows = vars(Transect_ID),
      cols = vars(shelf_loc_200),
      scales = "free_y",
      switch = "both"
    ) +
    
    scale_x_log10() +
    scale_y_log10(limits = c(1e-5, 1e11)) +
    
    scale_color_manual(values = depth_source_colors) +
    scale_linetype_manual(values = depth_linetypes) +
    
    labs(
      x = expression(Equivalent~spherical~diameter~(µm)),
      y = expression(Normalized~Biovolume~(µm^3~m^{-3}~µm^{-3})),
      title = paste("LISST–ISIIS Size Spectra Comparison:", title_suffix)
    ) +
    
    theme_classic(base_size = 20) +
    theme_classic(base_size = 20) +
    theme(
      axis.title.x = element_text(size = 22, margin = ggplot2::margin(t = 10)),
      axis.title.y = element_text(size = 22, margin = ggplot2::margin(r = 10)),
      
      # force x tick labels to show on EVERY facet row:
      axis.text.x  = element_text(size = 18, angle = 45, hjust = 1),
      axis.ticks.x = element_line(),
      
      axis.text.y  = element_text(size = 18),
      
      strip.text   = element_text(size = 20, face = "bold"),
      strip.placement = "outside",
      
      legend.text  = element_text(size = 18),
      legend.title = element_text(size = 20),
      legend.position = "top",
      
      panel.spacing = unit(1.5, "lines"),
      panel.spacing.y = unit(1.2, "lines")
    )
  
  
  print(p)
  
  #---------------------------------------------------
  # SAVE HIGH-RES PNG
  #---------------------------------------------------
  outfile <- file.path(output_dir,
                       paste0("LISST_ISIIS_panels_", file_suffix, ".png"))
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 18,
    height = 26,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  message("Saved plot to: ", outfile)
  return(p)
}




compare_LISST_ISIIS_panels(dat_ISIIS_full, "all", season_filter = "Winter")
compare_LISST_ISIIS_panels(dat_ISIIS_full, "all", season_filter = "Summer")

compare_LISST_ISIIS_panels(dat_ISIIS_detritus, "detritus", season_filter = "Winter")
compare_LISST_ISIIS_panels(dat_ISIIS_detritus, "detritus", season_filter = "Summer")

compare_LISST_ISIIS_panels(dat_ISIIS_plankton, "plankton", season_filter = "Winter")
compare_LISST_ISIIS_panels(dat_ISIIS_plankton, "plankton", season_filter = "Summer")



PSS_panel_plot_sources <- function(dat_full,
                                   dat_detritus,
                                   dat_plankton,
                                   shelf_loc = c("off_shelf", "on_shelf"),
                                   PSS_var = "slope",
                                   output_dir = output_dir_figs_NB) {
  
  library(tidyverse)
  
  # ----------------------------
  # Validate inputs
  # ----------------------------
  shelf_loc <- match.arg(shelf_loc)
  
  if (!PSS_var %in% c("slope", "intercept")) {
    stop("PSS_var must be 'slope' or 'intercept'")
  }
  
  yvar <- sym(PSS_var)
  
  # ------------------------------------------
  # Labels and fixed axis limits
  # ------------------------------------------
  if (PSS_var == "slope") {
    ylab <- "Mean slope (± SD)"
    y_limits <- c(-2.1, -0.5)
  } else {
    ylab <- "log10 intercept"
    y_limits <- c(5, 20)
  }
  
  # ------------------------------------------
  # Ordered levels
  # ------------------------------------------
  transect_levels <- c("GH","OB","CR","CM","NH","HH","RR")
  cruise_levels   <- c("w22", "s22", "w23", "s23")
  
  cruise_labels <- c(
    "w22" = "Winter 2022",
    "s22" = "Summer 2022",
    "w23" = "Winter 2023",
    "s23" = "Summer 2023"
  )
  
  # ----------------------------------------------------------
  # Helper to prepare each dataset (full, detritus, plankton)
  # ----------------------------------------------------------
  prep_fun <- function(df, source_label) {
    df %>%
      dplyr::mutate(
        source       = source_label,
        Transect_ID  = factor(Transect_ID, levels = transect_levels),
        Cruise       = factor(Cruise,      levels = cruise_levels),
        shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf","on_shelf")),
        
        depth_group = case_when(
          depth_mid < 50  ~ "3–50 m",
          depth_mid >= 50 ~ "50–100 m",
          TRUE ~ NA_character_
        ),
        
        depth_group = factor(depth_group, levels = c("3–50 m","50–100 m"))
      ) %>%
      drop_na(depth_group) %>%
      dplyr::filter(shelf_loc_200 == shelf_loc)
  }
  
  df_full     <- prep_fun(dat_full,     "All particles")
  df_det      <- prep_fun(dat_detritus, "Detritus")
  df_plankton <- prep_fun(dat_plankton, "Plankton")
  
  # Combine into a single long dataset
  dat_all <- bind_rows(df_full, df_det, df_plankton)
  
  # ----------------------------------------------------------
  # Summary statistics
  # ----------------------------------------------------------
  df_sum <- dat_all %>%
    group_by(Transect_ID, Cruise, depth_group, source) %>%
    summarise(
      mean_val = mean(!!yvar, na.rm = TRUE),
      sd_val   = sd(!!yvar, na.rm = TRUE),
      n        = n(),
      se_val   = sd_val / sqrt(n),
      .groups  = "drop"
    )
  
  # ----------------------------------------------------------
  # Aesthetics for data source
  # ----------------------------------------------------------
  source_colors <- c(
    "All particles" = "#1f78b4",  # blue
    "Detritus"      = "#b2df8a",  # greenish
    "Plankton"      = "#e31a1c"   # red
  )
  
  source_linetypes <- c(
    "All particles" = "solid",
    "Detritus"      = "dashed",
    "Plankton"      = "dotdash"
  )
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  dodge_amt <- 0.45
  
  p <- ggplot(df_sum, aes(x = Cruise, y = mean_val,
                          color = source, linetype = source)) +
    
    # connecting lines
    geom_line(
      aes(group = interaction(source, depth_group)),
      position = position_dodge(width = dodge_amt),
      linewidth = 1.2
    ) +
    
    # whiskers
    geom_errorbar(
      aes(ymin = mean_val - sd_val,
          ymax = mean_val + sd_val),
      width = 0.15,
      position = position_dodge(width = dodge_amt),
      linewidth = 0.8
    ) +
    
    # points
    geom_point(
      position = position_dodge(width = dodge_amt),
      size = 3.7,
      alpha = 0.9
    ) +
    
    facet_grid(
      rows = vars(Transect_ID),
      cols = vars(depth_group),
      switch = "y"
    ) +
    
    scale_color_manual(values = source_colors) +
    scale_linetype_manual(values = source_linetypes) +
    scale_x_discrete(labels = cruise_labels) +
    
    coord_cartesian(ylim = y_limits) +
    
    labs(
      title = paste0("NBSS ", PSS_var,
                     " across seasons and sources (", shelf_loc, ")"),
      x = "Season",
      y = ylab,
      color = "Data source",
      linetype = "Data source"
    ) +
    
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      legend.box      = "horizontal",
      
      strip.background = element_rect(fill = "gray90"),
      strip.text       = element_text(size = 18, face = "bold"),
      
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      
      plot.title  = element_text(size = 22, face = "bold"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.4, "lines")
    )
  
  print(p)
  
  # ----------------------------------------------------------
  # Save
  # ----------------------------------------------------------
  outfile <- file.path(
    output_dir,
    paste0("NBSS_", PSS_var, "_sources_", shelf_loc, ".png")
  )
  
  ggsave(outfile, p, width = 14, height = 22, dpi = 300, bg = "white")
  message("Saved plot to: ", outfile)
  
  return(p)
}

PSS_panel_plot_sources(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, shelf_loc = "on_shelf", PSS_var = "slope")
PSS_panel_plot_sources(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, shelf_loc = "off_shelf", PSS_var = "slope")

PSS_panel_plot_sources(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, shelf_loc = "on_shelf", PSS_var = "intercept")
PSS_panel_plot_sources(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, shelf_loc = "off_shelf", PSS_var = "intercept")


PSS_panel_plot_sources_test <- function(dat_full,
                                   dat_detritus,
                                   dat_plankton,
                                   dat_LISST,
                                   shelf_loc = c("off_shelf", "on_shelf"),
                                   PSS_var = "slope",
                                   output_dir = output_dir_figs_NB,
                                   save_table = TRUE) {
  

  # ----------------------------
  # Validate inputs
  # ----------------------------
  shelf_loc <- match.arg(shelf_loc)
  if (!PSS_var %in% c("slope", "intercept")) {
    stop("PSS_var must be 'slope' or 'intercept'")
  }
  yvar <- rlang::sym(PSS_var)
  
  # ------------------------------------------
  # Labels and fixed axis limits
  # ------------------------------------------
  if (PSS_var == "slope") {
    ylab <- "Mean slope (± SD)"
    y_limits <- c(-2.1, -0.5)
  } else {
    ylab <- "log10 intercept"
    y_limits <- c(5, 20)
  }
  
  # ------------------------------------------
  # Ordered levels
  # ------------------------------------------
  transect_levels <- c("GH","OB","CR","CM","NH","HH","RR")
  cruise_levels   <- c("w22", "s22", "w23", "s23")
  
  cruise_labels <- c(
    "w22" = "Winter 2022",
    "s22" = "Summer 2022",
    "w23" = "Winter 2023",
    "s23" = "Summer 2023"
  )
  
  # ----------------------------------------------------------
  # Helper to prepare each dataset (full, detritus, plankton, LISST)
  # ----------------------------------------------------------
  prep_fun <- function(df, source_label) {
    df %>%
      dplyr::mutate(
        source        = source_label,
        Transect_ID   = factor(Transect_ID, levels = transect_levels),
        Cruise        = factor(Cruise,      levels = cruise_levels),
        shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf","on_shelf")),
        depth_group = case_when(
          depth_mid < 50  ~ "3–50 m",
          depth_mid >= 50 ~ "50–100 m",
          TRUE ~ NA_character_
        ),
        depth_group = factor(depth_group, levels = c("3–50 m","50–100 m")),
        !!yvar := suppressWarnings(as.numeric(as.character(!!yvar)))  
      ) %>%
      drop_na(depth_group) %>%
      dplyr::filter(shelf_loc_200 == shelf_loc) %>%
      dplyr::filter(!is.na(!!yvar)) %>%
      dplyr::select(-dplyr::any_of("cast_id"))
  }
  
  
  df_full     <- prep_fun(dat_full,     "All particles")
  df_det      <- prep_fun(dat_detritus, "Detritus")
  df_plankton <- prep_fun(dat_plankton, "Plankton")
  df_lisst    <- prep_fun(dat_LISST,    "LISST")
  
  dat_all <- bind_rows(df_full, df_det, df_plankton, df_lisst)
  
  # ----------------------------------------------------------
  # Summary statistics (for plotting)
  # ----------------------------------------------------------
df_sum <- dat_all %>%
  dplyr::group_by(Transect_ID, Cruise, depth_group, source) %>%
  dplyr::summarise(
    mean_val = mean(.data[[.env$PSS_var]], na.rm = TRUE),
    sd_val   = sd(.data[[.env$PSS_var]], na.rm = TRUE),
    n        = dplyr::n(),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(se_val = sd_val / sqrt(n))
  
  # ----------------------------------------------------------
  # Panel-wise ANOVA + Tukey letters (compare SOURCE within each panel grouping)
  # Grouping: Transect_ID × depth_group × Cruise
  # ----------------------------------------------------------
  tukey_results <- dat_all %>%
    group_by(Transect_ID, depth_group, Cruise) %>%
    group_modify(function(.x, .g) {
      
      # Need at least 2 sources with data, and some replication
      if (dplyr::n_distinct(.x$source) < 2 || nrow(.x) < 4) {
        return(tibble(
          source = unique(.x$source),
          tukey_letter = NA_character_,
          anova_p = NA_real_
        ))
      }
      
      # ANOVA
      fit <- tryCatch(
        aov(reformulate("source", response = PSS_var), data = .x),
        error = function(e) NULL
        )
      
      if (is.null(fit)) {
        return(tibble(
          source = unique(.x$source),
          tukey_letter = NA_character_,
          anova_p = NA_real_
        ))
      }
      
      anova_p <- tryCatch(summary(fit)[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
      
      # Tukey letters (compact letter display)
      emm <- tryCatch(emmeans::emmeans(fit, ~ source), error = function(e) NULL)
      if (is.null(emm)) {
        return(tibble(
          source = unique(.x$source),
          tukey_letter = NA_character_,
          anova_p = anova_p
        ))
      }
      
      cld_tbl <- tryCatch(
        multcomp::cld(emm, Letters = letters, adjust = "tukey"),
        error = function(e) NULL
      )
      
      if (is.null(cld_tbl)) {
        return(tibble(
          source = unique(.x$source),
          tukey_letter = NA_character_,
          anova_p = anova_p
        ))
      }
      
      tibble(
        source = cld_tbl$source,
        tukey_letter = gsub("\\s+", "", cld_tbl$.group),
        anova_p = anova_p
      )
    }) %>%
    ungroup()
  
  # Join letters onto summary stats for plotting
  df_sum2 <- df_sum %>%
    left_join(tukey_results, by = c("Transect_ID","depth_group","Cruise","source"))
  
  # Position letters just above the errorbar/point
  df_sum2 <- df_sum2 %>%
    dplyr::mutate(label_y = pmin(mean_val + sd_val + 0.03 * diff(y_limits), y_limits[2]))
  
  # ----------------------------------------------------------
  # Aesthetics for data source
  # ----------------------------------------------------------
  source_colors <- c(
    "All particles" = "#1f78b4",
    "Detritus"      = "#b2df8a",
    "Plankton"      = "#e31a1c",
    "LISST"         = "#6a3d9a"
  )
  
  source_linetypes <- c(
    "All particles" = "solid",
    "Detritus"      = "dashed",
    "Plankton"      = "dotdash",
    "LISST"         = "twodash"
  )
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  dodge_amt <- 0.55
  
  p <- ggplot(df_sum2, aes(x = Cruise, y = mean_val,
                           color = source, linetype = source)) +
    
    geom_line(
      aes(group = source),
      position = position_dodge(width = dodge_amt),
      linewidth = 1.2
    ) +
    
    geom_errorbar(
      aes(ymin = mean_val - sd_val,
          ymax = mean_val + sd_val),
      width = 0.15,
      position = position_dodge(width = dodge_amt),
      linewidth = 0.8
    ) +
    
    geom_point(
      position = position_dodge(width = dodge_amt),
      size = 3.7,
      alpha = 0.9
    ) +
    
    # Tukey letters
    geom_text(
      aes(y = label_y, label = tukey_letter),
      position = position_dodge(width = dodge_amt),
      show.legend = FALSE,
      size = 4
    ) +
    
    facet_grid(
      rows = vars(Transect_ID),
      cols = vars(depth_group),
      switch = "y"
    ) +
    
    scale_color_manual(values = source_colors) +
    scale_linetype_manual(values = source_linetypes) +
    scale_x_discrete(labels = cruise_labels) +
    
    coord_cartesian(ylim = y_limits) +
    
    labs(
      title = paste0("NBSS ", PSS_var,
                     " across seasons and sources (", shelf_loc, ")"),
      x = "Season",
      y = ylab,
      color = "Data source",
      linetype = "Data source"
    ) +
    
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      legend.box      = "horizontal",
      
      strip.background = element_rect(fill = "gray90"),
      strip.text       = element_text(size = 18, face = "bold"),
      
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      
      plot.title  = element_text(size = 22, face = "bold"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.4, "lines")
    )
  
  print(p)
  
  # ----------------------------------------------------------
  # Save figure
  # ----------------------------------------------------------
  outfile <- file.path(
    output_dir,
    paste0("NBSS_", PSS_var, "_sources_plusLISST_", shelf_loc, ".png")
  )
  ggsave(outfile, p, width = 14, height = 22, dpi = 300, bg = "white")
  message("Saved plot to: ", outfile)
  
  # ----------------------------------------------------------
  # Build + save summary table (stats + ANOVA p + Tukey letters)
  # ----------------------------------------------------------
  stats_table <- df_sum2 %>%
    arrange(Transect_ID, depth_group, Cruise, source) %>%
    select(Transect_ID, depth_group, Cruise, source,
           n, mean_val, sd_val, se_val,
           anova_p, tukey_letter)
  
  if (save_table) {
    out_csv <- file.path(
      output_dir,
      paste0("NBSS_", PSS_var, "_sources_plusLISST_", shelf_loc, "_ANOVA_Tukey.csv")
    )
    readr::write_csv(stats_table, out_csv)
    message("Saved ANOVA/Tukey table to: ", out_csv)
  }
  
  return(list(plot = p, stats_table = stats_table))
}

# Example calls:

PSS_panel_plot_sources_test(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, dat_LISST, shelf_loc = "on_shelf", PSS_var = "slope")
PSS_panel_plot_sources_test(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, dat_LISST, shelf_loc = "off_shelf", PSS_var = "slope")

PSS_panel_plot_sources_test(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, dat_LISST, shelf_loc = "on_shelf", PSS_var = "intercept")
PSS_panel_plot_sources_test(dat_ISIIS_full, dat_ISIIS_detritus, dat_ISIIS_plankton, dat_LISST, shelf_loc = "off_shelf", PSS_var = "intercept")

# extra plot to see slope as function of chl

dat_ISIIS_plankton_2<- dat_ISIIS_plankton %>% filter(grouped_taxon=='phytoplankton'|
                                                       grouped_taxon=='diatoms'|
                                                       grouped_taxon=='protist')
mod <- lm(slope ~ chl_a_ug_l, data = dat_ISIIS_plankton_2)

# Extract coefficients
a  <- coef(mod)[2]   # slope
b  <- coef(mod)[1]   # intercept
r2 <- summary(mod)$r.squared

# Create formatted equation label
eq_label <- sprintf("y = %.3f x + %.3f\nR² = %.3f", a, b, r2)



p <- ggplot(dat_ISIIS_plankton_2, aes(x = chl_a_ug_l, y = slope)) +
  geom_point(alpha = 0.85) +
  
  # Linear regression line
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Equation text added to bottom-right corner
  geom_text(
    aes(
      x = max(chl_a_ug_l, na.rm = TRUE),
      y = min(slope, na.rm = TRUE),
      label = eq_label
    ),
    hjust = 1, vjust = 0,
    size = 15,
    inherit.aes = FALSE
  ) +
  

  labs(x = "chl a ugL^-1", y = "NBSS slope") +
  theme_classic() +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    aspect.ratio = 1
  ) +
  guides(color = guide_colorbar(
    barwidth  = unit(6, "cm"),
    barheight = unit(0.4, "cm")
  ))

p

out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_ISIIS_chl_all_phyto", ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(out_path, plot = p, width = 12, height = 8,
       bg = "white", dpi = 300)