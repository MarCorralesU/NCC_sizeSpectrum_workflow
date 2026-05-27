#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(lubridate)
  library(optparse)
  library(rlang)
  library(tidyr)
  library(grid)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

base_dir <- getwd()

option_list <- list(
  make_option(
    c("-e", "--envVars"),
    type = "character",
    default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l",
    help = "Comma-separated environmental variables",
    metavar = "character"
  ),
  make_option(
    c("--depthTol"),
    type = "double",
    default = 5,
    help = "Maximum allowed absolute depth difference (m) for pairing"
  ),
  make_option(
    c("--timeTolMin"),
    type = "double",
    default = 30,
    help = "Maximum allowed absolute time difference (minutes) for pairing"
  ),
  make_option(
    c("--latTol"),
    type = "double",
    default = 0.05,
    help = "Maximum allowed absolute latitude difference for pairing"
  ),
  make_option(
    c("--lonTol"),
    type = "double",
    default = 0.05,
    help = "Maximum allowed absolute longitude difference for pairing"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

env_vars <- strsplit(opt$envVars, ",")[[1]]

sizeSpectraDirPlankton <- file.path(base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded")
sizeSpectraDirAll      <- file.path(base_dir, "data", "ISIIS-LISST_mergedNBSS_5m_allParticles")

output_dir_figs <- file.path(base_dir, "figures", "size_spectra_paired_source_comparison")
output_dir_mods <- file.path(base_dir, "model_output", "size_spectra_paired_source_comparison")

dir.create(output_dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_mods, recursive = TRUE, showWarnings = FALSE)

df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_full     <- list.files(sizeSpectraDirAll,      pattern = "LISST-ISIIS_NBSS\\.csv$", full.names = TRUE, recursive = TRUE)

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
        dplyr::any_of(c("grouped_taxon", "season", "distance_sample_km")),
        dplyr::any_of(env_vars)
      )
    )
  )
  
  if (!"total_biovolume" %in% names(dat_ISIIS)) dat_ISIIS$total_biovolume <- NA_real_
  if (!"biovolume_um3"  %in% names(dat_ISIIS)) dat_ISIIS$biovolume_um3  <- NA_real_
  if (!"grouped_taxon"  %in% names(dat_ISIIS)) dat_ISIIS$grouped_taxon  <- NA_character_
  if (!"season"         %in% names(dat_ISIIS)) dat_ISIIS$season         <- NA_character_
  if (!"distance_sample_km" %in% names(dat_ISIIS)) dat_ISIIS$distance_sample_km <- NA_real_
  
  dat_ISIIS <- dat_ISIIS %>%
    mutate(
      total_biovolume = coalesce(total_biovolume, biovolume_um3)
    ) %>%
    select(-any_of("biovolume_um3")) %>%
    mutate(
      time_Pacific_orig = time_Pacific,
      mean_time_Pacific = time_Pacific - hours(1)
    ) %>%
    arrange(time_Pacific) %>%
    rename(group_id_old = group_id) %>%
    mutate(
      group_id = cumsum(coalesce(group_id_old != lag(group_id_old), TRUE)),
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
    group_by(group_id) %>%
    summarise(
      total_biovolume_sum = sum(total_biovolume, na.rm = TRUE),
      min_size = if (all(is.na(lower_size_edge))) NA_real_ else min(lower_size_edge, na.rm = TRUE),
      max_size = if (all(is.na(upper_size_edge))) NA_real_ else max(upper_size_edge, na.rm = TRUE),
      slope_g = first(slope),
      intercept_g = first(intercept),
      b_lin_g = first(b_lin),
      time_Pacific = first(time_Pacific),
      mean_time_Pacific = first(mean_time_Pacific),
      depth_mid = first(depth_mid),
      Latitude = first(Latitude),
      Longitude = first(Longitude),
      Cruise = first(Cruise),
      Transect_ID = first(Transect_ID),
      season = first(season),
      rel_dist = first(rel_dist),
      distance_sample_km = first(distance_sample_km),
      .groups = "drop"
    ) %>%
    mutate(
      integrated_biovol_all = if_else(
        is.finite(slope_g) & is.finite(b_lin_g),
        ifelse(
          abs(slope_g + 1) < 1e-12,
          b_lin_g * log(x2_full / x1_full),
          (b_lin_g / (slope_g + 1)) * (x2_full^(slope_g + 1) - x1_full^(slope_g + 1))
        ),
        NA_real_
      ),
      integrated_biovol_sample = if_else(
        is.finite(slope_g) & is.finite(b_lin_g) &
          is.finite(min_size) & is.finite(max_size) &
          min_size > 0 & max_size > 0 & max_size >= min_size,
        ifelse(
          abs(slope_g + 1) < 1e-12,
          b_lin_g * log(max_size / min_size),
          (b_lin_g / (slope_g + 1)) * (max_size^(slope_g + 1) - min_size^(slope_g + 1))
        ),
        NA_real_
      )
    ) %>%
    rename(
      slope = slope_g,
      intercept = intercept_g
    )
  
  group_summ
}

clean_ISIIS_dat <- function(dat_ISIIS) {
  dat_ISIIS <- dat_ISIIS %>%
    mutate(
      shelf_loc_200 = case_when(
        rel_dist >= 0 ~ "on_shelf",
        rel_dist <  0 ~ "off_shelf",
        TRUE          ~ NA_character_
      ),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      Cruise        = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      Transect_ID   = factor(Transect_ID, levels = c("GH", "OB", "CR", "CM", "NH", "HH", "RR")),
      season        = factor(season, levels = c("Winter", "Summer"))
    ) %>%
    filter(!(Cruise == "w22" & Transect_ID == "CM" & Longitude < -124.75))
  
  req_cols <- c("season", "distance_sample_km")
  miss <- req_cols[!req_cols %in% names(dat_ISIIS)]
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }
  
  if (all(is.na(dat_ISIIS$season))) stop("Column 'season' exists but is all NA.")
  if (all(is.na(dat_ISIIS$distance_sample_km))) stop("Column 'distance_sample_km' exists but is all NA.")
  
  dat_ISIIS
}

dat_all_raw      <- load_ISIIS_data(df_list_full, env_vars = env_vars)
dat_plankton_raw <- load_ISIIS_data(df_list_plankton, env_vars = env_vars)

dat_all      <- clean_ISIIS_dat(dat_all_raw)      %>% mutate(data_source = "allParticles")
dat_plankton <- clean_ISIIS_dat(dat_plankton_raw) %>% mutate(data_source = "plankton")

# ------------------------------------------------------------
# Pairing helper
# ------------------------------------------------------------
# Greedy nearest-neighbor matching:
# For each allParticles row, find the best unmatched plankton row
# within shared strata and tolerance windows.
# ------------------------------------------------------------

make_paired_matches <- function(dat_all,
                                dat_plankton,
                                depth_cut = 50,
                                response_var = "slope",
                                depth_tol = 5,
                                time_tol_min = 30,
                                lat_tol = 0.05,
                                lon_tol = 0.05) {
  
  needed <- c(
    response_var, "Cruise", "Transect_ID", "season", "shelf_loc_200",
    "depth_mid", "mean_time_Pacific", "Latitude", "Longitude",
    "distance_sample_km"
  )
  
  stopifnot(all(needed %in% names(dat_all)))
  stopifnot(all(needed %in% names(dat_plankton)))
  
  A <- dat_all %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      y = suppressWarnings(as.numeric(.data[[response_var]])),
      time_num = as.numeric(mean_time_Pacific)
    ) %>%
    filter(
      is.finite(y),
      !is.na(Cruise), !is.na(Transect_ID), !is.na(season), !is.na(shelf_loc_200),
      is.finite(depth_mid), is.finite(time_num),
      is.finite(Latitude), is.finite(Longitude)
    ) %>%
    arrange(Cruise, Transect_ID, season, shelf_loc_200, time_num, depth_mid) %>%
    mutate(row_id_all = row_number())
  
  P <- dat_plankton %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      y = suppressWarnings(as.numeric(.data[[response_var]])),
      time_num = as.numeric(mean_time_Pacific)
    ) %>%
    filter(
      is.finite(y),
      !is.na(Cruise), !is.na(Transect_ID), !is.na(season), !is.na(shelf_loc_200),
      is.finite(depth_mid), is.finite(time_num),
      is.finite(Latitude), is.finite(Longitude)
    ) %>%
    arrange(Cruise, Transect_ID, season, shelf_loc_200, time_num, depth_mid) %>%
    mutate(row_id_plankton = row_number())
  
  matched_plankton_ids <- integer(0)
  pair_list <- vector("list", nrow(A))
  pair_counter <- 0L
  
  for (i in seq_len(nrow(A))) {
    ai <- A[i, ]
    
    cand <- P %>%
      filter(
        !(row_id_plankton %in% matched_plankton_ids),
        Cruise == ai$Cruise,
        Transect_ID == ai$Transect_ID,
        season == ai$season,
        shelf_loc_200 == ai$shelf_loc_200
      ) %>%
      mutate(
        depth_diff = abs(depth_mid - ai$depth_mid),
        time_diff_min = abs(time_num - ai$time_num) / 60,
        lat_diff = abs(Latitude - ai$Latitude),
        lon_diff = abs(Longitude - ai$Longitude)
      ) %>%
      filter(
        depth_diff <= depth_tol,
        time_diff_min <= time_tol_min,
        lat_diff <= lat_tol,
        lon_diff <= lon_tol
      ) %>%
      mutate(
        # weighted match score; depth and time carry most of the weight
        match_score =
          (depth_diff / max(depth_tol, 1e-9)) +
          (time_diff_min / max(time_tol_min, 1e-9)) +
          0.5 * (lat_diff / max(lat_tol, 1e-9)) +
          0.5 * (lon_diff / max(lon_tol, 1e-9))
      ) %>%
      arrange(match_score, depth_diff, time_diff_min)
    
    if (nrow(cand) == 0) next
    
    best <- cand[1, ]
    matched_plankton_ids <- c(matched_plankton_ids, best$row_id_plankton)
    pair_counter <- pair_counter + 1L
    
    pair_list[[i]] <- tibble(
      pair_id = paste0("pair_", sprintf("%06d", pair_counter)),
      
      Cruise = as.character(ai$Cruise),
      Transect_ID = as.character(ai$Transect_ID),
      season = as.character(ai$season),
      shelf_loc_200 = as.character(ai$shelf_loc_200),
      
      depth_mid_all = ai$depth_mid,
      depth_mid_plankton = best$depth_mid,
      depth_diff = best$depth_diff,
      
      time_Pacific_all = ai$mean_time_Pacific,
      time_Pacific_plankton = best$mean_time_Pacific,
      time_diff_min = best$time_diff_min,
      
      Latitude_all = ai$Latitude,
      Latitude_plankton = best$Latitude,
      lat_diff = best$lat_diff,
      
      Longitude_all = ai$Longitude,
      Longitude_plankton = best$Longitude,
      lon_diff = best$lon_diff,
      
      distance_sample_km_all = ai$distance_sample_km,
      distance_sample_km_plankton = best$distance_sample_km,
      
      y_all = ai$y,
      y_plankton = best$y
    )
  }
  
  pairs_wide <- bind_rows(pair_list)
  
  if (nrow(pairs_wide) == 0) {
    stop("No matched pairs found under current tolerances.")
  }
  
  pairs_long <- pairs_wide %>%
    transmute(
      pair_id,
      Cruise,
      Transect_ID,
      season,
      shelf_loc_200,
      
      depth_mid_all,
      depth_mid_plankton,
      depth_diff,
      time_diff_min,
      lat_diff,
      lon_diff,
      
      data_source_all = "allParticles",
      data_source_plankton = "plankton",
      
      y_all,
      y_plankton,
      
      distance_sample_km_all,
      distance_sample_km_plankton
    ) %>%
    pivot_longer(
      cols = c(y_all, y_plankton,
               data_source_all, data_source_plankton,
               distance_sample_km_all, distance_sample_km_plankton),
      names_to = c(".value", "source_stub"),
      names_pattern = "(y|data_source|distance_sample_km)_(all|plankton)"
    ) %>%
    mutate(
      data_source = recode(source_stub,
                           all = "allParticles",
                           plankton = "plankton"),
      depth_mid = if_else(data_source == "allParticles", depth_mid_all, depth_mid_plankton),
      distance_sample_km = distance_sample_km,
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      Transect_ID = factor(Transect_ID, levels = c("GH", "OB", "CR", "CM", "NH", "HH", "RR")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      data_source = factor(data_source, levels = c("allParticles", "plankton"))
    ) %>%
    select(
      pair_id, Cruise, Transect_ID, season, shelf_loc_200,
      data_source, y, depth_mid, distance_sample_km,
      depth_diff, time_diff_min, lat_diff, lon_diff
    )
  
  list(
    pairs_wide = pairs_wide,
    pairs_long = pairs_long
  )
}

# ------------------------------------------------------------
# Paired mixed model
# ------------------------------------------------------------

run_paired_mm <- function(dat_all,
                          dat_plankton,
                          response_var = "slope",
                          depth_cut = 50,
                          log10_transform = FALSE,
                          depth_tol = 5,
                          time_tol_min = 30,
                          lat_tol = 0.05,
                          lon_tol = 0.05,
                          output_dir_figs,
                          output_dir_mods) {
  
  pair_obj <- make_paired_matches(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    depth_cut = depth_cut,
    response_var = response_var,
    depth_tol = depth_tol,
    time_tol_min = time_tol_min,
    lat_tol = lat_tol,
    lon_tol = lon_tol
  )
  
  df <- pair_obj$pairs_long
  
  if (log10_transform) {
    df <- df %>%
      mutate(y = if_else(y > 0, log10(y), NA_real_)) %>%
      filter(is.finite(y))
  }
  
  # fixed effects:
  # data_source + season + shelf + interactions + continuous covariates
  # random effects:
  # Cruise/Transect_ID and pair_id
  fml <- y ~ data_source * season * shelf_loc_200 +
    scale(depth_mid) + scale(distance_sample_km) +
    (1 | Cruise/Transect_ID) + (1 | pair_id)
  
  m1 <- lmerTest::lmer(fml, data = df, REML = TRUE)
  
  emm <- emmeans::emmeans(
    m1,
    ~ data_source * season * shelf_loc_200,
    at = list(
      `scale(depth_mid)` = 0,
      `scale(distance_sample_km)` = 0
    )
  )
  
  source_contrasts <- contrast(
    emm,
    method = "revpairwise",
    by = c("season", "shelf_loc_200")
  )
  
  emm_df <- as.data.frame(emm)
  contr_df <- as.data.frame(source_contrasts)
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(emm_df))) {
    emm_df <- emm_df %>%
      rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  } else if (all(c("LCL", "UCL") %in% names(emm_df))) {
    emm_df <- emm_df %>%
      rename(lower.CL = LCL, upper.CL = UCL)
  }
  
  # ------------------------
  # Plot
  # ------------------------
  season_colors <- c("Winter" = "#0072B2", "Summer" = "#E69F00")
  source_shapes <- c("allParticles" = 16, "plankton" = 17)
  
  if (response_var %in% c("integrated_biovol_all", "integrated_biovol_sample", "total_biovolume")) {
    y_lab <- if (log10_transform) expression(log[10](Biovolume~(mu*m^3~m^-3))) else expression(Biovolume~(mu*m^3~m^-3))
  } else if (response_var == "intercept") {
    y_lab <- expression(Intercept)
  } else if (response_var == "slope") {
    y_lab <- expression(Slope)
  } else {
    y_lab <- response_var
  }
  
  p1 <- ggplot(
    df,
    aes(
      x = shelf_loc_200,
      y = y,
      fill = season,
      color = data_source,
      group = interaction(shelf_loc_200, season, data_source)
    )
  ) +
    geom_violin(
      position = position_dodge(width = 0.75),
      trim = FALSE,
      alpha = 0.45,
      linewidth = 0.35
    ) +
    geom_boxplot(
      aes(group = interaction(shelf_loc_200, season, data_source)),
      position = position_dodge(width = 0.75),
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.30,
      linewidth = 0.30
    ) +
    scale_fill_manual(values = season_colors, name = "Season") +
    scale_color_manual(values = c("allParticles" = "black", "plankton" = "gray40"), name = "Source") +
    scale_x_discrete(labels = c("off_shelf" = "off shelf", "on_shelf" = "on shelf")) +
    labs(
      title = paste0(response_var, " — paired matched samples"),
      x = NULL,
      y = y_lab
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  p2 <- ggplot(
    emm_df,
    aes(
      x = shelf_loc_200,
      y = emmean,
      color = season,
      shape = data_source,
      group = interaction(season, data_source)
    )
  ) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      position = position_dodge(width = 0.60),
      width = 0.05,
      linewidth = 0.6
    ) +
    geom_point(
      position = position_dodge(width = 0.60),
      size = 3
    ) +
    scale_color_manual(values = season_colors, name = "Season") +
    scale_shape_manual(values = source_shapes, name = "Source") +
    scale_x_discrete(labels = c("off_shelf" = "off shelf", "on_shelf" = "on shelf")) +
    labs(
      title = "Estimated marginal means ± 95% CI",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  p1 <- p1 + annotate("text", x = -Inf, y = Inf, label = "a.", hjust = -0.6, vjust = 1.5, size = 6, fontface = "bold")
  p2 <- p2 + annotate("text", x = -Inf, y = Inf, label = "b.", hjust = -0.6, vjust = 1.5, size = 6, fontface = "bold")
  
  p_combined <- p1 + p2 + plot_layout(widths = c(1.35, 1))
  
  file_stub <- paste0(
    "paired_", response_var,
    if (log10_transform) "_log10" else "",
    "_depthlt", depth_cut,
    "_dTol", depth_tol,
    "_tTol", time_tol_min,
    "_latTol", lat_tol,
    "_lonTol", lon_tol
  )
  
  ggsave(
    filename = file.path(output_dir_figs, paste0(file_stub, ".png")),
    plot = p_combined,
    width = 12,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  write_csv(pair_obj$pairs_wide, file.path(output_dir_mods, paste0(file_stub, "_matched_pairs_wide.csv")))
  write_csv(df,                 file.path(output_dir_mods, paste0(file_stub, "_matched_pairs_long.csv")))
  write_csv(emm_df,             file.path(output_dir_mods, paste0(file_stub, "_emmeans.csv")))
  write_csv(contr_df,           file.path(output_dir_mods, paste0(file_stub, "_source_contrasts.csv")))
  
  capture.output(summary(m1), file = file.path(output_dir_mods, paste0(file_stub, "_summary.txt")))
  capture.output(anova(m1),   file = file.path(output_dir_mods, paste0(file_stub, "_anova.txt")))
  
  saveRDS(
    list(
      model = m1,
      pairs_wide = pair_obj$pairs_wide,
      pairs_long = df,
      emmeans = emm_df,
      contrasts = contr_df
    ),
    file.path(output_dir_mods, paste0(file_stub, ".rds"))
  )
  
  message("--------------------------------------------------")
  message("Response variable: ", response_var)
  message("Matched pairs: ", n_distinct(df$pair_id))
  message("Rows in long paired dataset: ", nrow(df))
  message("Saved outputs with stub: ", file_stub)
  message("--------------------------------------------------")
  
  invisible(list(
    model = m1,
    pairs_wide = pair_obj$pairs_wide,
    pairs_long = df,
    emmeans = emm_df,
    contrasts = contr_df,
    plot = p_combined
  ))
}

# ------------------------------------------------------------
# Run separate paired analyses
# ------------------------------------------------------------

res_slope <- run_paired_mm(
  dat_all = dat_all,
  dat_plankton = dat_plankton,
  response_var = "slope",
  depth_cut = 50,
  log10_transform = FALSE,
  depth_tol = opt$depthTol,
  time_tol_min = opt$timeTolMin,
  lat_tol = opt$latTol,
  lon_tol = opt$lonTol,
  output_dir_figs = output_dir_figs,
  output_dir_mods = output_dir_mods
)

res_intercept <- run_paired_mm(
  dat_all = dat_all,
  dat_plankton = dat_plankton,
  response_var = "intercept",
  depth_cut = 50,
  log10_transform = FALSE,
  depth_tol = opt$depthTol,
  time_tol_min = opt$timeTolMin,
  lat_tol = opt$latTol,
  lon_tol = opt$lonTol,
  output_dir_figs = output_dir_figs,
  output_dir_mods = output_dir_mods
)

res_intbio <- run_paired_mm(
  dat_all = dat_all,
  dat_plankton = dat_plankton,
  response_var = "integrated_biovol_all",
  depth_cut = 50,
  log10_transform = TRUE,
  depth_tol = opt$depthTol,
  time_tol_min = opt$timeTolMin,
  lat_tol = opt$latTol,
  lon_tol = opt$lonTol,
  output_dir_figs = output_dir_figs,
  output_dir_mods = output_dir_mods
)