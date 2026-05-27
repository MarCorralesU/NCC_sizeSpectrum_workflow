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
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt)

env_vars <- strsplit(opt$envVars, ",")[[1]] #env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]

sizeSpectraDirPlankton <- file.path(base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded")
sizeSpectraDirAll      <- file.path(base_dir, "data", "ISIIS-LISST_mergedNBSS_5m_allParticles")
sizeSpectraDirDetritus <- file.path(base_dir, "data", "size_spectra_detritus_5-10m-bins_thresholded")

output_dir_figs_NB <- file.path(base_dir, "figures", "size_spectra_mixed_models")
if (!dir.exists(output_dir_figs_NB)) {
  dir.create(output_dir_figs_NB, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_NB)
}

df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_full     <- list.files(sizeSpectraDirAll,      pattern = "LISST-ISIIS_NBSS\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_detritus <- list.files(sizeSpectraDirDetritus, pattern = "by_taxa\\.csv$", full.names = TRUE, recursive = TRUE)

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
    select(
      group_id,
      total_biovolume_sum,
      integrated_biovol_all,
      integrated_biovol_sample
    )
  
  dat_ISIIS <- dat_ISIIS %>%
    left_join(group_summ, by = "group_id") %>%
    mutate(
      slope_anomaly = slope + 1,
      NB = 10^(logNB),
      total_biovolume = total_biovolume_sum
    )
  
  return(dat_ISIIS)
}

dat_ISIIS_full_raw     <- load_ISIIS_data(df_list_full, env_vars = env_vars)
dat_ISIIS_plankton_raw <- load_ISIIS_data(df_list_plankton, env_vars = env_vars)
dat_ISIIS_detritus_raw <- load_ISIIS_data(df_list_detritus, env_vars = env_vars)




clean_ISIIS_dat <- function(dat_ISIIS) {
  dat_ISIIS <- dat_ISIIS %>%
    mutate(
      shelf_loc_200 = case_when(
        rel_dist >= 0 ~ "on_shelf",
        rel_dist <  0 ~ "off_shelf",
        TRUE          ~ NA_character_
      ),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf")),
      Cruise        = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      Transect_ID   = factor(Transect_ID, levels = c("GH", "OB", "CR", "CM", "NH", "HH", "RR")),
      season        = factor(season, levels = c("Winter", "Summer"))
    ) %>%
    filter(!(Cruise == "w22" & Transect_ID == "CM" & Longitude < -124.75))
  
  req_cols <- c("season", "distance_sample_km")
  miss <- req_cols[!req_cols %in% names(dat_ISIIS)]
  if (length(miss) > 0) {
    stop("Missing required columns for mixed models: ", paste(miss, collapse = ", "))
  }
  
  if (all(is.na(dat_ISIIS$season))) {
    stop("Column 'season' exists but is all NA.")
  }
  if (all(is.na(dat_ISIIS$distance_sample_km))) {
    stop("Column 'distance_sample_km' exists but is all NA.")
  }
  
  dat_ISIIS
}

dat_ISIIS_full     <- clean_ISIIS_dat(dat_ISIIS_full_raw)
dat_ISIIS_plankton <- clean_ISIIS_dat(dat_ISIIS_plankton_raw)
dat_ISIIS_detritus <- clean_ISIIS_dat(dat_ISIIS_detritus_raw)


PSS_violin_plus_mixedModel <- function(dat,
                                       string_ID,
                                       PSS_var = "slope",
                                       depth_cut = 100,
                                       log10_transform = FALSE,
                                       output_dir = output_dir_figs_NB,
                                       y_limits = NULL,
                                       random_effect = "(1 | Cruise/Transect_ID)",
                                       include_covariates = TRUE,
                                       season_colors = c("Winter" = "#0072B2",
                                                         "Summer" = "#E69F00"),
                                       gap = 0.35,
                                       ci_cap_width = 0.03,
                                       width = 12,
                                       height = 6,
                                       dpi = 300,
                                       save_plot = TRUE,
                                       save_model = TRUE,
                                       model_dir = output_dir_figs_NB,
                                       model_prefix = NULL) {
  
  needed <- c(PSS_var, "season", "shelf_loc_200", "depth_mid",
              "Cruise", "Transect_ID", "group_id", "distance_sample_km")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))
  
  yvar <- sym(PSS_var)
  
  df <- dat %>%
    filter(.data$depth_mid < depth_cut) %>%
    mutate(
      season = factor(.data$season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      Cruise = factor(.data$Cruise),
      Transect_ID = factor(.data$Transect_ID),
      .y_raw = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    filter(!is.na(.data$season), !is.na(.data$shelf_loc_200), is.finite(.data$.y_raw))
  
  if (nrow(df) == 0) stop("After filtering depth/NA/non-finite, no rows remain.")
  
  if (log10_transform) {
    df <- df %>%
      mutate(.y = if_else(.data$.y_raw > 0, log10(.data$.y_raw), NA_real_)) %>%
      filter(is.finite(.data$.y))
  } else {
    df <- df %>% mutate(.y = .data$.y_raw)
  }
  
  if (PSS_var %in% c("total_biovolume", "integrated_biovol_all", "integrated_biovol_sample")) {
    y_label_expr <- expression(Biovolume~(mu*m^3~m^-3))
  } else if (PSS_var == "intercept") {
    y_label_expr <- expression(Intercept~(mu*m^3~m^-3~mu*m^-3))
  } else if (PSS_var == "slope") {
    y_label_expr <- expression(Slope~(m^-3~mu*m^-3))
  } else {
    y_label_expr <- PSS_var
  }
  
  y_scale_obj <- if (log10_transform || PSS_var == "intercept") {
    ggplot2::scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", round(x, 2)))
    )
  } else {
    NULL
  }
  
  df_mod <- df %>%
    group_by(.data$Cruise, .data$Transect_ID, .data$group_id, .data$season, .data$shelf_loc_200) %>%
    summarise(
      y = first(.data$.y),
      depth_mid = first(.data$depth_mid),
      distance_sample_km = first(.data$distance_sample_km),
      .groups = "drop"
    ) %>%
    filter(is.finite(.data$y))
  
  fixed_part <- if (include_covariates) {
    "y ~ season * shelf_loc_200 + scale(depth_mid) + scale(distance_sample_km)"
  } else {
    "y ~ season * shelf_loc_200"
  }
  fml <- stats::as.formula(paste(fixed_part, "+", random_effect))
  
  m1 <- lmerTest::lmer(fml, data = df_mod, REML = TRUE)
  
  at_list <- if (include_covariates) {
    list(`scale(depth_mid)` = 0, `scale(distance_sample_km)` = 0)
  } else {
    list()
  }
  
  emm <- emmeans::emmeans(m1, ~ season * shelf_loc_200, at = at_list)
  pred_df <- as.data.frame(emm)
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$asymp.LCL, upper.CL = .data$asymp.UCL)
  } else if (all(c("LCL", "UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$LCL, upper.CL = .data$UCL)
  }
  
  dodge_w <- 0.75
  
  p_violin <- ggplot(df, aes(x = .data$shelf_loc_200, y = .data$.y)) +
    geom_violin(
      aes(fill = .data$season),
      position = position_dodge(width = dodge_w),
      trim = FALSE,
      alpha = 0.65,
      linewidth = 0.3
    ) +
    geom_boxplot(
      aes(group = interaction(.data$shelf_loc_200, .data$season)),
      position = position_dodge(width = dodge_w),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    scale_x_discrete(
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      )
    ) +
    scale_fill_manual(values = season_colors, name = "Season") +
    labs(
      title = paste0("Raw data (depth < ", depth_cut, " m) — ", string_ID),
      x = NULL,
      y = y_label_expr
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  pred_df2 <- pred_df %>%
    mutate(
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      season = factor(.data$season, levels = c("Winter", "Summer")),
      x_base_comp = 1 + (as.numeric(.data$shelf_loc_200) - 1) * gap,
      x = .data$x_base_comp + ifelse(.data$season == "Winter", -0.12, 0.12)
    )
  
  seg_df <- pred_df2 %>%
    select(.data$shelf_loc_200, .data$season, .data$x, .data$emmean) %>%
    pivot_wider(names_from = .data$season, values_from = c(.data$x, .data$emmean)) %>%
    transmute(
      shelf_loc_200,
      x = .data$x_Winter,
      xend = .data$x_Summer,
      y = .data$emmean_Winter,
      yend = .data$emmean_Summer
    )
  
  p_model <- ggplot() +
    geom_segment(
      data = seg_df,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      color = "black", linetype = "dashed", linewidth = 0.9
    ) +
    geom_errorbar(
      data = pred_df2,
      aes(x = .data$x, y = .data$emmean, ymin = .data$lower.CL, ymax = .data$upper.CL, color = .data$season),
      width = ci_cap_width, linewidth = 0.6
    ) +
    geom_point(
      data = pred_df2,
      aes(x = .data$x, y = .data$emmean, color = .data$season),
      size = 3
    ) +
    scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("off shelf", "on shelf")
    )+
    scale_color_manual(values = season_colors, name = "Season") +
    labs(
      title = "Mixed model (emmeans ± CI)",
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  if (is.null(y_limits)) {
    y_limits <- range(df$.y, pred_df2$lower.CL, pred_df2$upper.CL, na.rm = TRUE)
  }
  
  p_violin <- p_violin + coord_cartesian(ylim = y_limits)
  p_model  <- p_model  + coord_cartesian(ylim = y_limits)
  
  if (!is.null(y_scale_obj)) {
    p_violin <- p_violin + y_scale_obj
    p_model  <- p_model  + y_scale_obj
  }
  
  p_combined <- p_violin + p_model + patchwork::plot_layout(widths = c(1.3, 1))
  print(p_combined)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_", PSS_var,
             if (log10_transform) "_log10" else "",
             "_violin_plusMixedModel_depthlt", depth_cut, "_", string_ID, ".png")
    )
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p_combined, width = width, height = height, dpi = dpi, bg = "white")
    message("Saved plot to: ", out_path)
  }
  
  if (save_model) {
    if (is.null(model_prefix)) {
      model_prefix <- paste0(
        "SPECTRA_", PSS_var,
        if (log10_transform) "_log10" else "",
        "_mixedModel_depthlt", depth_cut, "_", string_ID
      )
    }
    
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    
    rds_path <- file.path(model_dir, paste0(model_prefix, ".rds"))
    save_obj <- list(
      model = m1,
      formula = fml,
      df_mod = df_mod,
      emmeans_raw = pred_df,
      emmeans_plot = pred_df2,
      seg_df = seg_df
    )
    saveRDS(save_obj, rds_path)
    
    sum_path   <- file.path(model_dir, paste0(model_prefix, "_summary.txt"))
    anova_path <- file.path(model_dir, paste0(model_prefix, "_anova.txt"))
    emm_path   <- file.path(model_dir, paste0(model_prefix, "_emmeans.csv"))
    
    capture.output(summary(m1), file = sum_path)
    capture.output(anova(m1), file = anova_path)
    utils::write.csv(pred_df, emm_path, row.names = FALSE)
  }
  
  invisible(list(
    plot = p_combined,
    model = m1,
    emmeans = pred_df2,
    df_mod = df_mod
  ))
}

run_mm_set <- function(dat, string_ID) {
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "slope",
    depth_cut = 100,
    output_dir = output_dir_figs_NB,
    y_limits = c(-2.5, 0),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "intercept",
    depth_cut = 100,
    output_dir = output_dir_figs_NB,
    y_limits = c(2, 20),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "total_biovolume",
    depth_cut = 100,
    log10_transform = TRUE,
    output_dir = output_dir_figs_NB,
    y_limits = c(7, 19),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "integrated_biovol_all",
    depth_cut = 100,
    log10_transform = TRUE,
    output_dir = output_dir_figs_NB,
    y_limits = c(7, 19),
    save_plot = TRUE,
    save_model = TRUE
  )
  
  PSS_violin_plus_mixedModel(
    dat = dat,
    string_ID = string_ID,
    PSS_var = "integrated_biovol_sample",
    depth_cut = 100,
    log10_transform = TRUE,
    output_dir = output_dir_figs_NB,
    y_limits = c(7, 19),
    save_plot = TRUE,
    save_model = TRUE
  )
}


prep_mm_for_compare <- function(dat,
                                string_ID,
                                PSS_var = "slope",
                                depth_cut = 100,
                                log10_transform = FALSE,
                                random_effect = "(1 | Cruise/Transect_ID)",
                                include_covariates = TRUE,
                                gap = 0.35) {
  
  needed <- c(PSS_var, "season", "shelf_loc_200", "depth_mid",
              "Cruise", "Transect_ID", "group_id", "distance_sample_km")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) stop("Missing columns in dat: ", paste(miss, collapse = ", "))
  
  yvar <- sym(PSS_var)
  
  df <- dat %>%
    filter(.data$depth_mid < depth_cut) %>%
    mutate(
      season = factor(.data$season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      Cruise = factor(.data$Cruise),
      Transect_ID = factor(.data$Transect_ID),
      .y_raw = suppressWarnings(as.numeric(as.character(!!yvar)))
    ) %>%
    filter(!is.na(.data$season), !is.na(.data$shelf_loc_200), is.finite(.data$.y_raw))
  
  if (nrow(df) == 0) stop("After filtering depth/NA/non-finite, no rows remain.")
  
  if (log10_transform) {
    df <- df %>%
      mutate(.y = if_else(.data$.y_raw > 0, log10(.data$.y_raw), NA_real_)) %>%
      filter(is.finite(.data$.y))
  } else {
    df <- df %>% mutate(.y = .data$.y_raw)
  }
  
  df_mod <- df %>%
    group_by(.data$Cruise, .data$Transect_ID, .data$group_id, .data$season, .data$shelf_loc_200) %>%
    summarise(
      y = first(.data$.y),
      depth_mid = first(.data$depth_mid),
      distance_sample_km = first(.data$distance_sample_km),
      .groups = "drop"
    ) %>%
    filter(is.finite(.data$y))
  
  fixed_part <- if (include_covariates) {
    "y ~ season * shelf_loc_200 + scale(depth_mid) + scale(distance_sample_km)"
  } else {
    "y ~ season * shelf_loc_200"
  }
  fml <- stats::as.formula(paste(fixed_part, "+", random_effect))
  
  m1 <- lmerTest::lmer(fml, data = df_mod, REML = TRUE)
  
  at_list <- if (include_covariates) {
    list(`scale(depth_mid)` = 0, `scale(distance_sample_km)` = 0)
  } else {
    list()
  }
  
  emm <- emmeans::emmeans(m1, ~ season * shelf_loc_200, at = at_list)
  pred_df <- as.data.frame(emm)
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$asymp.LCL, upper.CL = .data$asymp.UCL)
  } else if (all(c("LCL", "UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      rename(lower.CL = .data$LCL, upper.CL = .data$UCL)
  }
  
  pred_df2 <- pred_df %>%
    mutate(
      shelf_loc_200 = factor(.data$shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      season = factor(.data$season, levels = c("Winter", "Summer")),
      x_base_comp = 1 + (as.numeric(.data$shelf_loc_200) - 1) * gap
    )
  
  list(
    df = df,
    df_mod = df_mod,
    model = m1,
    pred_df = pred_df,
    pred_df2 = pred_df2,
    string_ID = string_ID,
    PSS_var = PSS_var,
    depth_cut = depth_cut,
    log10_transform = log10_transform
  )
}

plot_compare_sources_row <- function(dat_all,
                                     dat_plankton,
                                     PSS_var = "slope",
                                     depth_cut = 100,
                                     log10_transform = FALSE,
                                     y_limits = NULL,
                                     season_colors = c("Winter" = "#0072B2",
                                                       "Summer" = "#E69F00"),
                                     source_colors = c("allParticles" = "black",
                                                       "plankton" = "gray"),
                                     source_shapes = c("allParticles" = 16,
                                                       "plankton" = 17),
                                     violin_dodge = 0.75,
                                     emm_dodge = 0.22,
                                     ci_cap_width = 0.03,
                                     gap = 0.35) {
  
  obj_all <- prep_mm_for_compare(
    dat = dat_all,
    string_ID = "ISIIS_full",
    PSS_var = PSS_var,
    depth_cut = depth_cut,
    log10_transform = log10_transform,
    gap = gap
  )
  
  obj_plk <- prep_mm_for_compare(
    dat = dat_plankton,
    string_ID = "ISIIS_plankton",
    PSS_var = PSS_var,
    depth_cut = depth_cut,
    log10_transform = log10_transform,
    gap = gap
  )
  
  df_both <- bind_rows(
    obj_all$df %>% mutate(data_source = "allParticles"),
    obj_plk$df %>% mutate(data_source = "plankton")
  ) %>%
    mutate(
      data_source = factor(data_source, levels = c("allParticles", "plankton")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  pred_both <- bind_rows(
    obj_all$pred_df2 %>% mutate(data_source = "allParticles"),
    obj_plk$pred_df2 %>% mutate(data_source = "plankton")
  ) %>%
    mutate(
      data_source = factor(data_source, levels = c("allParticles", "plankton")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      
      # two source-centered subgroups within each shelf location
      source_offset = case_when(
        data_source == "allParticles" ~ -0.09,
        data_source == "plankton"     ~  0.09,
        TRUE ~ 0
      ),
      
      # winter first, summer second within each source subgroup
      season_offset = case_when(
        season == "Winter" ~ -0.03,
        season == "Summer" ~  0.03,
        TRUE ~ 0
      ),
      
      x = x_base_comp + source_offset + season_offset
    )
  
  seg_both <- pred_both %>%
    select(.data$shelf_loc_200, .data$data_source, .data$season, .data$x, .data$emmean) %>%
    pivot_wider(
      names_from = .data$season,
      values_from = c(.data$x, .data$emmean)
    ) %>%
    transmute(
      shelf_loc_200,
      data_source,
      x = .data$x_Winter,
      xend = .data$x_Summer,
      y = .data$emmean_Winter,
      yend = .data$emmean_Summer
    )
  
  if (PSS_var %in% c("total_biovolume", "integrated_biovol_all", "integrated_biovol_sample")) {
    y_label_expr <- expression(Biovolume~(mu*m^3~m^-3))
  } else if (PSS_var == "intercept") {
    y_label_expr <- expression(Intercept~(mu*m^3~m^-3~mu*m^-3))
  } else if (PSS_var == "slope") {
    y_label_expr <- expression(Slope~(m^-3~mu*m^-3))
  } else {
    y_label_expr <- PSS_var
  }
  
  y_scale_obj <- if (log10_transform || PSS_var == "intercept") {
    ggplot2::scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", round(x, 2)))
    )
  } else {
    NULL
  }
  
  if (is.null(y_limits)) {
    y_limits <- range(
      df_both$.y,
      pred_both$lower.CL,
      pred_both$upper.CL,
      na.rm = TRUE
    )
  }
  
  p_violin <- ggplot(
    df_both,
    aes(
      x = .data$shelf_loc_200,
      y = .data$.y,
      group = interaction(.data$shelf_loc_200, .data$season, .data$data_source)
    )
  ) +
    geom_violin(
      aes(fill = .data$season, color = .data$data_source),
      position = position_dodge(width = violin_dodge),
      trim = FALSE,
      alpha = 0.65,
      linewidth = 0.3
    ) +
    geom_boxplot(
      aes(
        fill = .data$season,
        color = .data$data_source,
        group = interaction(.data$shelf_loc_200, .data$season, .data$data_source)
      ),
      position = position_dodge(width = violin_dodge),
      width = 0.14,
      outlier.shape = NA,
      alpha = 0.35,
      linewidth = 0.3
    ) +
    scale_x_discrete(
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      )
    ) +
    scale_fill_manual(values = season_colors, name = "Season") +
    #scale_color_manual(values = season_colors, name = "Season") +
    scale_color_manual(
      values = source_colors,
      name = "Source",
      labels = c(
        allParticles = "Seston",
        plankton = "Plankton"
      )
    ) +
    labs(
      #title = paste0(PSS_var, " — raw data (depth < ", depth_cut, " m)"),
      x = NULL,
      y = y_label_expr
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  p_emm <- ggplot(
    pred_both,
    aes(
      x = .data$x,
      y = .data$emmean,
      color = .data$season,
      shape = .data$data_source
    )
  ) +
    geom_segment(
      data = seg_both,
      aes(
        x = .data$x,
        xend = .data$xend,
        y = .data$y,
        yend = .data$yend,
        group = interaction(.data$shelf_loc_200, .data$data_source)
      ),
      inherit.aes = FALSE,
      color = "black",
      #linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_errorbar(
      aes(ymin = .data$lower.CL, ymax = .data$upper.CL),
      width = ci_cap_width,
      linewidth = 0.6
    ) +
    geom_point(size = 3) +
    scale_x_continuous(
      breaks = c(1, 1 + gap),
      labels = c("off shelf", "on shelf")
    ) +
    scale_color_manual(values = season_colors, name = "Season") +
    scale_shape_manual(
      values = source_shapes,
      name = "Source",
      labels = c(
        allParticles = "Seston",
        plankton = "Plankton"
      )
    )+
    labs(
      #title = "Mixed model (emmeans ± CI)",
      x = NULL,
      y = NULL
    ) +
    coord_cartesian(ylim = y_limits) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18)
    )
  
  if (!is.null(y_scale_obj)) {
    p_violin <- p_violin + y_scale_obj
    p_emm    <- p_emm + y_scale_obj
  }
  
  p_violin + p_emm + patchwork::plot_layout(widths = c(1.3, 1))
}

make_allParticles_vs_plankton_3x2 <- function(dat_all,
                                              dat_plankton,
                                              output_dir = output_dir_figs_NB,
                                              width = 13,
                                              height = 16,
                                              dpi = 300) {
  
  row_slope <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "slope",
    depth_cut = 100,
    log10_transform = FALSE,
    y_limits = c(-2.5, 0)
  )
  
  row_intercept <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "intercept",
    depth_cut = 100,
    log10_transform = FALSE,
    y_limits = c(2, 20)
  )
  
  row_integrated <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "integrated_biovol_all",
    depth_cut = 100,
    log10_transform = TRUE,
    y_limits = c(7, 19)
  )
  
  p_final <- (row_slope / row_intercept / row_integrated) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(p_final)
  
  out_path <- file.path(
    output_dir,
    "SPECTRA_allParticles_vs_plankton_violin_plusMixedModel_3x2.png"
  )
  
  ggsave(
    out_path,
    p_final,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  message("Saved combined comparison figure to: ", out_path)
  
  invisible(p_final)
}

run_mm_set(dat_ISIIS_full,     "ISIIS_full")
run_mm_set(dat_ISIIS_plankton, "ISIIS_plankton")
run_mm_set(dat_ISIIS_detritus, "ISIIS_detritus")

make_allParticles_vs_plankton_3x2(
  dat_all = dat_ISIIS_full,
  dat_plankton = dat_ISIIS_plankton,
  output_dir = output_dir_figs_NB
)

# ============================================================
# 2x4 figure:
#   Y = log10(integrated_biovol_all)
#   X = slope
#   rows = Seston, Plankton
#   cols = w22, s22, w23, s23
#   separate linear fits for shelf_loc_200
# Also saves regression summary table.
# ============================================================
fit_biovol_vs_slope_models <- function(reg_dat) {
  
  split_dat <- reg_dat %>%
    group_by(source, Cruise, shelf_loc_200) %>%
    group_split()
  
  fit_tbl <- purrr::map_dfr(split_dat, function(df_sub) {
    
    source_i <- unique(df_sub$source)[1]
    cruise_i <- unique(df_sub$Cruise)[1]
    shelf_i  <- unique(df_sub$shelf_loc_200)[1]
    n_i      <- nrow(df_sub)
    
    if (n_i < 3 || dplyr::n_distinct(df_sub$slope) < 2) {
      return(tibble(
        source = source_i,
        Cruise = cruise_i,
        shelf_loc_200 = shelf_i,
        n = n_i,
        intercept = NA_real_,
        slope_coef = NA_real_,
        r_squared = NA_real_,
        f_statistic = NA_real_,
        df1 = NA_real_,
        df2 = NA_real_,
        p_value = NA_real_,
        equation_label = "Insufficient data"
      ))
    }
    
    fit <- lm(log10_integrated_biovol_all ~ slope, data = df_sub)
    smry <- summary(fit)
    
    intercept_i <- unname(coef(fit)[1])
    slope_i     <- unname(coef(fit)[2])
    r2_i        <- smry$r.squared
    
    f_i   <- unname(smry$fstatistic[1])
    df1_i <- unname(smry$fstatistic[2])
    df2_i <- unname(smry$fstatistic[3])
    p_i   <- pf(f_i, df1_i, df2_i, lower.tail = FALSE)
    
    eq_lab <- paste0(
      ifelse(shelf_i == "off_shelf", "Off shelf", "On shelf"),
      ": y = ",
      formatC(intercept_i, format = "f", digits = 2),
      ifelse(slope_i >= 0, " + ", " - "),
      formatC(abs(slope_i), format = "f", digits = 2),
      "x\nR² = ",
      formatC(r2_i, format = "f", digits = 2),
      ", F = ",
      formatC(f_i, format = "f", digits = 1),
      ", p = ",
      format.pval(p_i, digits = 2, eps = 1e-4)
    )
    
    tibble(
      source = source_i,
      Cruise = cruise_i,
      shelf_loc_200 = shelf_i,
      n = n_i,
      intercept = intercept_i,
      slope_coef = slope_i,
      r_squared = r2_i,
      f_statistic = f_i,
      df1 = df1_i,
      df2 = df2_i,
      p_value = p_i,
      equation_label = eq_lab
    )
  })
  
  fit_tbl
}

make_annotation_positions <- function(df_panel,
                                      fit_panel,
                                      x_limits = c(-2, NA),
                                      y_limits = c(NA, 20)) {
  
  x_data_min <- min(df_panel$slope, na.rm = TRUE)
  x_data_max <- max(df_panel$slope, na.rm = TRUE)
  y_data_min <- min(df_panel$log10_integrated_biovol_all, na.rm = TRUE)
  y_data_max <- max(df_panel$log10_integrated_biovol_all, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_max - 0.05 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_max - 0.22 * (y_max - y_min),
        TRUE ~ y_max
      )
    )
}

plot_biovol_vs_slope_panel <- function(df_panel,
                                       fit_tbl,
                                       source_i,
                                       cruise_i,
                                       x_limits = c(-2, NA),
                                       y_limits = c(NA, 20),
                                       point_size = 0.8,
                                       point_alpha = 0.18,
                                       show_x_title = FALSE,
                                       show_y_title = FALSE,
                                       show_legend = FALSE,
                                       strip_title = NULL,
                                       panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  p <- ggplot(
    df_panel,
    aes(x = slope, y = log10_integrated_biovol_all, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    scale_y_continuous(
      labels = function(x) parse(text = paste0("10^", x))
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(Slope~(m^-3~mu*m^-3)) else NULL,
      y = if (show_y_title) expression(Integrated~biovolume~(mu*m^3~m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

prep_biovol_vs_slope_data <- function(dat,
                                      source_name,
                                      depth_cut = 100) {
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    ) %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      integrated_biovol_all = first(integrated_biovol_all),
      .groups = "drop"
    ) %>%
    mutate(
      source = source_name,
      log10_integrated_biovol_all = if_else(
        is.finite(integrated_biovol_all) & integrated_biovol_all > 0,
        log10(integrated_biovol_all),
        NA_real_
      )
    ) %>%
    filter(
      is.finite(slope),
      is.finite(log10_integrated_biovol_all)
    )
}

plot_biovol_vs_slope_2x4 <- function(dat_all,
                                     dat_plankton,
                                     depth_cut = 100,
                                     output_dir = output_dir_figs_NB,
                                     width = 16,
                                     height = 8,
                                     dpi = 300,
                                     save_plot = TRUE,
                                     save_table = TRUE) {
  
  reg_all <- prep_biovol_vs_slope_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_biovol_vs_slope_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for biovolume-vs-slope plotting after filtering.")
  }
  
  fit_tbl <- fit_biovol_vs_slope_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_biovol_vs_slope_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(-2, NA)
  ylim_default <- c(NA, 20)
  
  xlim_s23_plk <- c(NA, NA)
  ylim_s23_plk <- c(NA, NA)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_s23_plk,
    y_limits = ylim_s23_plk,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_integrated_biovol_all_vs_slope_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved biovolume-vs-slope figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      base_dir,
      "tables",
      paste0("SPECTRA_integrated_biovol_all_vs_slope_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}

biovol_slope_out <- plot_biovol_vs_slope_2x4(
  dat_all = dat_ISIIS_full,
  dat_plankton = dat_ISIIS_plankton,
  depth_cut = 100,
  output_dir = output_dir_figs_NB,
  save_plot = TRUE,
  save_table = TRUE
)

#Table S6:
suppressPackageStartupMessages({
  library(gt)
})

# ------------------------------------------------
# Collapse to group-level observations
# ------------------------------------------------

prep_group_level <- function(dat, source_name, depth_cut = 100) {
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      intercept_untrans = 10^intercept) %>%
    filter(
      !is.na(season),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      intercept = first(intercept_untrans),
      integrated_biovolume = first(integrated_biovol_all),
      .groups = "drop"
    ) %>%
    mutate(source = source_name)
  
}

group_all <- prep_group_level(dat_ISIIS_full, "Seston")
group_plk <- prep_group_level(dat_ISIIS_plankton, "Plankton")

group_all_data <- bind_rows(group_all, group_plk)

# ------------------------------------------------
# Summary statistics
# ------------------------------------------------

summary_stats <- group_all_data %>%
  group_by(source, season, shelf_loc_200) %>%
  summarise(
    
    slope_med = median(slope, na.rm = TRUE),
    slope_sd  = sd(slope, na.rm = TRUE),
    
    int_med = median(intercept, na.rm = TRUE),
    int_q1  = quantile(intercept, 0.25, na.rm = TRUE, names = FALSE),
    int_q3  = quantile(intercept, 0.75, na.rm = TRUE, names = FALSE),
    
    bio_med = median(integrated_biovolume, na.rm = TRUE),
    bio_q1  = quantile(integrated_biovolume, 0.25, na.rm = TRUE, names = FALSE),
    bio_q3  = quantile(integrated_biovolume, 0.75, na.rm = TRUE, names = FALSE),
    
    .groups = "drop"
  )

# ------------------------------------------------
# Formatting helpers
# ------------------------------------------------

fmt_slope <- function(med, sd) {
  if (!is.finite(med)) return(NA_character_)
  if (!is.finite(sd)) {
    return(formatC(med, format = "f", digits = 2))
  }
  paste0(
    formatC(med, format = "f", digits = 2),
    " ± ",
    formatC(sd, format = "f", digits = 2)
  )
}

fmt_range <- function(med, q1, q3) {
  if (!is.finite(med)) return(NA_character_)
  if (!is.finite(q1) || !is.finite(q3)) {
    return(formatC(med, format = "e", digits = 2))
  }
  paste0(
    formatC(med, format = "e", digits = 2),
    " [",
    formatC(q1, format = "e", digits = 2),
    "–",
    formatC(q3, format = "e", digits = 2),
    "]"
  )
}


# ------------------------------------------------
# Convert to long metric rows
# ------------------------------------------------

table_long <- summary_stats %>%
  mutate(
    slope = mapply(fmt_slope, slope_med, slope_sd),
    intercept = mapply(fmt_range, int_med, int_q1, int_q3),
    integrated_biovolume = mapply(fmt_range, bio_med, bio_q1, bio_q3)
  ) %>%
  select(source, season, shelf_loc_200, slope, intercept, integrated_biovolume) %>%
  pivot_longer(
    cols = c(slope, intercept, integrated_biovolume),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    source = factor(source, levels = c("Seston", "Plankton")),
    metric = factor(metric, levels = c("slope", "intercept", "integrated_biovolume")),
    row_name = paste(
      recode(metric,
             slope = "Slope",
             intercept = "Intercept",
             integrated_biovolume = "Integrated biovolume"
      ),
      source,
      sep = " – "
    ),
    column_name = paste(
      as.character(season),
      ifelse(shelf_loc_200 == "off_shelf", "off shelf", "on shelf"),
      sep = " – "
    )
  )

# ------------------------------------------------
# Final wide table, ordered by metric first
# ------------------------------------------------

table_wide <- table_long %>%
  mutate(
    metric_order = factor(metric, levels = c("slope", "intercept", "integrated_biovolume")),
    source_order = factor(source, levels = c("Seston", "Plankton"))
  ) %>%
  arrange(metric_order, source_order) %>%
  select(row_name, column_name, value, metric_order, source_order) %>%
  pivot_wider(
    names_from = column_name,
    values_from = value
  ) %>%
  arrange(metric_order, source_order) %>%
  select(
    row_name,
    `Winter – off shelf`,
    `Winter – on shelf`,
    `Summer – off shelf`,
    `Summer – on shelf`
  )

# ------------------------------------------------
# Save CSV
# ------------------------------------------------

dir.create(file.path(base_dir,"tables"), showWarnings = FALSE)

readr::write_csv(
  table_wide,
  file.path(base_dir, "tables", "season_shelf_PSS_coefficients_table.csv"),
  na = ""
)
# ------------------------------------------------
# GT formatted table
# ------------------------------------------------

gt_table <- table_wide %>%
  gt(rowname_col = "row_name") %>%
  cols_align("center", everything()) %>%
  tab_options(
    table.font.size = px(13),
    column_labels.font.weight = "bold"
  )

gtsave(
  gt_table,
  file.path(base_dir,"tables","season_shelf_PSS_coefficients_table.html")
)

gt_table


# ------------------------------------------------
# Chla and slope
# ------------------------------------------------
prep_slope_vs_chla_data <- function(dat,
                                    source_name,
                                    depth_cut = 100) {
  
  needed <- c("slope", "chl_a_ug_l", "Cruise", "Transect_ID", "group_id", "season", "shelf_loc_200")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns for slope vs chl_a_ug_l plot: ",
         paste(miss, collapse = ", "))
  }
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      slope = as.numeric(slope),
      chl_a_ug_l = as.numeric(chl_a_ug_l)
    ) %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      chl_a_ug_l = first(chl_a_ug_l),
      .groups = "drop"
    ) %>%
    mutate(
      source = source_name,
      log10_chl_a_ug_l = if_else(
        is.finite(chl_a_ug_l) & chl_a_ug_l > 0,
        log10(chl_a_ug_l),
        NA_real_
      )
    ) %>%
    filter(
      is.finite(slope),
      is.finite(log10_chl_a_ug_l)
    )
}

fit_slope_vs_chla_models <- function(reg_dat) {
  
  split_dat <- reg_dat %>%
    group_by(source, Cruise, shelf_loc_200) %>%
    group_split()
  
  fit_tbl <- purrr::map_dfr(split_dat, function(df_sub) {
    
    source_i <- unique(df_sub$source)[1]
    cruise_i <- unique(df_sub$Cruise)[1]
    shelf_i  <- unique(df_sub$shelf_loc_200)[1]
    n_i      <- nrow(df_sub)
    
    if (n_i < 3 || dplyr::n_distinct(df_sub$log10_chl_a_ug_l) < 2) {
      return(tibble(
        source = source_i,
        Cruise = cruise_i,
        shelf_loc_200 = shelf_i,
        n = n_i,
        intercept = NA_real_,
        slope_coef = NA_real_,
        r_squared = NA_real_,
        f_statistic = NA_real_,
        df1 = NA_real_,
        df2 = NA_real_,
        p_value = NA_real_,
        equation_label = "Insufficient data"
      ))
    }
    
    fit <- lm(slope ~ log10_chl_a_ug_l, data = df_sub)
    smry <- summary(fit)
    
    intercept_i <- unname(coef(fit)[1])
    slope_i     <- unname(coef(fit)[2])
    r2_i        <- smry$r.squared
    
    f_i   <- unname(smry$fstatistic[1])
    df1_i <- unname(smry$fstatistic[2])
    df2_i <- unname(smry$fstatistic[3])
    p_i   <- pf(f_i, df1_i, df2_i, lower.tail = FALSE)
    
    eq_lab <- paste0(
      ifelse(shelf_i == "off_shelf", "Off shelf", "On shelf"),
      ": y = ",
      formatC(intercept_i, format = "f", digits = 2),
      ifelse(slope_i >= 0, " + ", " - "),
      formatC(abs(slope_i), format = "f", digits = 2),
      "x\nR² = ",
      formatC(r2_i, format = "f", digits = 2),
      ", F = ",
      formatC(f_i, format = "f", digits = 1),
      ", p = ",
      format.pval(p_i, digits = 2, eps = 1e-4)
    )
    
    tibble(
      source = source_i,
      Cruise = cruise_i,
      shelf_loc_200 = shelf_i,
      n = n_i,
      intercept = intercept_i,
      slope_coef = slope_i,
      r_squared = r2_i,
      f_statistic = f_i,
      df1 = df1_i,
      df2 = df2_i,
      p_value = p_i,
      equation_label = eq_lab
    )
  })
  
  fit_tbl
}

make_annotation_positions_slope_chla <- function(df_panel,
                                                 fit_panel,
                                                 x_limits = c(NA, NA),
                                                 y_limits = c(-2.5, 0)) {
  
  x_data_min <- min(df_panel$log10_chl_a_ug_l, na.rm = TRUE)
  x_data_max <- max(df_panel$log10_chl_a_ug_l, na.rm = TRUE)
  y_data_min <- min(df_panel$slope, na.rm = TRUE)
  y_data_max <- max(df_panel$slope, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_min + 0.28 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_min + 0.14 * (y_max - y_min),
        TRUE ~ y_min
      )
    )
}

plot_slope_vs_chla_panel <- function(df_panel,
                                     fit_tbl,
                                     source_i,
                                     cruise_i,
                                     x_limits = c(NA, NA),
                                     y_limits = c(-2.5, 0),
                                     point_size = 0.8,
                                     point_alpha = 0.18,
                                     show_x_title = FALSE,
                                     show_y_title = FALSE,
                                     show_legend = FALSE,
                                     strip_title = NULL,
                                     panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions_slope_chla(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  p <- ggplot(
    df_panel,
    aes(x = log10_chl_a_ug_l, y = slope, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(log[10]~chl[a]~(mu*g~L^-1)) else NULL,
      y = if (show_y_title) expression(Slope~(m^-3~mu*m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_slope_vs_chla_2x4 <- function(dat_all,
                                   dat_plankton,
                                   depth_cut = 100,
                                   output_dir = output_dir_figs_NB,
                                   width = 16,
                                   height = 8,
                                   dpi = 300,
                                   save_plot = TRUE,
                                   save_table = TRUE) {
  
  reg_all <- prep_slope_vs_chla_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_slope_vs_chla_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for slope-vs-chla plotting after filtering.")
  }
  
  fit_tbl <- fit_slope_vs_chla_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_slope_vs_chla_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(NA, NA)
  ylim_default <- c(-2.5, 0)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_slope_vs_log10_chla_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved slope-vs-chla figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      base_dir,
      "tables",
      paste0("SPECTRA_slope_vs_log10_chla_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}

slope_chla_out <- plot_slope_vs_chla_2x4(
  dat_all = dat_ISIIS_full,
  dat_plankton = dat_ISIIS_plankton,
  depth_cut = 100,
  output_dir = output_dir_figs_NB,
  save_plot = TRUE,
  save_table = TRUE
)

# ------------------------------------------------
# temperature and slope
# ------------------------------------------------

prep_slope_vs_temperature_data <- function(dat,
                                           source_name,
                                           depth_cut = 100) {
  
  needed <- c("slope", "Temperature", "Cruise", "Transect_ID", "group_id", "season", "shelf_loc_200")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns for slope vs Temperature plot: ",
         paste(miss, collapse = ", "))
  }
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      slope = as.numeric(slope),
      Temperature = as.numeric(Temperature)
    ) %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      Temperature = first(Temperature),
      .groups = "drop"
    ) %>%
    mutate(source = source_name) %>%
    filter(
      is.finite(slope),
      is.finite(Temperature)
    )
}

fit_slope_vs_temperature_models <- function(reg_dat) {
  
  split_dat <- reg_dat %>%
    group_by(source, Cruise, shelf_loc_200) %>%
    group_split()
  
  fit_tbl <- purrr::map_dfr(split_dat, function(df_sub) {
    
    source_i <- unique(df_sub$source)[1]
    cruise_i <- unique(df_sub$Cruise)[1]
    shelf_i  <- unique(df_sub$shelf_loc_200)[1]
    n_i      <- nrow(df_sub)
    
    if (n_i < 3 || dplyr::n_distinct(df_sub$Temperature) < 2) {
      return(tibble(
        source = source_i,
        Cruise = cruise_i,
        shelf_loc_200 = shelf_i,
        n = n_i,
        intercept = NA_real_,
        slope_coef = NA_real_,
        r_squared = NA_real_,
        f_statistic = NA_real_,
        df1 = NA_real_,
        df2 = NA_real_,
        p_value = NA_real_,
        equation_label = "Insufficient data"
      ))
    }
    
    fit <- lm(slope ~ Temperature, data = df_sub)
    smry <- summary(fit)
    
    intercept_i <- unname(coef(fit)[1])
    slope_i     <- unname(coef(fit)[2])
    r2_i        <- smry$r.squared
    
    f_i   <- unname(smry$fstatistic[1])
    df1_i <- unname(smry$fstatistic[2])
    df2_i <- unname(smry$fstatistic[3])
    p_i   <- pf(f_i, df1_i, df2_i, lower.tail = FALSE)
    
    eq_lab <- paste0(
      ifelse(shelf_i == "off_shelf", "Off shelf", "On shelf"),
      ": y = ",
      formatC(intercept_i, format = "f", digits = 2),
      ifelse(slope_i >= 0, " + ", " - "),
      formatC(abs(slope_i), format = "f", digits = 2),
      "x\nR² = ",
      formatC(r2_i, format = "f", digits = 2),
      ", F = ",
      formatC(f_i, format = "f", digits = 1),
      ", p = ",
      format.pval(p_i, digits = 2, eps = 1e-4)
    )
    
    tibble(
      source = source_i,
      Cruise = cruise_i,
      shelf_loc_200 = shelf_i,
      n = n_i,
      intercept = intercept_i,
      slope_coef = slope_i,
      r_squared = r2_i,
      f_statistic = f_i,
      df1 = df1_i,
      df2 = df2_i,
      p_value = p_i,
      equation_label = eq_lab
    )
  })
  
  fit_tbl
}

make_annotation_positions_slope_temperature <- function(df_panel,
                                                        fit_panel,
                                                        x_limits = c(NA, NA),
                                                        y_limits = c(-2.5, 0)) {
  
  x_data_min <- min(df_panel$Temperature, na.rm = TRUE)
  x_data_max <- max(df_panel$Temperature, na.rm = TRUE)
  y_data_min <- min(df_panel$slope, na.rm = TRUE)
  y_data_max <- max(df_panel$slope, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_min + 0.28 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_min + 0.14 * (y_max - y_min),
        TRUE ~ y_min
      )
    )
}

plot_slope_vs_temperature_panel <- function(df_panel,
                                            fit_tbl,
                                            source_i,
                                            cruise_i,
                                            x_limits = c(NA, NA),
                                            y_limits = c(-2.5, 0),
                                            point_size = 0.8,
                                            point_alpha = 0.18,
                                            show_x_title = FALSE,
                                            show_y_title = FALSE,
                                            show_legend = FALSE,
                                            strip_title = NULL,
                                            panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions_slope_temperature(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  p <- ggplot(
    df_panel,
    aes(x = Temperature, y = slope, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(Temperature~(degree*C)) else NULL,
      y = if (show_y_title) expression(Slope~(m^-3~mu*m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_slope_vs_temperature_2x4 <- function(dat_all,
                                          dat_plankton,
                                          depth_cut = 100,
                                          output_dir = output_dir_figs_NB,
                                          width = 16,
                                          height = 8,
                                          dpi = 300,
                                          save_plot = TRUE,
                                          save_table = TRUE) {
  
  reg_all <- prep_slope_vs_temperature_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_slope_vs_temperature_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for slope-vs-temperature plotting after filtering.")
  }
  
  fit_tbl <- fit_slope_vs_temperature_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_slope_vs_temperature_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(NA, NA)
  ylim_default <- c(-2.5, 0)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_slope_vs_temperature_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved slope-vs-temperature figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      base_dir,
      "tables",
      paste0("SPECTRA_slope_vs_temperature_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}

slope_temperature_out <- plot_slope_vs_temperature_2x4(
  dat_all = dat_ISIIS_full,
  dat_plankton = dat_ISIIS_plankton,
  depth_cut = 100,
  output_dir = output_dir_figs_NB,
  save_plot = TRUE,
  save_table = TRUE
)

# ------------------------------------------------
# oxygen and slope
# ------------------------------------------------

prep_slope_vs_oxygen_data <- function(dat,
                                           source_name,
                                           depth_cut = 100) {
  
  needed <- c("slope", "Oxygen_ml_l", "Cruise", "Transect_ID", "group_id", "season", "shelf_loc_200")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns for slope vs oxygen plot: ",
         paste(miss, collapse = ", "))
  }
  
  dat %>%
    filter(depth_mid < depth_cut) %>%
    mutate(
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      season = factor(season, levels = c("Winter", "Summer")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf")),
      slope = as.numeric(slope),
      Oxygen_ml_l = as.numeric(Oxygen_ml_l)
    ) %>%
    dplyr::filter(.data$Cruise != "s23") %>%   # skip incorrect Summer 2023 oxygen %>%
    filter(
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(group_id),
      !is.na(shelf_loc_200)
    ) %>%
    group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
    summarise(
      slope = first(slope),
      Oxygen_ml_l = first(Oxygen_ml_l),
      .groups = "drop"
    ) %>%
    mutate(source = source_name) %>%
    filter(
      is.finite(slope),
      is.finite(Oxygen_ml_l)
    )
}

fit_slope_vs_oxygen_models <- function(reg_dat) {
  
  split_dat <- reg_dat %>%
    group_by(source, Cruise, shelf_loc_200) %>%
    group_split()
  
  fit_tbl <- purrr::map_dfr(split_dat, function(df_sub) {
    
    source_i <- unique(df_sub$source)[1]
    cruise_i <- unique(df_sub$Cruise)[1]
    shelf_i  <- unique(df_sub$shelf_loc_200)[1]
    n_i      <- nrow(df_sub)
    
    if (n_i < 3 || dplyr::n_distinct(df_sub$Oxygen_ml_l) < 2) {
      return(tibble(
        source = source_i,
        Cruise = cruise_i,
        shelf_loc_200 = shelf_i,
        n = n_i,
        intercept = NA_real_,
        slope_coef = NA_real_,
        r_squared = NA_real_,
        f_statistic = NA_real_,
        df1 = NA_real_,
        df2 = NA_real_,
        p_value = NA_real_,
        equation_label = "Insufficient data"
      ))
    }
    
    fit <- lm(slope ~ Oxygen_ml_l, data = df_sub)
    smry <- summary(fit)
    
    intercept_i <- unname(coef(fit)[1])
    slope_i     <- unname(coef(fit)[2])
    r2_i        <- smry$r.squared
    
    f_i   <- unname(smry$fstatistic[1])
    df1_i <- unname(smry$fstatistic[2])
    df2_i <- unname(smry$fstatistic[3])
    p_i   <- pf(f_i, df1_i, df2_i, lower.tail = FALSE)
    
    eq_lab <- paste0(
      ifelse(shelf_i == "off_shelf", "Off shelf", "On shelf"),
      ": y = ",
      formatC(intercept_i, format = "f", digits = 2),
      ifelse(slope_i >= 0, " + ", " - "),
      formatC(abs(slope_i), format = "f", digits = 2),
      "x\nR² = ",
      formatC(r2_i, format = "f", digits = 2),
      ", F = ",
      formatC(f_i, format = "f", digits = 1),
      ", p = ",
      format.pval(p_i, digits = 2, eps = 1e-4)
    )
    
    tibble(
      source = source_i,
      Cruise = cruise_i,
      shelf_loc_200 = shelf_i,
      n = n_i,
      intercept = intercept_i,
      slope_coef = slope_i,
      r_squared = r2_i,
      f_statistic = f_i,
      df1 = df1_i,
      df2 = df2_i,
      p_value = p_i,
      equation_label = eq_lab
    )
  })
  
  fit_tbl
}

make_annotation_positions_slope_oxygen <- function(df_panel,
                                                        fit_panel,
                                                        x_limits = c(NA, NA),
                                                        y_limits = c(-2.5, 0)) {
  
  x_data_min <- min(df_panel$Oxygen_ml_l, na.rm = TRUE)
  x_data_max <- max(df_panel$Oxygen_ml_l, na.rm = TRUE)
  y_data_min <- min(df_panel$slope, na.rm = TRUE)
  y_data_max <- max(df_panel$slope, na.rm = TRUE)
  
  x_min <- if (!is.na(x_limits[1])) x_limits[1] else x_data_min
  x_max <- if (!is.na(x_limits[2])) x_limits[2] else x_data_max
  y_min <- if (!is.na(y_limits[1])) y_limits[1] else y_data_min
  y_max <- if (!is.na(y_limits[2])) y_limits[2] else y_data_max
  
  fit_panel %>%
    mutate(
      x = x_max - 0.04 * (x_max - x_min),
      y = case_when(
        shelf_loc_200 == "on_shelf"  ~ y_min + 0.28 * (y_max - y_min),
        shelf_loc_200 == "off_shelf" ~ y_min + 0.14 * (y_max - y_min),
        TRUE ~ y_min
      )
    )
}

plot_slope_vs_oxygen_panel <- function(df_panel,
                                            fit_tbl,
                                            source_i,
                                            cruise_i,
                                            x_limits = c(NA, NA),
                                            y_limits = c(-2.5, 0),
                                            point_size = 0.8,
                                            point_alpha = 0.18,
                                            show_x_title = FALSE,
                                            show_y_title = FALSE,
                                            show_legend = FALSE,
                                            strip_title = NULL,
                                            panel_letter = NULL) {
  
  fit_panel <- fit_tbl %>%
    filter(source == source_i, Cruise == cruise_i)
  
  ann_panel <- make_annotation_positions_slope_oxygen(
    df_panel = df_panel,
    fit_panel = fit_panel,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  p <- ggplot(
    df_panel,
    aes(x = Oxygen_ml_l, y = slope, color = shelf_loc_200)
  ) +
    geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9
    ) +
    geom_text(
      data = ann_panel,
      aes(x = x, y = y, label = equation_label, color = shelf_loc_200),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.6,
      lineheight = 1.05,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf  = "#D55E00"
      ),
      labels = c(
        off_shelf = "off shelf",
        on_shelf  = "on shelf"
      ),
      name = "Shelf location"
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      title = strip_title,
      x = if (show_x_title) expression(Oxygen~(mL~L^-1)) else NULL,
      y = if (show_y_title) expression(Slope~(m^-3~mu*m^-3)) else NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.direction = "horizontal",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12)
    )
  
  if (!is.null(panel_letter)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_letter,
        hjust = -0.3, vjust = 1.4,
        size = 6, fontface = "bold",
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_slope_vs_oxygen_2x4 <- function(dat_all,
                                          dat_plankton,
                                          depth_cut = 100,
                                          output_dir = output_dir_figs_NB,
                                          width = 16,
                                          height = 8,
                                          dpi = 300,
                                          save_plot = TRUE,
                                          save_table = TRUE) {
  
  reg_all <- prep_slope_vs_oxygen_data(
    dat = dat_all,
    source_name = "Seston",
    depth_cut = depth_cut
  )
  
  reg_plk <- prep_slope_vs_oxygen_data(
    dat = dat_plankton,
    source_name = "Plankton",
    depth_cut = depth_cut
  )
  
  reg_dat <- bind_rows(reg_all, reg_plk) %>%
    mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
  
  if (nrow(reg_dat) == 0) {
    stop("No rows available for slope-vs-oxygen plotting after filtering.")
  }
  
  fit_tbl <- fit_slope_vs_oxygen_models(reg_dat)
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  get_panel_data <- function(src, cru) {
    reg_dat %>% filter(source == src, Cruise == cru)
  }
  
  make_panel <- function(src,
                         cru,
                         x_limits,
                         y_limits,
                         show_x_title = FALSE,
                         show_y_title = FALSE,
                         show_legend = FALSE,
                         strip_title = NULL,
                         panel_letter = NULL) {
    
    df_panel <- get_panel_data(src, cru)
    
    if (nrow(df_panel) == 0) {
      return(
        ggplot() +
          theme_void() +
          labs(title = strip_title) +
          annotate("text", x = 0.5, y = 0.5, label = "", size = 5) +
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
          )
      )
    }
    
    plot_slope_vs_oxygen_panel(
      df_panel = df_panel,
      fit_tbl = fit_tbl,
      source_i = src,
      cruise_i = cru,
      x_limits = x_limits,
      y_limits = y_limits,
      point_size = 0.8,
      point_alpha = 0.18,
      show_x_title = show_x_title,
      show_y_title = show_y_title,
      show_legend = show_legend,
      strip_title = strip_title,
      panel_letter = panel_letter
    )
  }
  
  xlim_default <- c(NA, NA)
  ylim_default <- c(-2.5, 0)
  
  p11 <- make_panel(
    src = "Seston",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_legend = TRUE,
    show_y_title = FALSE,
    strip_title = cruise_labels[["w22"]],
    panel_letter = "a."
  )
  
  p12 <- make_panel(
    src = "Seston",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s22"]],
    panel_letter = "b."
  )
  
  p13 <- make_panel(
    src = "Seston",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["w23"]],
    panel_letter = "c."
  )
  
  p14 <- make_panel(
    src = "Seston",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    strip_title = cruise_labels[["s23"]],
    panel_letter = "d."
  )
  
  p21 <- make_panel(
    src = "Plankton",
    cru = "w22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    show_y_title = TRUE,
    panel_letter = "e."
  )
  
  p22 <- make_panel(
    src = "Plankton",
    cru = "s22",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = TRUE,
    panel_letter = "f."
  )
  
  p23 <- make_panel(
    src = "Plankton",
    cru = "w23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "g."
  )
  
  p24 <- make_panel(
    src = "Plankton",
    cru = "s23",
    x_limits = xlim_default,
    y_limits = ylim_default,
    show_x_title = FALSE,
    panel_letter = "h."
  )
  
  top_row <- p11 | p12 | p13 | p14
  bottom_row <- p21 | p22 | p23 | p24
  
  p_final <- (top_row / bottom_row) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      subtitle = paste0("Depth < ", depth_cut, " m")
    ) &
    theme(legend.position = "top")
  
  print(p_final)
  
  if (save_plot) {
    out_path <- file.path(
      output_dir,
      paste0("SPECTRA_slope_vs_oxygen_2x4_depthlt", depth_cut, ".png")
    )
    
    ggsave(
      out_path,
      p_final,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    message("Saved slope-vs-oxygen figure to: ", out_path)
  }
  
  if (save_table) {
    dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    table_path <- file.path(
      base_dir,
      "tables",
      paste0("SPECTRA_slope_vs_oxygen_regression_table_depthlt", depth_cut, ".csv")
    )
    
    readr::write_csv(fit_tbl, table_path)
    message("Saved regression table to: ", table_path)
  }
  
  invisible(list(
    plot = p_final,
    regression_table = fit_tbl,
    plot_data = reg_dat
  ))
}


slope_oxygen_out <- plot_slope_vs_oxygen_2x4(
  dat_all = dat_ISIIS_full,
  dat_plankton = dat_ISIIS_plankton,
  depth_cut = 100,
  output_dir = output_dir_figs_NB,
  save_plot = TRUE,
  save_table = TRUE
)


