#!/usr/bin/env Rscript
# SPECTRA: self-contained mixed-model analysis
# Run from the SPECTRA project root, e.g. Rscript SCRIPTS/<script>.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lubridate)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(patchwork)
  library(optparse)
  library(rlang)
  library(multcompView)
})
options <- list(
  optparse::make_option('--thinned', action = 'store_true', default = FALSE,
                        help = 'Thin on-shelf casts to every third cast')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = options))
setwd('/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
base_dir <- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'model_functions.R'))
sizeSpectraDirPlankton <- file.path(base_dir, 'data', 'size_spectra_plankton_5-10m-bins_thresholded')
sizeSpectraDirAll <- file.path(base_dir, 'data', 'ISIIS-LISST_mergedNBSS_5m_allParticles')
df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = 'by_community\\.csv$', full.names = TRUE, recursive = TRUE)
df_list_full <- list.files(sizeSpectraDirAll, pattern = 'LISST-ISIIS_NBSS\\.csv$', full.names = TRUE, recursive = TRUE)
env_vars <- c('Temperature', 'Salinity', 'Oxygen_ml_l', 'pH', 'PAR', 'chl_a_ug_l')
if (!length(df_list_full) || !length(df_list_plankton)) stop('No input files found. Run from SPECTRA project root.')
load_ISIIS_data <- function(df_list, env_vars = c("Temperature","Salinity","Oxygen_ml_l","pH","PAR","chl_a_ug_l"), x1_full = 12.78849, x2_full = 4.547442e13) {
  
  dat_ISIIS <- purrr::map2_dfr(
    df_list,
    seq_along(df_list),
    ~ readr::read_csv(
      .x,
      show_col_types = FALSE,
      col_select = c(
        "group_id",
        dplyr::any_of(c("total_biovolume", "biovolume_um3")),
        "biovol_um3_mid", "biovol_um3_width",
        "time_Pacific", "rel_dist", "depth_mid", "Latitude", "Longitude",
        "slope", "intercept", "Cruise", "Transect_ID", "cast_id",
        "logNB", "ECD_um_mid",
        dplyr::any_of(c("grouped_taxon", "season", "distance_sample_km")),
        dplyr::any_of(env_vars)
      )
    ) %>%
      dplyr::mutate(
        source_file = .x,
        source_id = .y
      )
  )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(time_Pacific = lubridate::ymd_hms(time_Pacific,quiet = TRUE)) %>%
    dplyr::arrange(Cruise,Transect_ID,cast_id,time_Pacific)
  
  cast_types <- dat_ISIIS %>%
    dplyr::select(Cruise, Transect_ID,cast_id,time_Pacific,depth_mid) %>%
    dplyr::filter(!is.na(Cruise),!is.na(Transect_ID),!is.na(cast_id),!is.na(time_Pacific),!is.na(depth_mid)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Cruise,Transect_ID,cast_id,time_Pacific) %>%
    dplyr::group_by(Cruise,Transect_ID,cast_id) %>%
    dplyr::summarise(
      n_depths = dplyr::n_distinct(depth_mid),
      first_depth = dplyr::first(depth_mid),
      last_depth = dplyr::last(depth_mid),
      depth_change = last_depth - first_depth,
      cast_type = dplyr::case_when(
        n_depths == 1 ~ "one_depth",
        depth_change > 0 ~ "down_cast",
        depth_change < 0 ~ "up_cast",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    )
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::left_join(cast_types %>%
    dplyr::select(Cruise, Transect_ID, cast_id, cast_type),
      by = c("Cruise","Transect_ID","cast_id"))
  
  if (!"total_biovolume" %in% names(dat_ISIIS)) dat_ISIIS$total_biovolume <- NA_real_
  if (!"biovolume_um3"  %in% names(dat_ISIIS)) dat_ISIIS$biovolume_um3  <- NA_real_
  if (!"grouped_taxon"  %in% names(dat_ISIIS)) dat_ISIIS$grouped_taxon  <- NA_character_
  if (!"season"         %in% names(dat_ISIIS)) dat_ISIIS$season         <- NA_character_
  if (!"distance_sample_km" %in% names(dat_ISIIS)) dat_ISIIS$distance_sample_km <- NA_real_
  
  dat_ISIIS <- dat_ISIIS %>%
    mutate(
      total_biovolume = coalesce(total_biovolume, biovolume_um3)
    ) %>%
    dplyr::select(-any_of("biovolume_um3")) %>%
    mutate(
      time_Pacific_orig = time_Pacific,
      mean_time_Pacific = time_Pacific - hours(1)
    ) %>%
    arrange(time_Pacific) %>%
    rename(group_id_old = group_id) %>%
    mutate(
      group_id = dplyr::cur_group_id(), .by = c(source_id, group_id_old),
      total_biovolume  = as.numeric(total_biovolume),
      biovolume_um3_m3 = as.numeric((10^logNB)*biovol_um3_width),
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
      total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE),
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
    dplyr::select(
      group_id,
      total_biovolume_sum,
      total_biovolume_um3_m3,
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

remove_duplicate_NBSS_samples <- function(dat, verbose = TRUE) {
  
  # ------------------------------------------------------------
  # 1. Create sample-level dataset
  # ------------------------------------------------------------
  
  sample_data <- dat %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(
      Cruise = dplyr::first(Cruise),
      Transect_ID = dplyr::first(Transect_ID),
      cast_id = dplyr::first(cast_id),
      depth_mid = dplyr::first(depth_mid),
      Longitude = dplyr::first(Longitude),
      time_Pacific = dplyr::first(time_Pacific),
      slope = dplyr::first(slope),
      intercept = dplyr::first(intercept),
      total_biovolume_um3_m3 = dplyr::first(total_biovolume_um3_m3),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------
  # 2. Identify duplicated NBSS samples
  # ------------------------------------------------------------
  
  duplicate_samples <- sample_data %>%
    dplyr::filter(
      !is.na(depth_mid),
      !is.na(Longitude),
      !is.na(slope),
      !is.na(intercept),
      !is.na(total_biovolume_um3_m3)
    ) %>%
    dplyr::group_by(
      Cruise, Transect_ID, depth_mid, Longitude,
      slope, intercept, total_biovolume_um3_m3
    ) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::arrange(time_Pacific, group_id, .by_group = TRUE) %>%
    dplyr::mutate(
      duplicate_rank = dplyr::row_number(),
      action = dplyr::if_else(
        duplicate_rank == 1,
        "retain",
        "remove"
      )
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------
  # 3. Identify group_ids to remove
  # ------------------------------------------------------------
  
  groups_to_remove <- duplicate_samples %>%
    dplyr::filter(action == "remove") %>%
    dplyr::pull(group_id)
  
  # ------------------------------------------------------------
  # 4. Remove all size-bin rows for duplicated samples
  # ------------------------------------------------------------
  
  dat_clean <- dat %>%
    dplyr::filter(!group_id %in% groups_to_remove)
  
  # ------------------------------------------------------------
  # 5. Print summary
  # ------------------------------------------------------------
  
  if (verbose) {
    cat("\n========================================\n")
    cat("NBSS DUPLICATE SAMPLE CLEANING\n")
    cat("========================================\n")
    cat("Original samples:", dplyr::n_distinct(dat$group_id), "\n")
    cat("Duplicate sets:", sum(duplicate_samples$duplicate_rank == 1), "\n")
    cat("Samples removed:", length(groups_to_remove), "\n")
    cat("Remaining samples:", dplyr::n_distinct(dat_clean$group_id), "\n")
    cat("========================================\n\n")
  }
  
  return(list(
    data = dat_clean,
    duplicate_report = duplicate_samples,
    removed_group_ids = groups_to_remove
  ))

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
    stop("Missing required columns for mixed models: ", paste(miss, collapse = ", "))
  }
  
  if (all(is.na(dat_ISIIS$season))) {
    stop("Column 'season' exists but is all NA.")
  }
  if (all(is.na(dat_ISIIS$distance_sample_km))) {
    stop("Column 'distance_sample_km' exists but is all NA.")
  }
  
  return(dat_ISIIS)
  
}

  
load_clean <- function(files) {
  x <- load_ISIIS_data(files, env_vars = env_vars)
  x <- remove_duplicate_NBSS_samples(x)$data
  clean_ISIIS_dat(x)
}
dat_ISIIS_full <- load_clean(df_list_full)
dat_ISIIS_plankton <- load_clean(df_list_plankton)

dat_ISIIS_full <- overwrite_distance_from_shore(dat_ISIIS_full)
dat_ISIIS_plankton <- overwrite_distance_from_shore(dat_ISIIS_plankton)

identify_coarse_casts <- function(x) {
  
  x <- x %>%
    dplyr::arrange(time_Pacific) %>%
    dplyr::mutate(cast_position = dplyr::row_number(),
                  keep_cast = (cast_position - 1) %% 3 == 0)
  
  x
}
thin_onshelf_casts <- function(dat) {
  
  cast_sequence <- dat %>%
    dplyr::filter(shelf_loc_200 == "on_shelf", cast_type %in% c("down_cast", "up_cast")) %>%
    dplyr::group_by(Cruise, Transect_ID, cast_id, cast_type) %>%
    dplyr::summarise(time_Pacific = min(time_Pacific, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Cruise, Transect_ID, time_Pacific) %>%
    dplyr::group_by(Cruise, Transect_ID) %>%
    dplyr::group_modify(~ identify_coarse_casts(.x)) %>%
    dplyr::ungroup()
  
  casts_to_keep <- cast_sequence %>%
    dplyr::filter(keep_cast) %>%
    dplyr::select(Cruise, Transect_ID, cast_id) %>%
    dplyr::distinct()
  
  dat_thinned <- dat %>%
    dplyr::filter(cast_type != "one_depth") %>%
    dplyr::filter(shelf_loc_200 != "on_shelf" |
                    paste(Cruise, Transect_ID, cast_id) %in%
                    paste(casts_to_keep$Cruise, casts_to_keep$Transect_ID, casts_to_keep$cast_id))
  
  return(list(
    data = dat_thinned,
    cast_sequence = cast_sequence
  ))
}

if (opt$thinned) {
  dat_ISIIS_full <- thin_onshelf_casts(dat_ISIIS_full)$data
  dat_ISIIS_plankton <- thin_onshelf_casts(dat_ISIIS_plankton)$data
  }


assemblages <- list(Seston = dat_ISIIS_full, Plankton = dat_ISIIS_plankton)


# Supplementary diagnostics table and figures for every fitted model.
model_diagnostics_summary <- function(fit, dat, model_id, output_dir) {
  out <- diagnose_mixed_model_residuals(fit, dat, model_id, output_dir)
  short <- out$correlations %>%
    dplyr::filter(is.finite(correlation), lag_mid == 2.5) %>%
    dplyr::select(direction, correlation, n_pairs) %>%
    tidyr::pivot_wider(names_from = direction,
                       values_from = c(correlation, n_pairs), names_sep = '_')
  checks <- out$model_checks
  checks$model_id <- model_id
  checks$rmse <- sqrt(mean(stats::residuals(fit)^2))
  checks$max_abs_std_residual <- max(abs(stats::residuals(fit)/stats::sigma(fit)))
  checks$short_horizontal_correlation <- if ('correlation_Horizontal' %in% names(short)) short$correlation_Horizontal[1] else NA_real_
  checks$short_vertical_correlation <- if ('correlation_Vertical' %in% names(short)) short$correlation_Vertical[1] else NA_real_
  checks$n_short_horizontal_pairs <- if ('n_pairs_Horizontal' %in% names(short)) short$n_pairs_Horizontal[1] else NA_integer_
  checks$n_short_vertical_pairs <- if ('n_pairs_Vertical' %in% names(short)) short$n_pairs_Vertical[1] else NA_integer_
  checks
}
fixed_effect_table <- function(fit, model_id) {
  co <- as.data.frame(summary(fit)$coefficients)
  co$term <- rownames(co)
  ci <- suppressMessages(stats::confint(fit, parm = 'beta_', method = 'Wald'))
  co %>%
    dplyr::mutate(model_id = model_id,
                  ci_low = ci[term, 1], ci_high = ci[term, 2]) %>%
    dplyr::rename(estimate = Estimate, standard_error = `Std. Error`,
                  t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    dplyr::select(model_id, term, estimate, standard_error, t_value,
                  p_value, ci_low, ci_high)
}
model_r2 <- function(fit) {
  if (!requireNamespace('performance', quietly = TRUE)) return(c(NA_real_, NA_real_))
  x <- tryCatch(suppressWarnings(performance::r2_nakagawa(fit)), error = function(e) NULL)
  if (is.null(x)) return(c(NA_real_, NA_real_))
  c(as.numeric(x$R2_marginal), as.numeric(x$R2_conditional))
}


prep_mm_for_compare <- function(dat,
                                string_ID,
                                PSS_var = "slope",
                                depth_cut = 100,
                                log10_transform = FALSE,
                                random_effect = "(1 | Cruise/Transect_ID/cast_id)",
                                include_covariates = TRUE,
                                gap = 0.35) {
  
  # ------------------------------------------------------------
  # 1. Identify the saved model
  # ------------------------------------------------------------
  
  source <- dplyr::case_when(
    string_ID == "ISIIS_full" ~ "Seston",
    string_ID == "ISIIS_plankton" ~ "Plankton",
    string_ID %in% c("Seston", "Plankton") ~ string_ID,
    TRUE ~ NA_character_
  )
  
  if (is.na(source)) {
    stop("Unrecognized data source: ", string_ID)
  }
  
  model_id <- paste(source, PSS_var, sep = "_")
  
  model_path <- file.path(
    output_dir,
    paste0(model_id, ".rds")
  )
  
  if (!file.exists(model_path)) {
    stop("Saved model not found: ", model_path)
  }
  
  # ------------------------------------------------------------
  # 2. Load the original fitted model and its dataset
  # ------------------------------------------------------------
  
  saved <- readRDS(model_path)
  
  m1 <- saved$model
  df_mod <- saved$data
  
  if (!inherits(m1, "merMod")) {
    stop("Saved object is not a fitted mixed model: ", model_id)
  }
  
  if (nrow(df_mod) != stats::nobs(m1)) {
    stop("Saved data and model have different observation counts: ", model_id)
  }
  
  # The saved model uses 'y' as its response.
  # Total biovolume was already log10-transformed
  # before the original model was fitted.
  
  df <- df_mod %>%
    dplyr::mutate(
      .y = y,
      season = factor(
        season,
        levels = c("Winter", "Summer")
      ),
      shelf_loc_200 = factor(
        shelf_loc_200,
        levels = c("off_shelf", "on_shelf")
      )
    )
  
  # ------------------------------------------------------------
  # 3. Calculate EMMs from the SAME fitted model
  # ------------------------------------------------------------
  
  emm <- emmeans::emmeans(
    m1,
    ~ season * shelf_loc_200
  )
  
  pred_df <- as.data.frame(emm)
  
  if (all(c("asymp.LCL", "asymp.UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      dplyr::rename(
        lower.CL = asymp.LCL,
        upper.CL = asymp.UCL
      )
  } else if (all(c("LCL", "UCL") %in% names(pred_df))) {
    pred_df <- pred_df %>%
      dplyr::rename(
        lower.CL = LCL,
        upper.CL = UCL
      )
  }
  
  # ------------------------------------------------------------
  # 4. Prepare EMM positions for the comparison figure
  # ------------------------------------------------------------
  
  pred_df2 <- pred_df %>%
    dplyr::mutate(
      shelf_loc_200 = factor(
        shelf_loc_200,
        levels = c("off_shelf", "on_shelf")
      ),
      season = factor(
        season,
        levels = c("Winter", "Summer")
      ),
      x_base_comp = 1 +
        (as.numeric(shelf_loc_200) - 1) * gap
    )
  
  # ------------------------------------------------------------
  # 5. Return objects expected by plotting function
  # ------------------------------------------------------------
  
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
    dplyr::select(.data$shelf_loc_200, .data$data_source, .data$season, .data$x, .data$emmean) %>%
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
  
  if (PSS_var %in% c("total_biovolume", "total_biovolume_um3_m3", "integrated_biovol_all", "integrated_biovol_sample")) {
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
    scale_color_manual(
      values = source_colors,
      name = "Source",
      labels = c(
        allParticles = "Seston",
        plankton = "Plankton"
      )
    ) +
    labs(
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

make_allParticles_vs_plankton_totalBiovol_3x2 <- function(dat_all,
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
  
  row_total_biovol <- plot_compare_sources_row(
    dat_all = dat_all,
    dat_plankton = dat_plankton,
    PSS_var = "total_biovolume_um3_m3",
    depth_cut = 100,
    log10_transform = TRUE,
    y_limits = c(7, 19)
  )
  
  p_final <- (row_slope / row_intercept / row_total_biovol) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(p_final)
  
  out_path <- file.path(
    output_dir,
    "SPECTRA_allParticles_vs_plankton_violin_plusMixedModel_totalBiovolume_3x2.png"
  )
  
  ggsave(
    out_path,
    p_final,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  message("Saved total-biovolume comparison figure to: ", out_path)
  
  invisible(p_final)
}

output_dir <- file.path(base_dir, 'figures', 'size_spectra_mixed_models_cast_ID',
                        if (opt$thinned) 'season_shelf_LMM_thinned' else 'season_shelf_LMM')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_rows <- list()
coefficient_rows <- list()
check_rows <- list()
emm_rows <- list()
for (source in names(assemblages)) {
  dat <- assemblages[[source]]
  for (response in c('slope', 'intercept', 'total_biovolume_um3_m3')) {
    d <- dat %>%
      dplyr::filter(depth_mid < 100, !is.na(season), !is.na(shelf_loc_200),
                    !is.na(Cruise), !is.na(Transect_ID), !is.na(cast_id),
                    is.finite(depth_mid), is.finite(distance_sample_km)) %>%
      dplyr::group_by(Cruise, Transect_ID, group_id, season, shelf_loc_200) %>%
      dplyr::summarise(y_raw = dplyr::first(.data[[response]]),
                       depth_mid = dplyr::first(depth_mid),
                       distance_sample_km = dplyr::first(distance_sample_km),
                       cast_id = dplyr::first(cast_id),
                       .groups = 'drop') %>%
      dplyr::mutate(y = if (response == 'total_biovolume_um3_m3')
                      ifelse(y_raw > 0, log10(y_raw), NA_real_) else y_raw) %>%
      dplyr::filter(is.finite(y)) %>%
      droplevels()
    d <- as.data.frame(d)
    rownames(d) <- NULL
    id <- paste(source, response, sep = '_')
    if (nrow(d) < 20 || dplyr::n_distinct(d$Cruise) < 2 ||
        dplyr::n_distinct(d$season) < 2 ||
        dplyr::n_distinct(d$shelf_loc_200) < 2) {
      warning('Skipping ', id, ': insufficient observations or factor levels')
      next
    }
    fit <- lmerTest::lmer(
      y ~ season * shelf_loc_200 + scale(depth_mid) +
        scale(distance_sample_km) + (1 | Cruise/Transect_ID/cast_id),
      data = d, REML = TRUE
    )
    saveRDS(list(model = fit, data = d), file.path(output_dir, paste0(id, '.rds')))
    capture.output(summary(fit), anova(fit),
                   file = file.path(output_dir, paste0(id, '_summary.txt')))
    coefficient_rows[[id]] <- fixed_effect_table(fit, id)
    check_rows[[id]] <- model_diagnostics_summary(fit, d, id, output_dir)
    r2 <- model_r2(fit)
    model_rows[[id]] <- tibble::tibble(model_id = id, assemblage = source,
      response = response, n = nrow(d), n_casts = nrow(dplyr::distinct(d, Cruise, Transect_ID, cast_id)),
      marginal_r2 = r2[1], conditional_r2 = r2[2],
      singular_fit = lme4::isSingular(fit))
    emm <- as.data.frame(emmeans::emmeans(fit, ~ season * shelf_loc_200))
    emm$assemblage <- source
    emm$response <- response
    emm$model_id <- id
    emm_rows[[id]] <- emm
    p <- ggplot2::ggplot(d, ggplot2::aes(x = shelf_loc_200, y = y, fill = season)) +
      ggplot2::geom_violin(position = ggplot2::position_dodge(0.8), alpha = 0.5, trim = FALSE) +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), width = 0.15,
                            outlier.shape = NA, alpha = 0.5) +
      ggplot2::labs(title = id, x = 'Shelf location',
        y = if (response == 'total_biovolume_um3_m3') 'log10 total biovolume (µm³ m⁻³)' else response) +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(output_dir, paste0(id, '_observed.png')),
                    p, width = 8, height = 5, dpi = 300)
  }
}
readr::write_csv(dplyr::bind_rows(model_rows), file.path(output_dir, 'model_summary.csv'))
readr::write_csv(dplyr::bind_rows(coefficient_rows), file.path(output_dir, 'fixed_effects.csv'))
readr::write_csv(dplyr::bind_rows(check_rows), file.path(output_dir, 'model_diagnostics.csv'))
readr::write_csv(dplyr::bind_rows(emm_rows), file.path(output_dir, 'estimated_marginal_means.csv'))

# ------------------------------------------------------------
# Cruise comparisons: slope
# ------------------------------------------------------------

cruise_levels <- c("w22", "s22", "w23", "s23")

cruise_labels <- c(
  "w22" = "Winter 2022",
  "s22" = "Summer 2022",
  "w23" = "Winter 2023",
  "s23" = "Summer 2023"
)

cruise_colors <- c(
  "Winter" = "#0072B2",
  "Summer" = "#E69F00"
)

for (source in names(assemblages)) {
  
  dat <- assemblages[[source]]
  
  d_cruise <- dat %>%
    dplyr::filter(
      depth_mid < 100,
      !is.na(Cruise),
      !is.na(Transect_ID),
      !is.na(cast_id),
      is.finite(slope),
      is.finite(depth_mid),
      is.finite(distance_sample_km)
    ) %>%
    dplyr::group_by(Cruise, Transect_ID, group_id) %>%
    dplyr::summarise(
      slope = dplyr::first(slope),
      depth_mid = dplyr::first(depth_mid),
      distance_sample_km = dplyr::first(distance_sample_km),
      cast_id = dplyr::first(cast_id),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Cruise = factor(as.character(Cruise), levels = cruise_levels),
      cruise_season = factor(
        dplyr::if_else(
          as.character(Cruise) %in% c("w22", "w23"),
          "Winter",
          "Summer"
        ),
        levels = c("Winter", "Summer")
      ),
      Transect_ID = factor(Transect_ID),
      cast_id = factor(cast_id)
    ) %>%
    droplevels() %>%
    as.data.frame()
  
 
  
  # ----------------------------------------------------------
  # Mixed-effects ANCOVA
  # ----------------------------------------------------------
  
  fit_cruise <- lmerTest::lmer(
    slope ~ Cruise +
      scale(depth_mid) +
      scale(distance_sample_km) +
      (1 | Cruise:Transect_ID) +
      (1 | Cruise:Transect_ID:cast_id),
    data = d_cruise,
    REML = TRUE
  )
  
  message(
    "Singular cruise model (", source, "): ",
    lme4::isSingular(fit_cruise)
  )
  
  saveRDS(
    list(model = fit_cruise, data = d_cruise),
    file.path(output_dir, paste0(source, "_slope_cruise_LMM.rds"))
  )
  
  # Overall cruise test
  cruise_anova <- as.data.frame(
    stats::anova(fit_cruise, type = 3, ddf = "Satterthwaite")
  )
  
  cruise_anova$term <- rownames(cruise_anova)
  
  readr::write_csv(
    cruise_anova,
    file.path(output_dir, paste0(source, "_slope_cruise_ANOVA.csv"))
  )
  
  # ----------------------------------------------------------
  # Tukey-adjusted comparisons
  # ----------------------------------------------------------
  
  cruise_emm <- emmeans::emmeans(
    fit_cruise,
    ~ Cruise,
    lmer.df = "satterthwaite",
    lmerTest.limit = 30000,
    disable.pbkrtest = TRUE
  )
  
  cruise_tukey <- emmeans::contrast(
    cruise_emm,
    method = "pairwise",
    adjust = "tukey"
  )
  
  tukey_df <- as.data.frame(cruise_tukey)
  
  p_values <- stats::setNames(
    tukey_df$p.value,
    gsub(" - ", "-", tukey_df$contrast, fixed = TRUE)
  )
  
  cruise_letters <- multcompView::multcompLetters(
    p_values,
    threshold = 0.05
  )$Letters
  
  letters_df <- data.frame(
    Cruise = factor(names(cruise_letters), levels = cruise_levels),
    letter = unname(cruise_letters)
  )
  
  
  
  letter_y <- max(d_cruise$slope, na.rm = TRUE) +
    0.08 * diff(range(d_cruise$slope, na.rm = TRUE))
  
  letters_df$y <- letter_y
  
  readr::write_csv(
    as.data.frame(cruise_tukey),
    file.path(output_dir, paste0(source, "_slope_cruise_Tukey.csv"))
  )
  
  readr::write_csv(
    as.data.frame(cruise_emm),
    file.path(output_dir, paste0(source, "_slope_cruise_EMMs.csv"))
  )
  
  capture.output(
    summary(fit_cruise),
    stats::anova(fit_cruise, type = 3, ddf = "Satterthwaite"),
    cruise_emm,
    cruise_tukey,
    file = file.path(output_dir, paste0(source, "_slope_cruise_summary.txt"))
  )
  
  # ----------------------------------------------------------
  # Observed slope distributions with Tukey significance letters
  # ----------------------------------------------------------
  
  p_cruise <- ggplot2::ggplot(
    d_cruise,
    ggplot2::aes(
      x = Cruise,
      y = slope,
      fill = cruise_season
    )
  ) +
    ggplot2::geom_violin(
      trim = FALSE,
      alpha = 0.6,
      linewidth = 0.4
    ) +
    ggplot2::geom_boxplot(
      width = 0.13,
      fill = "white",
      outlier.shape = NA,
      alpha = 0.8
    ) +
    ggplot2::geom_text(
      data = letters_df,
      ggplot2::aes(
        x = Cruise,
        y = y,
        label = letter
      ),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    ) +
    ggplot2::scale_x_discrete(
      labels = cruise_labels
    ) +
    ggplot2::scale_fill_manual(
      values = cruise_colors,
      name = "Season"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.15))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "NBSS slope",
      title = paste(source, "NBSS slopes across cruises")
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 14),
      legend.position = "top"
    )
  
  ggplot2::ggsave(
    file.path(output_dir, paste0(source, "_slope_by_cruise_violin.png")),
    p_cruise,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  readr::write_csv(
  letters_df,
  file.path(output_dir, paste0(source, "_slope_cruise_letters.csv"))
)
  message("Completed cruise slope comparison for: ", source)
}

# ------------------------------------------------------------
# Combined 3x2 figure: slope, intercept, total biovolume
# ------------------------------------------------------------

make_allParticles_vs_plankton_totalBiovol_3x2(
  dat_all = dat_ISIIS_full,
  dat_plankton = dat_ISIIS_plankton,
  output_dir = output_dir,
  width = 13,
  height = 16,
  dpi = 300
)

message('Completed season × shelf models: ', output_dir)
