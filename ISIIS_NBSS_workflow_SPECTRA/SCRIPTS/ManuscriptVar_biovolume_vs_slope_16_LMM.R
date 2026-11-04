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
})
options <- list(
  optparse::make_option('--thinned', action = 'store_true', default = FALSE,
                        help = 'Thin on-shelf casts to every third cast')
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = options))
setwd('/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
base_dir <- getwd()
source(file.path(base_dir, "SCRIPTS", "functions", "model_functions.R"))
sizeSpectraDirPlankton <- file.path(base_dir, 'data', 'size_spectra_plankton_5-10m-bins_thresholded')
sizeSpectraDirAll <- file.path(base_dir, 'data', 'ISIIS-LISST_mergedNBSS_5m_allParticles')
df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = 'by_community\\.csv$', full.names = TRUE, recursive = TRUE)
df_list_full <- list.files(sizeSpectraDirAll, pattern = 'LISST-ISIIS_NBSS\\.csv$', full.names = TRUE, recursive = TRUE)
env_vars <- c('Temperature', 'Salinity', 'Oxygen_ml_l', 'pH', 'PAR', 'chl_a_ug_l')
if (!length(df_list_full) || !length(df_list_plankton)) stop('No input files found. Run from SPECTRA project root.')
load_ISIIS_data <- function(df_list,
                            env_vars = c("Temperature","Salinity","Oxygen_ml_l","pH","PAR","chl_a_ug_l"),
                            x1_full = 12.78849,
                            x2_full = 4.547442e13) {
  
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
    dplyr::filter(is.finite(correlation), n_pairs >= 10) %>%
    dplyr::group_by(direction) %>%
    dplyr::slice_min(lag_mid, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(direction = factor(direction, levels = c("Horizontal", "Vertical"))) %>%
    dplyr::select(direction, correlation, n_pairs) %>%
    tidyr::pivot_wider(names_from = direction,
                       values_from = c(correlation, n_pairs),
                       names_sep = '_',
                       names_expand = TRUE)
  
  if (!nrow(short)) {
    short <- tibble::tibble(
      correlation_Horizontal = NA_real_,
      correlation_Vertical = NA_real_,
      n_pairs_Horizontal = NA_integer_,
      n_pairs_Vertical = NA_integer_
    )
  }
  
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

output_dir <- file.path(base_dir, 'figures', 'size_spectra_mixed_models_cast_ID',
                        if (opt$thinned) 'biovolume_vs_slope_16LMM_thinned' else 'biovolume_vs_slope_16LMM')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_rows <- list()
coefficient_rows <- list()
check_rows <- list()
plot_data <- list()
for (source in names(assemblages)) {
  dat <- assemblages[[source]]
  d0 <- dat %>%
    dplyr::filter(depth_mid < 100, !is.na(Cruise), !is.na(Transect_ID),
                  !is.na(cast_id), !is.na(shelf_loc_200)) %>%
    dplyr::group_by(Cruise, Transect_ID, group_id, shelf_loc_200) %>%
    dplyr::summarise(slope = dplyr::first(slope),
                     biovolume = dplyr::first(total_biovolume_um3_m3),
                     depth_mid = dplyr::first(depth_mid),
                     distance_sample_km = dplyr::first(distance_sample_km),
                     cast_id = dplyr::first(cast_id), .groups = 'drop') %>%
    dplyr::mutate(log10_biovolume = ifelse(biovolume > 0, log10(biovolume), NA_real_)) %>%
    dplyr::filter(is.finite(log10_biovolume), is.finite(slope),
                  is.finite(depth_mid), is.finite(distance_sample_km))
  for (cruise in c('w22', 's22', 'w23', 's23')) {
    for (shelf in c('off_shelf', 'on_shelf')) {
      d <- d0 %>% dplyr::filter(as.character(Cruise) == cruise,
                                as.character(shelf_loc_200) == shelf) %>% droplevels()
      d <- as.data.frame(d)
      rownames(d) <- NULL
      id <- paste(source, cruise, shelf, sep = '_')
      n_casts <- nrow(dplyr::distinct(d, Transect_ID, cast_id))
      n_transects <- dplyr::n_distinct(d$Transect_ID)
      if (nrow(d) < 10 || n_transects < 2 || n_casts < 3 ||
          nrow(d) <= n_casts + 3 || dplyr::n_distinct(d$slope) < 3 ||
          dplyr::n_distinct(d$depth_mid) < 2 ||
          dplyr::n_distinct(d$distance_sample_km) < 2) {
        warning('Skipping ', id, ': insufficient replication or variation')
        model_rows[[id]] <- tibble::tibble(model_id = id, assemblage = source,
          Cruise = cruise, shelf = shelf, n = nrow(d), n_casts = n_casts,
          status = 'insufficient replication', marginal_r2 = NA_real_,
          conditional_r2 = NA_real_)
        next
      }
      fit <- tryCatch(lmerTest::lmer(
        log10_biovolume ~ slope + scale(depth_mid) + scale(distance_sample_km) +
          (1 | Transect_ID/cast_id), data = d, REML = TRUE
      ), error = function(e) e)
      if (inherits(fit, 'error')) {
        warning(id, ': ', conditionMessage(fit))
        model_rows[[id]] <- tibble::tibble(model_id = id, assemblage = source,
          Cruise = cruise, shelf = shelf, n = nrow(d), n_casts = n_casts,
          status = conditionMessage(fit), marginal_r2 = NA_real_,
          conditional_r2 = NA_real_)
        next
      }
      saveRDS(list(model = fit, data = d), file.path(output_dir, paste0(id, '.rds')))
      capture.output(summary(fit), file = file.path(output_dir, paste0(id, '_summary.txt')))
      coefficient_rows[[id]] <- fixed_effect_table(fit, id)
      check_rows[[id]] <- model_diagnostics_summary(fit, d, id, output_dir)
      r2 <- model_r2(fit)
      baseline_fit <- lmerTest::lmer(
        log10_biovolume ~ scale(depth_mid) + scale(distance_sample_km) +
          (1 | Transect_ID/cast_id),
        data = d,
        REML = TRUE
      )
      
      baseline_r2 <- model_r2(baseline_fit)
      delta_marginal_r2 <- r2[1] - baseline_r2[1]
      model_rows[[id]] <- tibble::tibble(model_id = id, assemblage = source,
        Cruise = cruise, shelf = shelf, n = nrow(d), n_casts = n_casts,
        status = if (lme4::isSingular(fit)) 'singular' else 'OK',
        marginal_r2 = r2[1],
        conditional_r2 = r2[2],
        baseline_marginal_r2 = baseline_r2[1],
        delta_marginal_r2 = delta_marginal_r2)
      # Population-level fit at mean depth/distance, random effects = zero.
      grid <- data.frame(slope = seq(min(d$slope), max(d$slope), length.out = 100),
                         depth_mid = mean(d$depth_mid),
                         distance_sample_km = mean(d$distance_sample_km),
                         Transect_ID = d$Transect_ID[1], cast_id = d$cast_id[1])
      grid$fit <- stats::predict(fit, newdata = grid, re.form = NA)
      d$source <- source
      d$Cruise <- cruise
      d$shelf_loc_200 <- shelf
      grid$source <- source
      grid$Cruise <- cruise
      grid$shelf_loc_200 <- shelf
      plot_data[[id]] <- list(observed = d, fitted = grid)
    }
  }
}
readr::write_csv(dplyr::bind_rows(model_rows), file.path(output_dir, 'model_summary_16.csv'))
readr::write_csv(dplyr::bind_rows(coefficient_rows), file.path(output_dir, 'fixed_effects_16.csv'))
readr::write_csv(dplyr::bind_rows(check_rows), file.path(output_dir, 'model_diagnostics_16.csv'))
if (length(plot_data)) {
  
  observed <- dplyr::bind_rows(lapply(plot_data, `[[`, "observed"))
  fitted <- dplyr::bind_rows(lapply(plot_data, `[[`, "fitted"))
  
  cruise_labels <- c(
    w22 = "Winter 2022",
    s22 = "Summer 2022",
    w23 = "Winter 2023",
    s23 = "Summer 2023"
  )
  
  shelf_colors <- c(
    off_shelf = "#0072B2",
    on_shelf = "#D55E00"
  )
  
  observed <- observed %>%
    dplyr::mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = names(cruise_labels),
                      labels = cruise_labels),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("off_shelf", "on_shelf"))
    )
  
  fitted <- fitted %>%
    dplyr::mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = names(cruise_labels),
                      labels = cruise_labels),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("off_shelf", "on_shelf")),
      biovolume_fit = 10^fit
    )
  
  annotations <- dplyr::bind_rows(model_rows) %>%
    dplyr::filter(is.finite(delta_marginal_r2)) %>%
    dplyr::mutate(
      source = factor(assemblage, levels = c("Seston", "Plankton")),
      Cruise = factor(Cruise, levels = names(cruise_labels),
                      labels = cruise_labels),
      shelf = factor(shelf, levels = c("off_shelf", "on_shelf")),
      label = sprintf("\u0394R\u00b2m = %.3f", delta_marginal_r2),
      annotation_row = ifelse(shelf == "off_shelf", 1.4, 2.8)
    )
  
  panel_labels <- tidyr::expand_grid(
    source = factor(c("Seston", "Plankton"), levels = c("Seston", "Plankton")),
    Cruise = factor(unname(cruise_labels), levels = unname(cruise_labels))
  ) %>%
    dplyr::mutate(
      panel_letter = paste0(letters[1:8], ".")
    )
  
  p <- ggplot2::ggplot(
    observed,
    ggplot2::aes(
      x = slope,
      y = biovolume,
      color = shelf_loc_200
    )
  ) +
    ggplot2::geom_point(alpha = 0.18, size = 0.7) +
    ggplot2::geom_line(
      data = fitted,
      ggplot2::aes(
        x = slope,
        y = biovolume_fit,
        color = shelf_loc_200
      ),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    ggplot2::geom_text(
      data = annotations,
      ggplot2::aes(
        x = Inf,
        y = Inf,
        label = label,
        color = shelf
      ),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = annotations$annotation_row,
      size = 5.3,
      show.legend = FALSE
    ) +
    ggplot2::facet_grid(
      source ~ Cruise,
      scales = "free_x",
      axes = "all_x",
      axis.labels = "all_x"
    ) +
    ggplot2::geom_text(
      data = panel_labels,
      ggplot2::aes(x = -Inf, y = Inf, label = panel_letter),
      inherit.aes = FALSE,
      hjust = -0.2,
      vjust = 1.5,
      size = 5.5,
      fontface = "bold"
    ) +
    ggplot2::scale_y_log10(
      breaks = 10^seq(8, 20, by = 4),
      labels = scales::label_math(10^.x, format = log10)
    ) +
    ggplot2::scale_color_manual(
      values = shelf_colors,
      breaks = c("off_shelf", "on_shelf"),
      labels = c("Off shelf", "On shelf"),
      name = "Shelf location"
    ) +
    ggplot2::labs(
      x = "NBSS slope",
      y = expression("Total biovolume"~(mu*m^3~m^-3))
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 17),
      legend.text = ggplot2::element_text(size = 15),
      legend.key.width = grid::unit(1.5, "cm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 13, face = "bold"),
      axis.title = ggplot2::element_text(size = 15),
      panel.spacing = grid::unit(0.6, "lines")
    )
  
  ggplot2::ggsave(
    file.path(output_dir, "SPECTRA_biovolume_vs_slope_2x4.png"),
    p,
    width = 16,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}
message('Completed biovolume vs slope models: ', output_dir)
