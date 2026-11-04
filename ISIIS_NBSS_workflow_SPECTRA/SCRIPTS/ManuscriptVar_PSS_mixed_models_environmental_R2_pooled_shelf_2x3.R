#!/usr/bin/env Rscript
# SPECTRA: six environmental LMMs of NBSS slope (3 predictors x 2 assemblages).
# One model per predictor for seston and plankton, pooled across cruises/shelf regions.
# All three predictors are compared on identical complete-case samples per assemblage.

setwd('/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(lubridate)
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(patchwork)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option('--thinned', action = 'store_true', default = FALSE,
                        help = 'Thin on-shelf casts (every third cast) before fitting')
)))
base_dir <- getwd()

source(file.path(base_dir, 'SCRIPTS', 'functions', 'model_functions.R'))
sizeSpectraDirPlankton <- file.path(base_dir, 'data', 'size_spectra_plankton_5-10m-bins_thresholded')
sizeSpectraDirAll <- file.path(base_dir, 'data', 'ISIIS-LISST_mergedNBSS_5m_allParticles')
output_dir_figs_NB <- file.path(base_dir, 'figures', 'size_spectra_mixed_models_cast_ID')
analysis_output_dir <- if (opt$thinned) file.path(output_dir_figs_NB, 'revisions') else output_dir_figs_NB

# Load only columns needed for deduplication and environmental regressions.
load_slope_data <- function(directory, pattern) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop('No matching input files in: ', directory)
  purrr::map2_dfr(files, seq_along(files), function(path, source_id) {
    readr::read_csv(path, show_col_types = FALSE, col_select = c(
      group_id, time_Pacific, rel_dist, depth_mid, Latitude, Longitude, slope, intercept,
      Cruise, Transect_ID, cast_id, distance_sample_km,
      Temperature, chl_a_ug_l, Oxygen_ml_l,
      logNB, biovol_um3_width
    )) %>%
      dplyr::mutate(source_id = source_id)
  }) %>%
    dplyr::mutate(
      time_Pacific = lubridate::ymd_hms(time_Pacific, quiet = TRUE),
      group_id_old = group_id,
      biovolume_um3_m3 = (10^logNB) * biovol_um3_width
    ) %>%
    dplyr::mutate(
      group_id = dplyr::cur_group_id(),
      .by = c(source_id, group_id_old)
    )  %>%
    dplyr::group_by(group_id) %>%
    dplyr::mutate(total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE)) %>%
    dplyr::ungroup()
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



dat_ISIIS_full <- remove_duplicate_NBSS_samples(
  load_slope_data(sizeSpectraDirAll, 'LISST-ISIIS_NBSS\\.csv$')
)$data

dat_ISIIS_plankton <- remove_duplicate_NBSS_samples(
  load_slope_data(sizeSpectraDirPlankton, 'by_community\\.csv$')
)$data

dat_ISIIS_full <- overwrite_distance_from_shore(dat_ISIIS_full)
dat_ISIIS_plankton <- overwrite_distance_from_shore(dat_ISIIS_plankton)

clean_slope_data <- function(dat) {
  dat %>%
    dplyr::mutate(
      shelf_loc_200 = dplyr::case_when(
        rel_dist >= 0 ~ 'on_shelf',
        rel_dist < 0 ~ 'off_shelf',
        TRUE ~ NA_character_
      ),
      shelf_loc_200 = factor(shelf_loc_200, levels = c('on_shelf', 'off_shelf')),
      Cruise = factor(Cruise, levels = c('w22', 's22', 'w23', 's23')),
      Transect_ID = factor(Transect_ID, levels = c('GH', 'OB', 'CR', 'CM', 'NH', 'HH', 'RR'))
    ) %>%
    dplyr::filter(!(Cruise == 'w22' & Transect_ID == 'CM' & Longitude < -124.75))
}

dat_ISIIS_full <- clean_slope_data(dat_ISIIS_full)
dat_ISIIS_plankton <- clean_slope_data(dat_ISIIS_plankton)

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
  add_cast_type <- function(dat) {
    types <- dat %>%
      dplyr::filter(!is.na(Cruise), !is.na(Transect_ID), !is.na(cast_id),
                    !is.na(time_Pacific), !is.na(depth_mid)) %>%
      dplyr::distinct(Cruise, Transect_ID, cast_id, time_Pacific, depth_mid) %>%
      dplyr::arrange(Cruise, Transect_ID, cast_id, time_Pacific) %>%
      dplyr::group_by(Cruise, Transect_ID, cast_id) %>%
      dplyr::summarise(
        n_depths = dplyr::n_distinct(depth_mid),
        depth_change = dplyr::last(depth_mid) - dplyr::first(depth_mid),
        cast_type = dplyr::case_when(
          n_depths == 1 ~ 'one_depth',
          depth_change > 0 ~ 'down_cast',
          depth_change < 0 ~ 'up_cast',
          TRUE ~ NA_character_
        ),
        .groups = 'drop'
      )
    dplyr::left_join(dat, dplyr::select(types, Cruise, Transect_ID, cast_id, cast_type),
                     by = c('Cruise', 'Transect_ID', 'cast_id'))
  }
  dat_ISIIS_full <- thin_onshelf_casts(add_cast_type(dat_ISIIS_full))$data
  dat_ISIIS_plankton <- thin_onshelf_casts(add_cast_type(dat_ISIIS_plankton))$data
}

dat_ISIIS_full_analysis <- dat_ISIIS_full
dat_ISIIS_plankton_analysis <- dat_ISIIS_plankton

# ============================================================
# NEW: Environmental predictors of NBSS slope
#
# Three independent LMMs:
#   slope ~ z(Temperature) + z(depth) + z(distance) + nested random intercepts
#   slope ~ z(log10(chl-a)) + z(depth) + z(distance) + nested random intercepts
#   slope ~ z(Oxygen) + z(depth) + z(distance) + nested random intercepts
#
# Models pooled across cruises to estimate:
#   (1 | Cruise/Transect_ID/cast_id)
# Separate model comparisons for Seston and Plankton, pooling both shelf regions.
# All three models in a comparison use identical complete-case samples.
# Shelf location is included as a fixed-effect adjustment (no interaction).
#
# ============================================================

if (!requireNamespace("performance", quietly = TRUE)) {
  stop("Package 'performance' is required for marginal/conditional R2. Install with install.packages('performance').")
}

environmental_output_dir <- file.path(analysis_output_dir, "envVars_slope_LMM")
dir.create(environmental_output_dir, recursive = TRUE, showWarnings = FALSE)

environmental_vars <- c("Temperature", "log10_chl_a_ug_l", "Oxygen_ml_l")
environmental_labels <- c(
  Temperature = "Temperature",
  log10_chl_a_ug_l = "log10 chlorophyll-a",
  Oxygen_ml_l = "Dissolved oxygen"
)

# Keep one observation per NBSS group_id, not one observation per size class.
prepare_environmental_samples <- function(dat, source_name, depth_cut = 100) {
  required <- c(
    "group_id", "slope", "Temperature", "chl_a_ug_l", "Oxygen_ml_l",
    "depth_mid", "distance_sample_km", "Cruise", "Transect_ID",
    "cast_id", "shelf_loc_200"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop("Missing columns for ", source_name, ": ", paste(missing, collapse = ", "))
  }

  dat %>%
    dplyr::filter(
      is.finite(depth_mid), depth_mid < depth_cut,
      !is.na(Cruise), !is.na(Transect_ID),
      !is.na(cast_id), !is.na(shelf_loc_200)
    ) %>%
    dplyr::group_by(Cruise, Transect_ID, group_id, shelf_loc_200) %>%
    dplyr::summarise(
      slope = dplyr::first(slope),
      Temperature = dplyr::first(Temperature),
      chl_a_ug_l = dplyr::first(chl_a_ug_l),
      Oxygen_ml_l = dplyr::first(Oxygen_ml_l),
      depth_mid = dplyr::first(depth_mid),
      distance_sample_km = dplyr::first(distance_sample_km),
      cast_id = dplyr::first(cast_id),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      source = source_name,
      log10_chl_a_ug_l = dplyr::if_else(
        is.finite(chl_a_ug_l) & chl_a_ug_l > 0,
        log10(chl_a_ug_l), NA_real_
      )
    ) %>%
    dplyr::filter(
      is.finite(slope),
      is.finite(Temperature),
      is.finite(log10_chl_a_ug_l),
      is.finite(Oxygen_ml_l),
      is.finite(depth_mid),
      is.finite(distance_sample_km)
    ) %>%
    dplyr::mutate(
      Cruise = droplevels(factor(Cruise)),
      Transect_ID = droplevels(factor(Transect_ID)),
      cast_id = factor(cast_id),
      shelf_loc_200 = factor(shelf_loc_200, levels = c("off_shelf", "on_shelf"))
    )
}

env_samples <- dplyr::bind_rows(
  prepare_environmental_samples(dat_ISIIS_full_analysis, "Seston"),
  prepare_environmental_samples(dat_ISIIS_plankton_analysis, "Plankton")
)

readr::write_csv(
  env_samples,
  file.path(environmental_output_dir, "SPECTRA_environmental_common_complete_cases.csv")
)

standardize_column <- function(x, name) {
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) {
    stop("Cannot standardize ", name, ": zero or undefined standard deviation.")
  }
  as.numeric((x - mean(x)) / s)
}


fit_environmental_comparison <- function(d, source_i) {
  
  z_predictors <- c(
    Temperature = "z_temperature",
    log10_chl_a_ug_l = "z_chla",
    Oxygen_ml_l = "z_oxygen"
  )
  
  fit_one <- function(predictor) {
    
    # Exclude Summer 2023 from oxygen models only.
    d_model <- d
    
    if (predictor == "Oxygen_ml_l") {
      d_model <- d_model %>%
        dplyr::filter(as.character(Cruise) != "s23")
    }
    
    d_model <- droplevels(d_model)
    
    n_shelf_regions <- dplyr::n_distinct(d_model$shelf_loc_200)
    n_cruises <- dplyr::n_distinct(d_model$Cruise)
    n_transects <- dplyr::n_distinct(
      interaction(d_model$Cruise, d_model$Transect_ID, drop = TRUE)
    )
    n_casts <- dplyr::n_distinct(
      interaction(d_model$Cruise, d_model$Transect_ID,
                  d_model$cast_id, drop = TRUE)
    )
    
    if (nrow(d_model) < 20 || n_shelf_regions < 2 ||
        n_cruises < 2 || n_transects < 3 ||
        n_casts < 5 || nrow(d_model) <= n_casts + 4) {
      warning("Skipping ", source_i, " / ", predictor,
              ": insufficient replication.")
      return(NULL)
    }
    
    # Standardize using each model's actual dataset.
    d_model <- d_model %>%
      dplyr::mutate(
        z_depth = standardize_column(depth_mid, "depth_mid"),
        z_distance = standardize_column(distance_sample_km, "distance_sample_km"),
        z_temperature = standardize_column(Temperature, "Temperature"),
        z_chla = standardize_column(log10_chl_a_ug_l, "log10_chl_a_ug_l"),
        z_oxygen = standardize_column(Oxygen_ml_l, "Oxygen_ml_l")
      )
    
    z_name <- unname(z_predictors[[predictor]])
    
    baseline_formula <- slope ~ shelf_loc_200 + z_depth + z_distance +
      (1 | Cruise/Transect_ID/cast_id)
    
    baseline <- lmerTest::lmer(
      baseline_formula,
      data = d_model,
      REML = TRUE
    )
    
    baseline_r2 <- suppressWarnings(
      performance::r2_nakagawa(baseline)
    )
    
    baseline_marginal <- as.numeric(baseline_r2$R2_marginal)
    
    formula <- stats::as.formula(paste(
      "slope ~", z_name,
      "+ shelf_loc_200 + z_depth + z_distance",
      "+ (1 | Cruise/Transect_ID/cast_id)"
    ))
    
    model <- lmerTest::lmer(
      formula,
      data = d_model,
      REML = TRUE
    )
    
    coefs <- summary(model)$coefficients
    
    ci <- suppressMessages(
      stats::confint(model, parm = z_name, method = "Wald")
    )
    
    r2 <- suppressWarnings(
      performance::r2_nakagawa(model)
    )
    
    prefix <- paste(
      "SPECTRA", source_i, "pooled_shelf", predictor,
      sep = "_"
    )
    
    # ------------------------------------------------------------
    # Residual diagnostics: horizontal and vertical autocorrelation
    # ------------------------------------------------------------
    
    diagnose_mixed_model_residuals(
      model = model,
      df_mod = d_model,
      model_prefix = paste(source_i, predictor, sep = "_"),
      output_dir = environmental_output_dir
    )
    
    saveRDS(
      list(
        model = model,
        data = d_model,
        formula = formula,
        baseline = baseline
      ),
      file.path(environmental_output_dir,
                paste0(prefix, "_model.rds"))
    )
    
    capture.output(
      summary(model),
      file = file.path(environmental_output_dir,
                       paste0(prefix, "_summary.txt"))
    )
    
    tibble::tibble(
      source = source_i,
      shelf_scope = "pooled_on_and_off_shelf",
      predictor = predictor,
      excluded_cruises = if (predictor == "Oxygen_ml_l") "s23" else "None",
      n = nrow(d_model),
      n_shelf_regions = n_shelf_regions,
      n_cruises = n_cruises,
      n_transects = n_transects,
      n_casts = n_casts,
      beta_standardized_predictor = unname(lme4::fixef(model)[z_name]),
      standard_error = coefs[z_name, "Std. Error"],
      p_value = coefs[z_name, "Pr(>|t|)"],
      ci_low = unname(ci[1]),
      ci_high = unname(ci[2]),
      baseline_marginal_r2 = baseline_marginal,
      marginal_r2 = as.numeric(r2$R2_marginal),
      conditional_r2 = as.numeric(r2$R2_conditional),
      delta_marginal_r2 = as.numeric(r2$R2_marginal) - baseline_marginal,
      singular_fit = lme4::isSingular(model),
      convergence_messages = paste(
        model@optinfo$conv$lme4$messages,
        collapse = "; "
      )
    )
  }
  
  results <- purrr::map_dfr(names(z_predictors), fit_one)
  
  results <- results %>%
    dplyr::arrange(dplyr::desc(marginal_r2)) %>%
    dplyr::mutate(rank_by_marginal_r2 = dplyr::row_number())
  
  results
}


environmental_results <- env_samples %>%
  dplyr::group_by(source) %>%
  dplyr::group_split(.keep = TRUE) %>%
  purrr::map_dfr(function(d) {
    fit_environmental_comparison(
      d,
      source_i = as.character(d$source[1])
    )
  })

readr::write_csv(
  environmental_results,
  file.path(environmental_output_dir, "SPECTRA_environmental_model_R2_comparison.csv")
)

if (nrow(environmental_results) > 0) {
  p_r2 <- environmental_results %>%
    dplyr::mutate(
      predictor = factor(
        predictor,
        levels = c("Temperature", "log10_chl_a_ug_l", "Oxygen_ml_l"),
        labels = c("Temperature", "log10 chlorophyll-a", "Dissolved oxygen")
      )
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = predictor, y = marginal_r2)) +
    ggplot2::geom_col(width = 0.7, fill = "#0072B2") +
    ggplot2::facet_wrap(~ source, nrow = 1) +
    ggplot2::labs(
      x = NULL,
      y = expression("Marginal " * R^2),
      title = "Environmental predictors of NBSS slope",
      subtitle = "Both shelf regions pooled; adjusted for shelf, standardized depth and distance"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  ggplot2::ggsave(
    file.path(environmental_output_dir, "SPECTRA_environmental_marginal_R2_comparison.png"),
    p_r2, width = 11, height = 7, dpi = 300, bg = "white"
  )

  # Build explicitly so the factor labeling is preserved.
  p_delta <- environmental_results %>%
    dplyr::mutate(
      predictor = factor(
        predictor,
        levels = c("Temperature", "log10_chl_a_ug_l", "Oxygen_ml_l"),
        labels = c("Temperature", "log10 chlorophyll-a", "Dissolved oxygen")
      )
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = predictor, y = delta_marginal_r2)) +
    ggplot2::geom_col(width = 0.7, fill = "#D55E00") +
    ggplot2::facet_wrap(~ source, nrow = 1) +
    ggplot2::labs(
      x = NULL,
      y = expression(Delta * " marginal " * R^2),
      title = "Additional explanatory value beyond shelf location, depth and distance"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

  ggplot2::ggsave(
    file.path(environmental_output_dir, "SPECTRA_environmental_delta_marginal_R2.png"),
    p_delta, width = 11, height = 7, dpi = 300, bg = "white"
  )
}

message("Environmental R2 comparisons saved in: ", environmental_output_dir)


# ============================================================
# NEW: Six-panel environmental LMM visualization and two tables
# ============================================================
# Each panel displays observed NBSS slopes and a population-level
# (random effects excluded) model estimate, averaged equally over
# the two shelf regions with depth and distance held at their means.
# Shading is a 95% CI for the fixed-effect mean, not a prediction interval.

make_environmental_lmm_outputs <- function(results, output_dir) {
  if (nrow(results) != 6) {
    warning("Expected six environmental models; found ", nrow(results),
            ". Available models will be plotted and summarized.")
  }

  predictors <- c("Temperature", "log10_chl_a_ug_l", "Oxygen_ml_l")
  z_variables <- c(
    Temperature = "z_temperature",
    log10_chl_a_ug_l = "z_chla",
    Oxygen_ml_l = "z_oxygen"
  )
  predictor_labels <- c(
    Temperature = "Temperature (°C)",
    log10_chl_a_ug_l = "log10 chlorophyll-a (µg L⁻¹)",
    Oxygen_ml_l = "Dissolved oxygen (mL L⁻¹)"
  )

  plot_data <- list()
  plot_fits <- list()
  diagnostics <- list()

  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    source_i <- as.character(r$source)
    predictor_i <- as.character(r$predictor)
    prefix <- paste("SPECTRA", source_i, "pooled_shelf", predictor_i, sep = "_")
    object <- readRDS(file.path(output_dir, paste0(prefix, "_model.rds")))
    model <- object$model
    d <- object$data
    z_name <- unname(z_variables[[predictor_i]])

    # Raw points on original environmental-variable scale.
    plot_data[[i]] <- tibble::tibble(
      source = source_i,
      predictor = predictor_i,
      x = d[[predictor_i]],
      slope = d$slope,
      shelf_loc_200 = d$shelf_loc_200
    )

    # New observations at the means of depth and distance; use both
    # shelf regions and average their fixed-effect predictions equally.
    x_seq <- seq(min(d[[predictor_i]]), max(d[[predictor_i]]), length.out = 150)
    z_seq <- (x_seq - mean(d[[predictor_i]])) / stats::sd(d[[predictor_i]])
    newdata <- expand.grid(
      x_index = seq_along(x_seq),
      shelf_loc_200 = levels(d$shelf_loc_200),
      KEEP.OUT.ATTRS = FALSE
    )
    newdata$shelf_loc_200 <- factor(
      newdata$shelf_loc_200,
      levels = levels(d$shelf_loc_200)
    )
    newdata$z_depth <- 0
    newdata$z_distance <- 0
    newdata$z_temperature <- 0
    newdata$z_chla <- 0
    newdata$z_oxygen <- 0
    newdata[[z_name]] <- z_seq[newdata$x_index]

    fixed_formula <- lme4::nobars(stats::formula(model))
    design <- stats::model.matrix(
      stats::delete.response(stats::terms(fixed_formula)),
      data = newdata
    )
    beta <- lme4::fixef(model)
    vc <- as.matrix(stats::vcov(model))
    design <- design[, names(beta), drop = FALSE]

    averaged_design <- vapply(
      seq_along(x_seq),
      function(j) colMeans(design[newdata$x_index == j, , drop = FALSE]),
      numeric(length(beta))
    )
    averaged_design <- t(averaged_design)
    colnames(averaged_design) <- names(beta)
    fitted_mean <- as.numeric(averaged_design %*% beta)
    fitted_se <- sqrt(rowSums((averaged_design %*% vc) * averaged_design))

    plot_fits[[i]] <- tibble::tibble(
      source = source_i,
      predictor = predictor_i,
      x = x_seq,
      fit = fitted_mean,
      lower = fitted_mean - 1.96 * fitted_se,
      upper = fitted_mean + 1.96 * fitted_se
    )

    resid <- as.numeric(stats::residuals(model))
    fitted <- as.numeric(stats::fitted(model))
    std_resid <- resid / stats::sigma(model)
    conv <- model@optinfo$conv$lme4$messages
    correlations <- c(
      residual_fitted = suppressWarnings(stats::cor(resid, fitted)),
      abs_residual_fitted = suppressWarnings(stats::cor(abs(resid), fitted))
    )
    diagnostics[[i]] <- tibble::tibble(
      source = source_i,
      predictor = predictor_i,
      n = stats::nobs(model),
      residual_sd = stats::sigma(model),
      residual_mean = mean(resid),
      residual_rmse = sqrt(mean(resid^2)),
      standardized_residual_mean = mean(std_resid),
      standardized_residual_sd = stats::sd(std_resid),
      standardized_residual_min = min(std_resid),
      standardized_residual_max = max(std_resid),
      n_abs_standardized_residual_gt_3 = sum(abs(std_resid) > 3),
      residual_fitted_correlation = unname(correlations["residual_fitted"]),
      abs_residual_fitted_correlation = unname(correlations["abs_residual_fitted"]),
      singular_fit = lme4::isSingular(model),
      convergence_messages = if (length(conv)) paste(conv, collapse = "; ") else "None"
    )
  }

  # Table 1: all six environmental model results, including marginal and
  # conditional R2 and incremental marginal R2 against the common baseline.
  model_results_table <- results %>%
    dplyr::mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      predictor = factor(predictor, levels = predictors)
    ) %>%
    dplyr::arrange(source, predictor)
  readr::write_csv(
    model_results_table,
    file.path(output_dir, "SPECTRA_environmental_LMM_model_results_table.csv")
  )

  # Table 2: numerical model diagnostics for the same six models.
  model_diagnostics_table <- dplyr::bind_rows(diagnostics) %>%
    dplyr::mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      predictor = factor(predictor, levels = predictors)
    ) %>%
    dplyr::arrange(source, predictor)
  readr::write_csv(
    model_diagnostics_table,
    file.path(output_dir, "SPECTRA_environmental_LMM_diagnostics_table.csv")
  )

  observations <- dplyr::bind_rows(plot_data) %>%
    dplyr::mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      predictor = factor(predictor, levels = predictors, labels = predictor_labels)
    )
  estimates <- dplyr::bind_rows(plot_fits) %>%
    dplyr::mutate(
      source = factor(source, levels = c("Seston", "Plankton")),
      predictor = factor(predictor, levels = predictors, labels = predictor_labels)
    )

  
  # Panel labels in reading order.
  panel_labels <- tidyr::expand_grid(
    source = factor(c("Seston", "Plankton"),
                    levels = c("Seston", "Plankton")),
    predictor = factor(
      predictor_labels,
      levels = predictor_labels
    )
  ) %>%
    dplyr::mutate(
      panel = paste0(letters[1:6], ".")
    )
  
  
  
  p_2x3 <- ggplot2::ggplot(
    observations,
    ggplot2::aes(x = x, y = slope)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = shelf_loc_200),
      alpha = 0.20,
      size = 0.65
    ) +
    ggplot2::geom_line(
      data = estimates,
      ggplot2::aes(x = x, y = fit),
      inherit.aes = FALSE,
      color = "#111111",
      linewidth = 1
    ) +
    ggplot2::geom_text(
      data = panel_labels,
      ggplot2::aes(x = -Inf, y = Inf, label = panel),
      inherit.aes = FALSE,
      hjust = -0.15,
      vjust = -0.5,
      fontface = "bold",
      size = 5
    ) +
    ggplot2::facet_grid(
      source ~ predictor,
      scales = "free_x",
      switch = "x"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        off_shelf = "#0072B2",
        on_shelf = "#D55E00"
      ),
      labels = c(
        off_shelf = "Off shelf",
        on_shelf = "On shelf"
      ),
      name = "Shelf region"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = expression(NBSS~Slope~(m^-3~mu*m^-3))
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text.x.bottom = ggplot2::element_text(
        face = "bold",
        size = 12,
        margin = ggplot2::margin(t = 8)
      ),
      strip.text.y = ggplot2::element_text(
        face = "bold",
        size = 11
      ),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.margin = ggplot2::margin(15, 20, 15, 15)
    )
  

  ggplot2::ggsave(
    file.path(output_dir, "SPECTRA_environmental_LMM_slope_2x3.png"),
    p_2x3, width = 15, height = 9, dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, "SPECTRA_environmental_LMM_slope_2x3.pdf"),
    p_2x3, width = 15, height = 9, bg = "white"
  )

  invisible(list(
    plot = p_2x3,
    model_results = model_results_table,
    diagnostics = model_diagnostics_table
  ))
}

environmental_lmm_outputs <- make_environmental_lmm_outputs(
  results = environmental_results,
  output_dir = environmental_output_dir
)
