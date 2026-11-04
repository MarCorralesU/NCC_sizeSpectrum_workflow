library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)
library(NLP)
library(here)
library(lubridate)
library(runner)
library(suncalc)
library(purrr)
library(slider)
library(rlang)


collapse_taxa_cols <- function(
  df,
  start_name = "appendicularian",
  end_name   = "tunicate_salp",
  ignore_case = FALSE,
  detritus_handling = c("keep", "remove")  # <-- NEW
) {
  detritus_handling <- match.arg(detritus_handling)

  nm <- names(df)
  if (ignore_case) {
    nm_l <- tolower(nm)
    s <- match(tolower(start_name), nm_l)
    e <- match(tolower(end_name),   nm_l)
  } else {
    s <- match(start_name, nm)
    e <- match(end_name,   nm)
  }
  if (is.na(s) || is.na(e) || s > e) {
    stop(
      "Couldn't locate a valid taxa block between '",
      start_name, "' and '", end_name, "'."
    )
  }

  # Columns to collapse (between anchors)
  data_cols <- nm[s:e]

  # Define category patterns
  patterns <- list(
    detritus             = "detritus|art[ie]fact",  # allows both
    protist              = "protist",
    diatom               = "diatom",
    copepod              = "copepod",
    crustacean           = "crustacean",
    predatory_gelatinous = "(hydromedusae|chaetognath|siphonophore)",
    grazer_gelatinous    = "(appendicularian|tunicate|pyrosome)",
    fish                 = "fish"
  )

  # Logical match vectors
  idx <- lapply(patterns, function(p) grepl(p, data_cols, ignore.case = TRUE))
  
  # Prevent diatoms from also being counted as phytoplankton, not needed, all phyto in this case are diatoms

  # Taxa buckets
  detritus_cols <- data_cols[idx$detritus]

  matched_non_det <- Reduce(
    `|`,
    idx[c(
      "protist","diatom","copepod","crustacean",  #"phytoplankton",
      "predatory_gelatinous","grazer_gelatinous","fish"
    )]
  )

  matched_cols <- data_cols[matched_non_det]

  # Any remaining go to "other", NEVER detritus
  other_cols <- setdiff(data_cols, c(matched_cols, detritus_cols))

  # Safe row sums
  rs <- function(cols) {
    if (!length(cols)) return(rep(0, nrow(df)))
    m <- lapply(df[, cols, drop = FALSE],
                function(z) suppressWarnings(as.numeric(z)))
    rowSums(as.data.frame(m), na.rm = TRUE)
  }

  # Build grouped columns
  grouped <- tibble::tibble(
    protist              = rs(data_cols[idx$protist]),
    diatom               = rs(data_cols[idx$diatom]),
    copepod              = rs(data_cols[idx$copepod]),
    crustacean           = rs(data_cols[idx$crustacean]),
    predatory_gelatinous = rs(data_cols[idx$predatory_gelatinous]),
    grazer_gelatinous    = rs(data_cols[idx$grazer_gelatinous]),
    fish                 = rs(data_cols[idx$fish]),
    other                = rs(other_cols)
  )

  # LEFT & RIGHT dataframe segments
  keep_left  <- if (s > 1)        df[, 1:(s-1), drop = FALSE]        else df[, 0, drop = FALSE]
  keep_right <- if (e < ncol(df)) df[, (e+1):ncol(df), drop = FALSE] else df[, 0, drop = FALSE]

  # detritus output option
  if (detritus_handling == "keep" && length(detritus_cols)) {
    keep_det <- df[, detritus_cols, drop = FALSE]
  } else {
    keep_det <- NULL
  }

  dplyr::bind_cols(keep_left, grouped, keep_det, keep_right)

}


diagnose_mixed_model_residuals <- function(model, df_mod, model_prefix, output_dir,
                                           horizontal_breaks = seq(0, 100, by = 5),
                                           vertical_breaks = seq(0, 100, by = 5),
                                           max_neighbors = 20,
                                           min_pairs = 10) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ------------------------------------------------------------
  # 1. Extract observations actually used by the fitted model
  # ------------------------------------------------------------
  
  model_rows <- suppressWarnings(as.integer(rownames(stats::model.frame(model))))
  
  if (anyNA(model_rows) || length(model_rows) != stats::nobs(model)) {
    stop("Could not align model observations with df_mod.")
  }
  
  diag_dat <- df_mod[model_rows, , drop = FALSE] %>%
    dplyr::mutate(
      residual = as.numeric(stats::residuals(model)),
      fitted_value = as.numeric(stats::fitted(model)),
      standardized_residual = residual / stats::sigma(model)
    ) %>%
    dplyr::filter(
      is.finite(standardized_residual),
      is.finite(fitted_value)
    )
  
  if (!"cast_id" %in% names(diag_dat)) {
    stop("cast_id is required for vertical residual diagnostics.")
  }
  
  # ------------------------------------------------------------
  # 2. Calculate neighboring residual pairs without constructing
  #    every possible pair of observations
  # ------------------------------------------------------------
  
  make_lag_pairs <- function(dat, coordinate, breaks, direction) {
    
    dat <- dat[is.finite(dat[[coordinate]]), , drop = FALSE]
    dat <- dat[order(dat[[coordinate]]), , drop = FALSE]
    
    n <- nrow(dat)
    
    if (n < 2) {
      return(tibble::tibble())
    }
    
    pair_list <- vector("list", min(max_neighbors, n - 1))
    
    for (k in seq_len(min(max_neighbors, n - 1))) {
      
      i <- seq_len(n - k)
      j <- i + k
      
      pair_list[[k]] <- tibble::tibble(
        lag = abs(dat[[coordinate]][j] - dat[[coordinate]][i]),
        residual_1 = dat$standardized_residual[i],
        residual_2 = dat$standardized_residual[j]
      )
    }
    
    dplyr::bind_rows(pair_list) %>%
      dplyr::filter(
        is.finite(lag),
        lag > 0,
        lag <= max(breaks)
      ) %>%
      dplyr::mutate(
        direction = direction,
        lag_bin = cut(
          lag,
          breaks = breaks,
          include.lowest = TRUE,
          right = FALSE),
        lag_mid = ((breaks[-1] + breaks[-length(breaks)]) / 2)[as.integer(lag_bin)]
      )%>%
      dplyr::filter(!is.na(lag_bin))
  }
  
  # ------------------------------------------------------------
  # 3. Horizontal autocorrelation
  #    Compare observations at the same depth within transects
  # ------------------------------------------------------------
  
  horizontal_pairs <- diag_dat %>%
    dplyr::filter(
      is.finite(distance_sample_km),
      is.finite(depth_mid)
    ) %>%
    dplyr::group_by(Cruise, Transect_ID, depth_mid) %>%
    dplyr::group_modify(
      ~ make_lag_pairs(
        .x,
        coordinate = "distance_sample_km",
        breaks = horizontal_breaks,
        direction = "Horizontal"
      )
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------
  # 4. Vertical autocorrelation
  #    Compare depths within individual casts
  # ------------------------------------------------------------
  
  vertical_pairs <- diag_dat %>%
    dplyr::filter(
      !is.na(cast_id),
      is.finite(depth_mid)
    ) %>%
    dplyr::group_by(Cruise, Transect_ID, cast_id) %>%
    dplyr::group_modify(
      ~ make_lag_pairs(
        .x,
        coordinate = "depth_mid",
        breaks = vertical_breaks,
        direction = "Vertical"
      )
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------
  # 5. Summarize residual correlations by separation distance
  # ------------------------------------------------------------
  
  correlation_data <- dplyr::bind_rows(
    horizontal_pairs,
    vertical_pairs
  ) %>%
    dplyr::group_by(direction, lag_bin, lag_mid) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      correlation = if (
        dplyr::n() >= min_pairs &&
        stats::sd(residual_1) > 0 &&
        stats::sd(residual_2) > 0
      ) {
        stats::cor(residual_1, residual_2)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  
  readr::write_csv(
    correlation_data,
    file.path(output_dir, paste0(model_prefix, "_residual_correlations.csv"))
  )
  
  # Summarize short-distance residual autocorrelation
  
  short_lag_summary <- correlation_data %>%
    dplyr::filter(
      (direction == "Horizontal" & lag_mid == 2.5) |
        (direction == "Vertical" & lag_mid == 7.5)
    ) %>%
    dplyr::select(
      direction,
      lag_mid,
      n_pairs,
      correlation
    )
  
  cat("\nSHORT-DISTANCE RESIDUAL AUTOCORRELATION\n")
  cat("Model:", model_prefix, "\n")
  print(short_lag_summary)
  cat("========================================\n\n")
  
  # ------------------------------------------------------------
  # 6. Panel A: Horizontal residual correlation
  # ------------------------------------------------------------
  
  p_horizontal <- correlation_data %>%
    dplyr::filter(direction == "Horizontal") %>%
    ggplot2::ggplot(
      ggplot2::aes(x = lag_mid, y = correlation)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_point(color = "#0072B2", size = 2) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "A. Horizontal residual correlation",
      x = "Along-transect separation (km)",
      y = "Pearson correlation"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 7. Panel B: Vertical residual correlation
  # ------------------------------------------------------------
  
  p_vertical <- correlation_data %>%
    dplyr::filter(direction == "Vertical") %>%
    ggplot2::ggplot(
      ggplot2::aes(x = lag_mid, y = correlation)
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_line(color = "#D55E00", linewidth = 0.8) +
    ggplot2::geom_point(color = "#D55E00", size = 2) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = "B. Vertical residual correlation",
      x = "Depth separation (m)",
      y = "Pearson correlation"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 8. Panel C: Standardized residuals versus fitted values
  # ------------------------------------------------------------
  
  p_fitted <- ggplot2::ggplot(
    diag_dat,
    ggplot2::aes(x = fitted_value, y = standardized_residual)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_point(
      alpha = 0.2, size = 0.7, color = "#0072B2"
    ) +
    ggplot2::geom_smooth(
      method = "loess", se = FALSE, color = "#D55E00",
      linewidth = 0.9
    ) +
    ggplot2::labs(
      title = "C. Residuals vs. fitted values",
      x = "Fitted values",
      y = "Standardized residuals"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 9. Panel D: Normal Q-Q plot
  # ------------------------------------------------------------
  
  p_qq <- ggplot2::ggplot(
    diag_dat,
    ggplot2::aes(sample = standardized_residual)
  ) +
    ggplot2::stat_qq(
      alpha = 0.35, size = 0.8, color = "#0072B2"
    ) +
    ggplot2::stat_qq_line(
      color = "#D55E00", linewidth = 0.9
    ) +
    ggplot2::labs(
      title = "D. Normal Q-Q plot",
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    ) +
    ggplot2::theme_classic(base_size = 12)
  
  # ------------------------------------------------------------
  # 10. Assemble and save the four-panel diagnostic figure
  # ------------------------------------------------------------
  
  p_diagnostics <- (p_horizontal | p_vertical) /
    (p_fitted | p_qq)
  
  p_diagnostics <- p_diagnostics +
    patchwork::plot_annotation(
      title = paste("Mixed-effects model diagnostics:", model_prefix)
    )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(model_prefix, "_residual_diagnostics_4panel.png")
    ),
    plot = p_diagnostics,
    width = 13,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  # ------------------------------------------------------------
  # 11. Export model residuals and additional diagnostics
  # ------------------------------------------------------------
  
  readr::write_csv(
    diag_dat,
    file.path(output_dir, paste0(model_prefix, "_model_residuals.csv"))
  )
  
  convergence_messages <- model@optinfo$conv$lme4$messages
  
  model_checks <- tibble::tibble(
    n_observations = stats::nobs(model),
    residual_sd = stats::sigma(model),
    singular_fit = lme4::isSingular(model),
    convergence_messages = if (length(convergence_messages) == 0) {
      "None"
    } else {
      paste(convergence_messages, collapse = "; ")
    }
  )
  
  readr::write_csv(
    model_checks,
    file.path(output_dir, paste0(model_prefix, "_model_checks.csv"))
  )
  
  message("Saved residual diagnostics for: ", model_prefix)
  
  invisible(list(
    residual_data = diag_dat,
    correlations = correlation_data,
    model_checks = model_checks,
    plot = p_diagnostics
  ))
}


overwrite_distance_from_shore <- function(
    dat,
    land_ref_path = "ancillary/transect_shelf_land.csv"
) {

  land_ref <- readr::read_csv(
    land_ref_path,
    show_col_types = FALSE
  ) %>%
    dplyr::select(Transect_ID, land_long)

  if (!all(c("Transect_ID", "Longitude", "Latitude") %in% names(dat))) {
    stop("Missing Transect_ID, Longitude, or Latitude.")
  }

  dat <- dat %>%
    dplyr::select(-dplyr::any_of("land_long")) %>%
    dplyr::left_join(
      land_ref,
      by = "Transect_ID",
      relationship = "many-to-one"
    )

  if (any(is.na(dat$land_long))) {
    stop("Some samples have no matching coastal reference longitude.")
  }

  dat <- dat %>%
    dplyr::mutate(
      distance_sample_km = geosphere::distGeo(
        cbind(Longitude, Latitude),
        cbind(land_long, Latitude)
      ) / 1000
    ) %>%
    dplyr::select(-land_long)

  return(dat)
}
