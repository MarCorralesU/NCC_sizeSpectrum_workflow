#!/usr/bin/env Rscript
# Run after ManuscriptVar_PSS_mixed_models_only_sensitivity_cast_ID.R.
# Reads existing residual-correlation and model-check CSV outputs; does not refit models.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path(getwd(), "figures", "size_spectra_mixed_models_cast_ID", "season_shelf_LMM")
output_dir <- if (length(args) >= 2) args[[2]] else file.path(input_dir, "supplementary_diagnostics")
model_family <- basename(normalizePath(input_dir))

response_levels <- switch(
  model_family,
  season_shelf_LMM = c("Slope", "Intercept", "Log10 total biovolume"),
  envVars_slope_LMM = c("Temperature", "Chlorophyll-a", "Dissolved oxygen"),
  biovolume_vs_slope_16LMM = c(
    "w22 off shelf", "w22 on shelf",
    "s22 off shelf", "s22 on shelf",
    "w23 off shelf", "w23 on shelf",
    "s23 off shelf", "s23 on shelf"
  ),
  biovolume_vs_slope_16LMM_thinned = c(
    "w22 off shelf", "w22 on shelf",
    "s22 off shelf", "s22 on shelf",
    "w23 off shelf", "w23 on shelf",
    "s23 off shelf", "s23 on shelf"
  )
)

if (is.null(response_levels)) {
  stop("Unsupported model family: ", model_family)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cor_files <- list.files(input_dir, pattern = "_residual_correlations\\.csv$", full.names = TRUE)
if (!length(cor_files)) stop("No *_residual_correlations.csv files in: ", input_dir)


parse_model <- function(filename) {
  
  stem <- sub(
    "_residual_correlations\\.csv$",
    "",
    basename(filename)
  )
  
  pattern <- switch(
    model_family,
    season_shelf_LMM = "^(Seston|Plankton)_(slope|intercept|total_biovolume_um3_m3)$",
    envVars_slope_LMM = "^(Seston|Plankton)_(Temperature|log10_chl_a_ug_l|Oxygen_ml_l)$",
    biovolume_vs_slope_16LMM = "^(Seston|Plankton)_(w22|s22|w23|s23)_(off_shelf|on_shelf)$",
    biovolume_vs_slope_16LMM_thinned = "^(Seston|Plankton)_(w22|s22|w23|s23)_(off_shelf|on_shelf)$"
  )
  
  parts <- stringr::str_match(stem, pattern)
  
  tibble::tibble(
    model_id = stem,
    dataset = parts[, 2],
    response = dplyr::case_when(
      parts[, 3] == "slope" ~ "Slope",
      parts[, 3] == "intercept" ~ "Intercept",
      parts[, 3] == "total_biovolume_um3_m3" ~ "Log10 total biovolume",
      parts[, 3] == "Temperature" ~ "Temperature",
      parts[, 3] == "log10_chl_a_ug_l" ~ "Chlorophyll-a",
      parts[, 3] == "Oxygen_ml_l" ~ "Dissolved oxygen",
      model_family %in% c(
        "biovolume_vs_slope_16LMM",
        "biovolume_vs_slope_16LMM_thinned"
      ) & !is.na(parts[, 3]) ~ paste(
        parts[, 3],
        gsub("_", " ", parts[, ncol(parts)])
      ),
      TRUE ~ NA_character_
    )
  )
}





model_index <- map_dfr(cor_files, parse_model) %>%
  filter(!is.na(dataset), !is.na(response)) %>%
  distinct(model_id, .keep_all = TRUE)

if (!nrow(model_index)) {
  stop(
    "No matching primary model names found in: ",
    input_dir,
    "\nFiles detected:\n",
    paste(basename(cor_files), collapse = "\n")
  )
}



correlations <- map_dfr(seq_len(nrow(model_index)), function(i) {
  meta <- model_index[i, ]
  file <- file.path(input_dir, paste0(meta$model_id, "_residual_correlations.csv"))
  read_csv(file, show_col_types = FALSE) %>%
    mutate(model_id = meta$model_id, dataset = meta$dataset, response = meta$response)
}) %>%
  mutate(
    dataset = factor(dataset, levels = c("Seston", "Plankton")),
    response = factor(response, levels = response_levels),
    direction = factor(direction, levels = c("Horizontal", "Vertical"))
  )

checks <- map_dfr(seq_len(nrow(model_index)), function(i) {
  meta <- model_index[i, ]
  file <- file.path(input_dir, paste0(meta$model_id, "_model_checks.csv"))
  if (!file.exists(file)) return(tibble())
  read_csv(file, show_col_types = FALSE) %>%
    mutate(model_id = meta$model_id, dataset = meta$dataset, response = meta$response)
})

# Select the shortest separation bin with a finite correlation and at least 10 pairs.
# Horizontal separation is measured in km; vertical separation is measured in m.
# Correlations are descriptive: overlapping pairs are not independent.

short_lag <- correlations %>%
  dplyr::filter(
    is.finite(correlation),
    n_pairs >= 10
  ) %>%
  dplyr::group_by(model_id, direction) %>%
  dplyr::slice_min(lag_mid, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(direction = factor(direction, levels = c("Horizontal", "Vertical"))) %>%
  dplyr::select(model_id, direction, lag_mid, n_pairs, correlation) %>%
  tidyr::pivot_wider(
    id_cols = model_id,
    names_from = direction,
    values_from = c(lag_mid, n_pairs, correlation),
    names_glue = "{direction}_{.value}",
    names_expand = TRUE
  )


diagnostics_file <- file.path(
  input_dir,
  if (model_family %in% c(
    "biovolume_vs_slope_16LMM",
    "biovolume_vs_slope_16LMM_thinned"
  )) "model_diagnostics_16.csv" else "model_diagnostics.csv"
)

if (file.exists(diagnostics_file)) {
  
  additional_diagnostics <- readr::read_csv(
    diagnostics_file,
    show_col_types = FALSE
  ) %>%
    dplyr::select(model_id, rmse, max_abs_std_residual)
  
} else {
  
  additional_diagnostics <- purrr::map_dfr(model_index$model_id, function(id) {
    
    file <- file.path(input_dir, paste0(id, "_model_residuals.csv"))
    dat <- readr::read_csv(file, show_col_types = FALSE)
    
    tibble::tibble(
      model_id = id,
      rmse = sqrt(mean(dat$residual^2, na.rm = TRUE)),
      max_abs_std_residual = max(abs(dat$standardized_residual), na.rm = TRUE)
    )
  })
  
}



summary_table <- model_index %>%
  dplyr::left_join(short_lag, by = "model_id") %>%
  dplyr::left_join(checks %>%dplyr::select(model_id,n_observations,residual_sd,singular_fit,convergence_messages), by = "model_id") %>%
  dplyr::left_join(additional_diagnostics, by = "model_id") %>%
  dplyr::mutate(
    response = factor(response, levels = response_levels),
    dataset = factor(dataset, levels = c("Seston", "Plankton"))
  ) %>%
  dplyr::arrange(response, dataset) %>%
  dplyr::mutate(
    horizontal_status = dplyr::case_when(
      is.na(Horizontal_correlation) ~ "Not estimable",
      abs(Horizontal_correlation) >= 0.3 ~ "Inspect residual pattern",
      TRUE ~ "Low short-lag correlation"
    ),
    vertical_status = dplyr::case_when(
      is.na(Vertical_correlation) ~ "Not estimable",
      abs(Vertical_correlation) >= 0.3 ~ "Inspect residual pattern",
      TRUE ~ "Low short-lag correlation"
    )
  ) %>%
  dplyr::select(
    dataset,
    response,
    n_observations,
    residual_sd,
    rmse,
    max_abs_std_residual,
    singular_fit,
    convergence_messages,
    Horizontal_lag_mid,
    Horizontal_n_pairs,
    Horizontal_correlation,
    horizontal_status,
    Vertical_lag_mid,
    Vertical_n_pairs,
    Vertical_correlation,
    vertical_status
  )



write_csv(summary_table, file.path(output_dir, paste0(model_family, "_model_diagnostics_summary.csv")))
write_csv(correlations, file.path(output_dir, paste0(model_family, "_all_residual_correlations.csv")))

is_biovol <- model_family %in% c(
  "biovolume_vs_slope_16LMM",
  "biovolume_vs_slope_16LMM_thinned"
)

p <- ggplot(correlations, aes(x = lag_mid, y = correlation, color = direction, group = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.65, na.rm = TRUE) +
  geom_point(size = 1.25, na.rm = TRUE) +
  {
    if (is_biovol) {
      facet_wrap(~ response + dataset, nrow = 4, ncol = 4)
    } else {
      facet_grid(response ~ dataset)
    }
  } +
  scale_color_manual(values = c("Horizontal" = "#0072B2", "Vertical" = "#D55E00")) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "Separation distance (km horizontal; m vertical)",
       y = "Pearson correlation of conditional residual pairs",
       color = "Direction") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold"), panel.spacing = grid::unit(0.8, "lines"))

figure_width <- if (is_biovol) 16 else 12
figure_height <- if (is_biovol) 14 else 10

ggsave(file.path(output_dir, paste0(model_family, "_Figure_Sx_cast_model_residual_correlations.png")),
       p, width = figure_width, height = figure_height, dpi = 300, bg = "white")
ggsave(file.path(output_dir, paste0(model_family, "_Figure_Sx_cast_model_residual_correlations.pdf")),
       p, width = figure_width, height = figure_height, bg = "white")
print(summary_table, width = Inf)
message("Saved supplementary outputs to: ", normalizePath(output_dir))
