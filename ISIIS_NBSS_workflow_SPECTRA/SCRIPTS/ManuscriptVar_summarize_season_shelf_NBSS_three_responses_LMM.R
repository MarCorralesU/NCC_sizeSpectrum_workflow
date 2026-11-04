#!/usr/bin/env Rscript
# Summarize the six existing SPECTRA season x shelf REML-LMMs.
# Usage: Rscript SCRIPTS/summarize_season_shelf_NBSS_three_responses_LMM.R [input_dir]
# Run from the SPECTRA project root; no models are refitted.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[[1]] else file.path('figures', 'size_spectra_mixed_models_cast_ID', 'season_shelf_LMM')
if (!dir.exists(input_dir)) stop('Input directory not found: ', input_dir)
output_dir <- file.path(input_dir, 'summary_tables')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fractions <- c('Seston', 'Plankton')
responses <- c(slope = 'Slope', intercept = 'Intercept', total_biovolume_um3_m3 = 'Total Biovolume')
seasons <- c('Winter', 'Summer')
shelves <- c('off_shelf', 'on_shelf')

# Four significant figures; compact scientific notation for |x| < 0.001 or > 1000.
format_four_sig <- function(x) {
  vapply(x, function(v) {
    if (is.na(v) || !is.finite(v)) return('NA')
    if (v == 0) return('0')
    if (abs(v) < 0.001 || abs(v) > 1000) {
      return(sub('e([+-])0+([0-9]+)$', 'e\\1\\2', formatC(v, format = 'e', digits = 3)))
    }
    format(signif(v, 4), scientific = FALSE, trim = TRUE)
  }, character(1))
}

# Return original-scale response values: slope unchanged, intercept and biovolume back-transformed.
to_original_scale <- function(x, response) {
  if (response %in% c('intercept', 'total_biovolume_um3_m3')) 10^x else x
}

safe_stat <- function(x, fun) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  fun(x)
}

predictor_label <- function(term) {
  if (grepl('^season[^:]*:shelf_loc_200|^shelf_loc_200[^:]*:season', term)) return('Season x Shelf')
  if (grepl('^season', term)) return('Season')
  if (grepl('^shelf_loc_200', term)) return('Shelf')
  if (grepl('^scale\\(depth_mid\\)', term)) return('Depth')
  if (grepl('^scale\\(distance_sample_km\\)', term)) return('Distance from Shore')
  NA_character_
}

extract_tables <- function(fraction, response) {
  id <- paste(fraction, response, sep = '_')
  path <- file.path(input_dir, paste0(id, '.rds'))
  if (!file.exists(path)) stop('Missing expected model: ', path)
  saved <- readRDS(path)
  fit <- saved$model
  dat <- saved$data
  if (!inherits(fit, 'merMod') || is.null(dat)) stop('Expected list(model, data) in ', path)
  if (nrow(dat) != stats::nobs(fit)) stop('Observation-count mismatch in ', path)
  needed <- c('season', 'shelf_loc_200', 'y_raw', 'y')
  if (!all(needed %in% names(dat))) stop('Missing model-data columns in ', path, ': ', paste(setdiff(needed, names(dat)), collapse = ', '))

  # Table 1: summaries of observations actually used by each fitted model.
  raw <- dat %>%
    dplyr::mutate(Season = as.character(season), Shelf_location = as.character(shelf_loc_200),
                  original_value = to_original_scale(y, response)) %>%
    dplyr::group_by(Season, Shelf_location) %>%
    dplyr::summarise(n = dplyr::n(),
      Median = safe_stat(original_value, median),
      SD = safe_stat(original_value, sd),
      Q1 = safe_stat(original_value, function(z) unname(quantile(z, 0.25))),
      Q3 = safe_stat(original_value, function(z) unname(quantile(z, 0.75))),
      .groups = 'drop') %>%
    dplyr::mutate(Variable = unname(responses[[response]]), Fraction = fraction,
      SD = if (response == 'slope') SD else NA_real_,
      Q1 = if (response == 'slope') NA_real_ else Q1,
      Q3 = if (response == 'slope') NA_real_ else Q3) %>%
    dplyr::select(Variable, Fraction, Season, Shelf_location, n, Median, SD, Q1, Q3)

  # Table 2: coefficient estimates (excluding the reference-category intercept).
  co <- as.data.frame(coef(summary(fit)))
  required <- c('Estimate', 'Std. Error', 'df', 't value', 'Pr(>|t|)')
  if (!all(required %in% names(co))) stop('Missing lmerTest coefficient columns in ', path)
  co$term <- rownames(co)
  fixed <- tibble::tibble(Variable = unname(responses[[response]]), Fraction = fraction,
    Model_predictors = vapply(co$term, predictor_label, character(1)),
    Estimate = co[['Estimate']], Std_error = co[['Std. Error']],
    df = co[['df']], t_value = co[['t value']], p = co[['Pr(>|t|)']]) %>%
    dplyr::filter(!is.na(Model_predictors))
  if (nrow(fixed) != 5) warning('Expected five non-intercept fixed effects in ', id, '; found ', nrow(fixed))

  # Table 3: variance components, including residual variation.
  vc <- as.data.frame(lme4::VarCorr(fit))
  random <- tibble::tibble(Variable = unname(responses[[response]]), Fraction = fraction,
    Random_effects = ifelse(vc$grp == 'Residual', 'Residual', paste0(vc$grp, ': (Intercept)')),
    Variance = vc$vcov, SD = vc$sdcor)

  # Table 4: population-level EMMs at reference-grid means of continuous covariates.
  # emmeans handles scale(depth_mid) and scale(distance_sample_km) in the original formula.
  em <- emmeans::emmeans(fit, ~ season * shelf_loc_200, data = dat, lmer.df = "asymptotic")
  em_df <- as.data.frame(summary(em, infer = c(TRUE, FALSE), level = 0.95))
  lo <- intersect(c('lower.CL', 'asymp.LCL', 'LCL'), names(em_df))
  hi <- intersect(c('upper.CL', 'asymp.UCL', 'UCL'), names(em_df))
  if (!length(lo) || !length(hi)) stop('Cannot identify 95% CI columns in emmeans for ', id)
  emm <- tibble::tibble(Variable = unname(responses[[response]]),
                        Season = as.character(em_df$season), Fraction = fraction,
                        Shelf_location = as.character(em_df$shelf_loc_200),
                        EMM = to_original_scale(em_df$emmean, response),
                        Lower_CI = to_original_scale(em_df[[lo[1]]], response),
                        Upper_CI = to_original_scale(em_df[[hi[1]]], response))

  list(descriptive = raw, fixed = fixed, random = random, emm = emm)
}

model_keys <- tidyr::expand_grid(Fraction = fractions, response = names(responses))
results <- purrr::map2(model_keys$Fraction, model_keys$response, extract_tables)

order_table <- function(dat) {
  dat %>%
    dplyr::mutate(Variable = factor(Variable, levels = unname(responses)),
      Fraction = factor(Fraction, levels = fractions),
      dplyr::across(dplyr::any_of('Season'), ~ factor(.x, levels = seasons)),
      dplyr::across(dplyr::any_of('Shelf_location'), ~ factor(.x, levels = shelves))) %>%
    dplyr::arrange(Variable, Fraction, dplyr::across(dplyr::any_of(c('Season', 'Shelf_location')))) %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character))
}

write_formatted <- function(dat, filename) {
  dat <- order_table(dat) %>%
    dplyr::mutate(dplyr::across(where(is.numeric) & !dplyr::any_of('n'), format_four_sig))
  path <- file.path(output_dir, filename)
  readr::write_csv(dat, path, na = 'NA')
  message('Saved: ', path)
}

write_formatted(purrr::map_dfr(results, 'descriptive'), '01_NBSS_descriptive_statistics_4sig.csv')
write_formatted(purrr::map_dfr(results, 'fixed'), '02_NBSS_fixed_effects_4sig.csv')
write_formatted(purrr::map_dfr(results, 'random'), '03_NBSS_random_effects_4sig.csv')
write_formatted(purrr::map_dfr(results, 'emm'), '04_NBSS_estimated_marginal_means_4sig.csv')
message('Completed summaries for ', nrow(model_keys), ' saved models (no refitting).')
