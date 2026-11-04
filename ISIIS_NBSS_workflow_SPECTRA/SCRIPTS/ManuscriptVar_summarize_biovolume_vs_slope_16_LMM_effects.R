#!/usr/bin/env Rscript
# Summarize existing SPECTRA biovolume-vs-slope LMMs; only baseline models are fitted.
# Run from the SPECTRA project root. Optional argument: input directory.
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(lme4)
  library(lmerTest)
})
if (!requireNamespace('performance', quietly = TRUE)) stop("Install package 'performance'.")
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[[1]] else file.path('.', 'figures', 'size_spectra_mixed_models_cast_ID', 'biovolume_vs_slope_16LMM')
if (!dir.exists(input_dir)) stop('Input directory does not exist: ', input_dir)
output_dir <- file.path(input_dir, 'summary_tables')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_files <- list.files(input_dir, pattern = '^(Seston|Plankton)_(w22|s22|w23|s23)_(on_shelf|off_shelf)\\.rds$', full.names = TRUE)
if (!length(model_files)) stop('No matching model .rds files in ', input_dir)
if (length(model_files) != 16) warning('Expected 16 models; found ', length(model_files))

format_four_sig <- function(x) {
  vapply(x, function(value) {
    if (is.na(value) || !is.finite(value)) return(NA_character_)
    if (value == 0) return('0')
    if (abs(value) < 0.001 || abs(value) > 1000) {
      return(sub('e([+-])0+([0-9]+)$', 'e\\1\\2', formatC(value, format = 'e', digits = 3)))
    }
    format(signif(value, 4), scientific = FALSE, trim = TRUE)
  }, character(1))
}

safe_r2 <- function(fit) {
  out <- tryCatch(suppressWarnings(performance::r2_nakagawa(fit)), error = function(e) NULL)
  if (is.null(out)) return(c(marginal = NA_real_, conditional = NA_real_))
  c(marginal = as.numeric(out$R2_marginal), conditional = as.numeric(out$R2_conditional))
}

extract_model <- function(path) {
  id <- sub('\\.rds$', '', basename(path))
  parts <- stringr::str_match(id, '^(Seston|Plankton)_(w22|s22|w23|s23)_(on_shelf|off_shelf)$')
  fraction <- parts[1, 2]
  cruise <- parts[1, 3]
  shelf <- parts[1, 4]
  obj <- readRDS(path)
  fit <- if (inherits(obj, 'merMod')) obj else obj$model
  if (!inherits(fit, 'merMod')) stop('No mixed model found in ', path)
  if (is.null(obj$data)) stop('Saved model data missing in ', path, '; needed for comparable baseline.')

  coefs <- as.data.frame(coef(summary(fit)))
  required <- c('Estimate', 'Std. Error', 'df', 't value', 'Pr(>|t|)')
  if (!all(required %in% names(coefs))) stop('Missing lmerTest coefficient columns in ', path)
  coefs$term <- rownames(coefs)

  full_r2 <- safe_r2(fit)
  baseline_formula <- stats::update(stats::formula(fit), . ~ . - slope)
  baseline_fit <- tryCatch(
    lmerTest::lmer(baseline_formula, data = obj$data, REML = TRUE),
    error = function(e) { warning('Baseline model failed for ', id, ': ', conditionMessage(e)); NULL }
  )
  baseline_marginal <- if (is.null(baseline_fit)) NA_real_ else unname(safe_r2(baseline_fit)['marginal'])
  delta_marginal <- unname(full_r2['marginal']) - baseline_marginal

  fixed <- tibble::tibble(
    Fraction = fraction,
    Cruise = cruise,
    Shelf_region = shelf,
    Fixed_effect_predictor = coefs$term,
    Estimate = coefs[['Estimate']],
    Std_error = coefs[['Std. Error']],
    df = coefs[['df']],
    t_value = coefs[['t value']],
    p = coefs[['Pr(>|t|)']],
    Baseline_marginal_R2 = baseline_marginal,
    Marginal_R2 = unname(full_r2['marginal']),
    Delta_marginal_R2 = delta_marginal,
    Conditional_R2 = unname(full_r2['conditional'])
  )
  vc <- as.data.frame(lme4::VarCorr(fit))
  random <- tibble::tibble(
    Fraction = fraction,
    Cruise = cruise,
    Shelf_region = shelf,
    Random_effect = ifelse(vc$grp == 'Residual', 'Residual', paste0(vc$grp, ': ', vc$var1)),
    Variance = vc$vcov,
    SD = vc$sdcor
  )
  list(fixed = fixed, random = random)
}

results <- purrr::map(model_files, extract_model)
fixed_table <- purrr::map_dfr(results, 'fixed')
random_table <- purrr::map_dfr(results, 'random')
order_rows <- function(x) {
  x %>%
    dplyr::mutate(
      Fraction = factor(Fraction, levels = c('Seston', 'Plankton')),
      Cruise = factor(Cruise, levels = c('w22', 's22', 'w23', 's23')),
      Shelf_region = factor(Shelf_region, levels = c('off_shelf', 'on_shelf'))
    ) %>%
    dplyr::arrange(Fraction, Cruise, Shelf_region) %>%
    dplyr::mutate(dplyr::across(c(Fraction, Cruise, Shelf_region), as.character))
}
fixed_table <- order_rows(fixed_table)
random_table <- order_rows(random_table)
fixed_numeric <- c('Estimate', 'Std_error', 'df', 't_value', 'p', 'Baseline_marginal_R2', 'Marginal_R2', 'Delta_marginal_R2', 'Conditional_R2')
fixed_table <- fixed_table %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(fixed_numeric), format_four_sig),
    p = dplyr::if_else(!is.na(p) & p == '0', '<0.001', p)
  )
random_table <- random_table %>%
  dplyr::mutate(dplyr::across(c('Variance', 'SD'), format_four_sig))
fixed_file <- file.path(output_dir, 'SPECTRA_biovolume_vs_slope_fixed_effects_4sig.csv')
random_file <- file.path(output_dir, 'SPECTRA_biovolume_vs_slope_random_effects_4sig.csv')
readr::write_csv(fixed_table, fixed_file, na = 'NA')
readr::write_csv(random_table, random_file, na = 'NA')
message('Read ', length(model_files), ' models.')
message('Fixed effects: ', fixed_file)
message('Random effects: ', random_file)
