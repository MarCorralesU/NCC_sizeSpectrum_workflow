#!/usr/bin/env Rscript
# Summarize fitted environmental LMMs without refitting.
# Run from the SPECTRA project root.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(lme4)
  library(lmerTest)
})

if (!requireNamespace("performance", quietly = TRUE)) {
  stop("Install the 'performance' package to calculate marginal R2.")
}

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else file.path(
  ".", "figures", "size_spectra_mixed_models_cast_ID", "envVars_slope_LMM"
)
output_dir <- file.path(input_dir, 'summary_tables')
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(
  input_dir,
  pattern = "^SPECTRA_(Seston|Plankton)_pooled_shelf_(Temperature|log10_chl_a_ug_l|Oxygen_ml_l)_model\\.rds$",
  full.names = TRUE
)
if (!length(model_files)) stop("No matching environmental LMM .rds files found in: ", input_dir)
if (length(model_files) != 6) warning("Expected six models; found ", length(model_files))

predictor_labels <- c(
  Temperature = "Temperature",
  log10_chl_a_ug_l = "log10 chlorophyll-a",
  Oxygen_ml_l = "Dissolved oxygen"
)

# Four significant figures. Scientific notation only outside [0.001, 1000].
format_four_sig <- function(x) {
  vapply(x, function(value) {
    if (is.na(value) || !is.finite(value)) return(NA_character_)
    if (value == 0) return("0")
    if (abs(value) < 0.001 || abs(value) > 1000) {
      out <- formatC(value, format = "e", digits = 3)
      return(sub("e([+-])0+([0-9]+)$", "e\\1\\2", out))
    }
    format(signif(value, 4), scientific = FALSE, trim = TRUE)
  }, character(1))
}

extract_model <- function(path) {
  filename <- basename(path)
  stem <- sub("_model\\.rds$", "", filename)
  fraction <- sub("^SPECTRA_(Seston|Plankton)_pooled_shelf_.*$", "\\1", stem)
  predictor <- sub("^SPECTRA_(Seston|Plankton)_pooled_shelf_", "", stem)
  obj <- readRDS(path)
  fit <- if (inherits(obj, "merMod")) obj else obj$model
  if (!inherits(fit, "merMod")) stop("No lme4 mixed model found in: ", path)

  coefs <- as.data.frame(coef(summary(fit)))
  coefs$Fixed_effect_predictor <- rownames(coefs)
  required <- c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")
  if (!all(required %in% names(coefs))) {
    stop("Missing lmerTest coefficient columns in ", filename,
         ". Ensure lmerTest is installed and the model supports Satterthwaite df.")
  }
  r2 <- suppressWarnings(performance::r2_nakagawa(fit))
  marginal_r2 <- as.numeric(r2$R2_marginal)
  
  baseline_fit <- obj$baseline
  baseline_r2 <- suppressWarnings(performance::r2_nakagawa(baseline_fit))
  baseline_marginal_r2 <- as.numeric(baseline_r2$R2_marginal)
  
  delta_marginal_r2 <- marginal_r2 - baseline_marginal_r2

  fixed <- tibble::tibble(
    Environmental_predictor_variable = unname(predictor_labels[[predictor]]),
    Fraction = fraction,
    Fixed_effect_predictor = coefs$Fixed_effect_predictor,
    Estimate = coefs$Estimate,
    Std_error = coefs[["Std. Error"]],
    df = coefs[["df"]],
    t_value = coefs[["t value"]],
    p = coefs[["Pr(>|t|)"]],
    Baseline_marginal_R2 = baseline_marginal_r2,
    Marginal_R2 = marginal_r2,
    Delta_marginal_R2 = delta_marginal_r2
  )

  vc <- as.data.frame(lme4::VarCorr(fit))
  random <- tibble::tibble(
    Environmental_predictor_variable = unname(predictor_labels[[predictor]]),
    Fraction = fraction,
    Random_effect = ifelse(vc$grp == "Residual", "Residual", paste0(vc$grp, ": ", vc$var1)),
    Variance = vc$vcov,
    SD = vc$sdcor
  )

  list(fixed = fixed, random = random)
}

results <- purrr::map(model_files, extract_model)
fixed_table <- purrr::map_dfr(results, "fixed")
random_table <- purrr::map_dfr(results, "random")

# Order the output by fraction, environmental predictor, and effect.
order_rows <- function(dat) {
  dat %>%
    dplyr::mutate(
      Fraction = factor(Fraction, levels = c("Seston", "Plankton")),
      Environmental_predictor_variable = factor(
        Environmental_predictor_variable,
        levels = unname(predictor_labels)
      )
    ) %>%
    dplyr::arrange(Fraction, Environmental_predictor_variable) %>%
    dplyr::mutate(
      Fraction = as.character(Fraction),
      Environmental_predictor_variable = as.character(Environmental_predictor_variable)
    )
}

fixed_table <- order_rows(fixed_table)
random_table <- order_rows(random_table)

fixed_numeric <- c(
  "Estimate", "Std_error", "df", "t_value", "p",
  "Baseline_marginal_R2", "Marginal_R2", "Delta_marginal_R2"
)
random_numeric <- c("Variance", "SD")
fixed_table <- fixed_table %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(fixed_numeric), format_four_sig),
    p = dplyr::if_else(p == "0", "<0.001", p)
  )

random_table <- random_table %>% dplyr::mutate(dplyr::across(dplyr::all_of(random_numeric), format_four_sig))

fixed_file <- file.path(output_dir, "SPECTRA_environmental_LMM_fixed_effects_4sig.csv")
random_file <- file.path(output_dir, "SPECTRA_environmental_LMM_random_effects_4sig.csv")
readr::write_csv(fixed_table, fixed_file, na = "NA")
readr::write_csv(random_table, random_file, na = "NA")
message("Read ", length(model_files), " fitted models.")
message("Saved fixed effects: ", fixed_file)
message("Saved random effects: ", random_file)
