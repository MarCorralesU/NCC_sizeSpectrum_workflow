#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(optparse)
  library(gt)
})

setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
base_dir <- getwd()

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------
option_list <- list(
  make_option(
    c("-s", "--sizeVar"),
    type = "character",
    default = "contour",
    help = "size variable dir used in ISIIS size spectra path"
  ),
  make_option(
    c("-o", "--outputDir"),
    type = "character",
    default = NULL,
    help = "Optional output directory"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

output_dir <- if (is.null(opt$outputDir)) {
  file.path(base_dir, "tables", "size_spectra_coefficients")
} else {
  opt$outputDir
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Constants: preserve figure order
# ------------------------------------------------------------
cruise_levels   <- c("w22", "s22", "w23", "s23")
transect_levels <- c("GH", "CR", "CM", "NH", "HH", "RR")
source_levels   <- c("allParticles", "plankton")

x1_full <- 12.78849
x2_full <- 4.547442e13

# ------------------------------------------------------------
# Helper formatting
# ------------------------------------------------------------
fmt_slope <- function(med, sd) {
  ifelse(
    is.na(med),
    NA_character_,
    paste0(formatC(med, format = "f", digits = 2),
           " \u00B1 ",
           formatC(sd, format = "f", digits = 2))
  )
}

fmt_intercept <- function(med, q1, q3) {
  ifelse(
    is.na(med),
    NA_character_,
    paste0(formatC(med, format = "e", digits = 2),
           " [",
           formatC(q1, format = "e", digits = 2),
           "\u2013",
           formatC(q3, format = "e", digits = 2),
           "]")
  )
}

fmt_biovol <- function(med, q1, q3) {
  ifelse(
    is.na(med),
    NA_character_,
    paste0(
      formatC(med, format = "e", digits = 2),
      " [",
      formatC(q1, format = "e", digits = 2),
      "\u2013",
      formatC(q3, format = "e", digits = 2),
      "]"
    )
  )
}

# ------------------------------------------------------------
# Helper to safely read and bind files
# ------------------------------------------------------------
read_many_csv <- function(files) {
  
  if (length(files) == 0) return(tibble())
  
  purrr::map2_dfr(
    files,
    seq_along(files),
    ~ readr::read_csv(.x, show_col_types = FALSE) %>%
      dplyr::mutate(source_file_id = .y)
  ) %>%
    dplyr::rename(group_id_old = group_id) %>%
    dplyr::group_by(source_file_id, group_id_old) %>%
    dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
    dplyr::ungroup()
}

calculate_total_biovolume <- function(dat) {
  
  required <- c(
    "group_id", "logNB", "biovol_um3_width",
    "slope", "intercept"
  )
  
  if (nrow(dat) == 0 || !all(required %in% names(dat))) {
    return(dat)
  }
  
  dat %>%
    dplyr::mutate(
      logNB = as.numeric(logNB),
      biovol_um3_width = as.numeric(biovol_um3_width),
      slope = as.numeric(slope),
      intercept = as.numeric(intercept),
      biovolume_um3_m3 = (10^logNB) * biovol_um3_width
    ) %>%
    dplyr::group_by(group_id) %>%
    dplyr::summarise(
      slope = dplyr::first(slope),
      intercept = dplyr::first(intercept),
      total_biovolume_um3_m3 = sum(biovolume_um3_m3, na.rm = TRUE),
      .groups = "drop"
    ) 
}

# ------------------------------------------------------------
# One-panel coefficient loader
# ------------------------------------------------------------
load_panel_coeffs <- function(Cruise, Transect, size_var = "contour") {
  
  # ---- plankton files (same logic as figure script)
  plankton_dir <- file.path(
    base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded",
    size_var, Cruise, Transect
  )
  
  plankton_files <- list.files(
    plankton_dir,
    pattern = "by_community\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  dat_plankton <- read_many_csv(plankton_files)
  
  # ---- allParticles coeffs used in Figure 3 come from merged LISST-ISIIS data
  merged_dir <- file.path(
    base_dir, "data", "ISIIS-LISST_mergedNBSS_5m_allParticles"
  )
  
  merged_files_all <- list.files(
    merged_dir,
    pattern = "LISST-ISIIS_NBSS\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Try to restrict to this panel using filename
  merged_files_panel <- list.files(
    merged_dir,
    pattern = paste0("^", Cruise, "_", Transect, "_LISST-ISIIS_NBSS\\.csv$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  dat_allparticles <- read_many_csv(merged_files_panel)
  
  # Fallback: if filename filtering failed but columns exist, subset after reading all files
  if (nrow(dat_allparticles) == 0 && length(merged_files_all) > 0) {
    dat_allparticles <- read_many_csv(merged_files_all)
    
    possible_cruise_cols <- c("Cruise", "cruise", "cruise_id")
    possible_transect_cols <- c("Transect_ID", "transect_id", "Transect", "transect")
    
    cruise_col <- possible_cruise_cols[possible_cruise_cols %in% names(dat_allparticles)][1]
    transect_col <- possible_transect_cols[possible_transect_cols %in% names(dat_allparticles)][1]
    
    if (!is.na(cruise_col) && !is.na(transect_col)) {
      dat_allparticles <- dat_allparticles %>%
        filter(.data[[cruise_col]] == Cruise, .data[[transect_col]] == Transect)
    } else {
      dat_allparticles <- tibble()
    }
    
  }

  dat_plankton <- calculate_total_biovolume(dat_plankton)
  
  dat_allparticles <- calculate_total_biovolume(dat_allparticles)
  
  stopifnot(
    !anyDuplicated(dat_plankton$group_id),
    !anyDuplicated(dat_allparticles$group_id)
  )
  
  # ----------------------------------------------------------
  # Build summaries
  # NOTE:
  #   intercept is assumed to be in log10 space, as in your figure script,
  #   so "untransformed intercept" = 10^intercept
  # ----------------------------------------------------------
  summarize_coeffs <- function(dat, source_name) {
    if (nrow(dat) == 0 || !all(c("slope", "intercept") %in% names(dat))) {
      return(tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        source = source_name,
        slope_median = NA_real_,
        slope_sd = NA_real_,
        intercept_median_untrans = NA_real_,
        intercept_q1_untrans = NA_real_,
        intercept_q3_untrans = NA_real_,
        total_biovol_median =  NA_real_,
        total_biovol_q1 = NA_real_,
        total_biovol_q3 = NA_real_,
        biovol_median = NA_real_,
        biovol_q1 = NA_real_,
        biovol_q3 = NA_real_
      ))
    }
    
    slope_vals <- suppressWarnings(as.numeric(dat$slope))
    int_vals   <- suppressWarnings(as.numeric(dat$intercept))
    
    keep <- is.finite(slope_vals) & is.finite(int_vals)
    slope_vals <- slope_vals[keep]
    int_vals   <- int_vals[keep]
    
    int_untrans <- 10^int_vals
    
    # theoretical integrated biovolume across full common size range
    biovol_vals <- rep(NA_real_, length(slope_vals))
    if (length(slope_vals) == length(int_vals) && length(slope_vals) > 0) {
      b_lin <- 10^int_vals
      
      biovol_vals <- ifelse(
        abs(slope_vals + 1) < 1e-12,
        b_lin * log(x2_full / x1_full),
        (b_lin / (slope_vals + 1)) * (x2_full^(slope_vals + 1) - x1_full^(slope_vals + 1))
      )
      
      biovol_vals <- biovol_vals[is.finite(biovol_vals)]
    } else {
      biovol_vals <- numeric(0)
    }
    
    total_biovol_vals <- suppressWarnings(
      as.numeric(dat$total_biovolume_um3_m3)
    )
    
    total_biovol_vals <- total_biovol_vals[keep]
    total_biovol_vals <- total_biovol_vals[is.finite(total_biovol_vals)]
    
    tibble(
      Cruise = Cruise,
      Transect_ID = Transect,
      source = source_name,
      slope_median = if (length(slope_vals) > 0) median(slope_vals, na.rm = TRUE) else NA_real_,
      slope_sd = if (length(slope_vals) > 1) sd(slope_vals, na.rm = TRUE) else NA_real_,
      intercept_median_untrans = if (length(int_untrans) > 0) median(int_untrans, na.rm = TRUE) else NA_real_,
      intercept_q1_untrans = if (length(int_untrans) > 0) quantile(int_untrans, 0.25, na.rm = TRUE, names = FALSE) else NA_real_,
      intercept_q3_untrans = if (length(int_untrans) > 0) quantile(int_untrans, 0.75, na.rm = TRUE, names = FALSE) else NA_real_,
      total_biovol_median = if (length(total_biovol_vals) > 0) median(total_biovol_vals) else NA_real_,
      total_biovol_q1 = if (length(total_biovol_vals) > 0) quantile(total_biovol_vals, 0.25, names = FALSE) else NA_real_,
      total_biovol_q3 = if (length(total_biovol_vals) > 0) quantile(total_biovol_vals, 0.75, names = FALSE) else NA_real_,
      biovol_median = if (length(biovol_vals) > 0) median(biovol_vals, na.rm = TRUE) else NA_real_,
      biovol_q1 = if (length(biovol_vals) > 0) quantile(biovol_vals, 0.25, na.rm = TRUE, names = FALSE) else NA_real_,
      biovol_q3 = if (length(biovol_vals) > 0) quantile(biovol_vals, 0.75, na.rm = TRUE, names = FALSE) else NA_real_
    )
  }
  
  bind_rows(
    summarize_coeffs(dat_allparticles, "allParticles"),
    summarize_coeffs(dat_plankton, "plankton")
  )
}

# ------------------------------------------------------------
# Run across all panels in figure order
# ------------------------------------------------------------
panel_grid <- expand.grid(
  Cruise = cruise_levels,
  Transect_ID = transect_levels,
  stringsAsFactors = FALSE
)

  coef_summary_long <- purrr::pmap_dfr(
    list(panel_grid$Cruise, panel_grid$Transect_ID),
    ~ load_panel_coeffs(..1, ..2, size_var = opt$sizeVar)
  ) %>%
    mutate(
      Cruise = factor(Cruise, levels = cruise_levels),
      Transect_ID = factor(Transect_ID, levels = transect_levels),
      source = factor(source, levels = source_levels),
      slope_summary = fmt_slope(slope_median, slope_sd),
      intercept_summary = fmt_intercept(
        intercept_median_untrans,
        intercept_q1_untrans,
        intercept_q3_untrans
      ),
      total_biovol_summary = fmt_biovol(
        total_biovol_median,
        total_biovol_q1,
        total_biovol_q3
      ),
      biovol_summary = fmt_biovol(
        biovol_median,
        biovol_q1,
        biovol_q3
      )
    ) %>%
    arrange(Transect_ID, source, Cruise)

# ------------------------------------------------------------
# Wide manuscript-style table
# rows = transect + source subrows
# cols = cruises (slope/intercept under each)
# ------------------------------------------------------------
coef_summary_wide <- coef_summary_long %>%
  select(
    Transect_ID, source, Cruise,
    slope_summary, intercept_summary,
    biovol_summary, total_biovol_summary
  ) %>%
  pivot_wider(
    names_from = Cruise,
    values_from = c(
      slope_summary, intercept_summary,
      biovol_summary, total_biovol_summary
    ),
    names_glue = "{Cruise}_{.value}"
  )

# ------------------------------------------------------------
# Optional: nicer column names
# ------------------------------------------------------------
  coef_summary_wide_pretty <- coef_summary_wide %>%
    rename(
      Source = source,
      `Winter 2022 slope` = w22_slope_summary,
      `Winter 2022 intercept` = w22_intercept_summary,
      `Winter 2022 total biovolume` = w22_total_biovol_summary,
      `Winter 2022 integrated biovolume` = w22_biovol_summary,

      `Summer 2022 slope` = s22_slope_summary,
      `Summer 2022 intercept` = s22_intercept_summary,
      `Summer 2022 total biovolume` = s22_total_biovol_summary,
      `Summer 2022 integrated biovolume` = s22_biovol_summary,
      
      `Winter 2023 slope` = w23_slope_summary,
      `Winter 2023 intercept` = w23_intercept_summary,
      `Winter 2023 total biovolume` = w23_total_biovol_summary,
      `Winter 2023 integrated biovolume` = w23_biovol_summary,
      
      `Summer 2023 slope` = s23_slope_summary,
      `Summer 2023 intercept` = s23_intercept_summary,
      `Summer 2023 total biovolume` = s23_total_biovol_summary,
      `Summer 2023 integrated biovolume` = s23_biovol_summary
    )

##gt table:
# ------------------------------------------------------------
# Create GT table
# ------------------------------------------------------------
  coef_gt <- coef_summary_wide_pretty %>%
    gt(
      rowname_col = "Transect_ID",
      groupname_col = "Source"
    ) %>%
    tab_spanner(
      label = "Winter 2022",
      columns = c(
        `Winter 2022 slope`,
        `Winter 2022 intercept`,
        `Winter 2022 total biovolume`,
        `Winter 2022 integrated biovolume`
      )
    ) %>%
    tab_spanner(
      label = "Summer 2022",
      columns = c(
        `Summer 2022 slope`,
        `Summer 2022 intercept`,
        `Summer 2022 total biovolume`,
        `Summer 2022 integrated biovolume`
      )
    ) %>%
    tab_spanner(
      label = "Winter 2023",
      columns = c(
        `Winter 2023 slope`,
        `Winter 2023 intercept`,
        `Winter 2023 total biovolume`,
        `Winter 2023 integrated biovolume`
      )
    ) %>%
    tab_spanner(
      label = "Summer 2023",
      columns = c(
        `Summer 2023 slope`,
        `Summer 2023 intercept`,
        `Summer 2023 total biovolume`,
        `Summer 2023 integrated biovolume`
      )
    ) %>%
    cols_align(
      align = "center",
      everything()
    ) %>%
    cols_width(
      everything() ~ px(120)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) %>%
    tab_options(
      table.font.size = px(13),
      column_labels.font.weight = "bold"
    )



# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
readr::write_csv(
  coef_summary_long,
  file.path(output_dir, "size_spectra_coefficients_summary_long.csv"),
  na = ""
)

readr::write_csv(
  coef_summary_wide_pretty,
  file.path(output_dir, "size_spectra_coefficients_summary_wide.csv"),
  na = ""
)

gtsave(
  coef_gt,
  file.path(output_dir, "size_spectra_coefficients_table.html")
)

message("Wrote:")
message("  ", file.path(output_dir, "size_spectra_coefficients_summary_long.csv"))
message("  ", file.path(output_dir, "size_spectra_coefficients_summary_wide.csv"))