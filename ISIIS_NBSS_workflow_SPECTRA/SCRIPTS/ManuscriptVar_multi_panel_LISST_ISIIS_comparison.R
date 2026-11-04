#!/usr/bin/env Rscript

# ============================================================
# Multi-panel LISST vs ISIIS size-spectrum comparison
#
# Outputs:
#   1) LISST vs ISIIS_allParticles
#   2) LISST vs ISIIS_plankton
#   3) median-only comparison with:
#        LISST
#        ISIIS_allParticles
#        ISIIS_plankton
#
# Rows    = Transect_ID: GH, CR, CM, NH, HH, RR
# Columns = Cruise:      w22, s22, w23, s23
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(optparse)
  library(lubridate)
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
  ),
  make_option(
    c("--withEq"),
    action = "store_true",
    default = FALSE,
    help = "Add per-panel fitted-equation annotations to the two 2-source figures"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ------------------------------------------------------------
# Constants
# ------------------------------------------------------------
cruise_levels   <- c("w22", "s22", "w23", "s23")
transect_levels <- c("GH", "CR", "CM", "NH", "HH", "RR")

cruise_label_map <- c(
  w22 = "Winter 2022",
  s22 = "Summer 2022",
  w23 = "Winter 2023",
  s23 = "Summer 2023"
)

transect_label_map <- c(
  GH = "Grays Harbor",
  CR = "Columbia River",
  CM = "Cape Meares",
  NH = "Newport Line",
  HH = "Heceta Head",
  RR = "Rogue River"
)

src_cols <- c(
  LISST = "seagreen4",
  ISIIS_allParticles = "royalblue3",
  ISIIS_plankton = "darkorange"
)

output_dir <- if (is.null(opt$outputDir)) {
  file.path(base_dir, "figures", "multi_panel_size_spectra_comparison")
} else {
  opt$outputDir
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper formatting
# ------------------------------------------------------------
fmt_num <- function(x) formatC(x, format = "fg", digits = 3)
fmt_plotmath_sci <- function(x, digits = 2) {
  if (!is.finite(x) || x <= 0) return(NA_character_)
  expn <- floor(log10(x))
  mant <- signif(x / 10^expn, digits = digits)
  paste0(mant, "%*%10^", expn)
}

biovol_to_ecd <- function(v_um3) {
  ((6 * v_um3) / pi)^(1/3)
}

# ------------------------------------------------------------
# One-panel loader/prepper
# ------------------------------------------------------------
load_panel_data <- function(Cruise, Transect, size_var = "contour") {
  
  isiis_dir_particles <- file.path(
    base_dir, "data", "size_spectra_allParticles_5-10m-bins_thresholded",
    size_var, Cruise, Transect
  )
  
  isiis_dir_particles_LISST_ISIIS <- file.path(
    base_dir, "data", "ISIIS-LISST_mergedNBSS_5m_allParticles")
  
  isiis_dir_plankton <- file.path(
    base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded",
    size_var, Cruise, Transect
  )
  
  isiis_files_particles <- list.files(
    isiis_dir_particles,
    pattern = "by_community\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  isiis_files_particles_LISST_ISIIS <- list.files(
    isiis_dir_particles_LISST_ISIIS,
    pattern = paste0("^", Cruise, "_", Transect, "_LISST-ISIIS_NBSS\\.csv$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  message("Merged file(s) for ", Cruise, "/", Transect, ":")
  print(basename(isiis_files_particles_LISST_ISIIS))
  
  message("----")
  message("Panel: ", Cruise, " / ", Transect)
  message("Merged LISST-ISIIS files found: ", length(isiis_files_particles_LISST_ISIIS))
  print(basename(isiis_files_particles_LISST_ISIIS)[1:min(10, length(isiis_files_particles_LISST_ISIIS))])
  
  isiis_files_plankton <- list.files(
    isiis_dir_plankton,
    pattern = "by_community\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  lisst_file <- file.path(
    base_dir, "data", "LISST_data", "LISST_size_spectra_thresholded",
    "5m_depthBin",
    paste0("SPECTRA_", Cruise, "_", Transect, "_5m-depthBin_LISST_size_spectra.csv")
  )
  
  if ((length(isiis_files_particles) == 0 && length(isiis_files_plankton) == 0) || !file.exists(lisst_file)) {
    message("Skipping missing panel: ", Cruise, " / ", Transect)
    return(NULL)
  }
  
  # -------------------------
  # Load ISIIS allParticles
  # -------------------------
  dat_ISIIS_particles <- if (length(isiis_files_particles) > 0) {
    purrr::map_dfr(isiis_files_particles, ~ readr::read_csv(.x, show_col_types = FALSE)) %>%
      mutate(
        time_Pacific_orig = time_Pacific,
        mean_time_Pacific = time_Pacific - hours(1)
      ) %>%
      dplyr::rename(group_id_old = group_id) %>%
      dplyr::mutate(
        group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE))
      )
  } else {
    tibble()
  }
  
  # -------------------------
  # Load LISST-ISIIS allParticles
  # -------------------------
  dat_ISIIS_particles_LISST_ISIIS <- if (length(isiis_files_particles_LISST_ISIIS) > 0) {
    purrr::map_dfr(isiis_files_particles_LISST_ISIIS, ~ readr::read_csv(.x, show_col_types = FALSE)) %>%
      mutate(
        time_Pacific_orig = time_Pacific,
        total_biovolume = biovolume_um3
      ) %>%
      dplyr::rename(group_id_old = group_id) %>%
      dplyr::mutate(
        group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE))
      )
  } else {
    tibble()
  }
  
  # -------------------------
  # Load ISIIS plankton
  # -------------------------
  dat_ISIIS_plankton <- if (length(isiis_files_plankton) > 0) {
    purrr::map_dfr(isiis_files_plankton, ~ readr::read_csv(.x, show_col_types = FALSE)) %>%
      mutate(
        time_Pacific_orig = time_Pacific,
        mean_time_Pacific = time_Pacific - hours(1)
      ) %>%
      dplyr::rename(group_id_old = group_id) %>%
      dplyr::mutate(
        group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE))
      )
  } else {
    tibble()
  }
  
  # -------------------------
  # Load LISST
  # -------------------------
  dat_LISST <- readr::read_csv(lisst_file, show_col_types = FALSE) %>%
    dplyr::filter(slope < 1)
  
  # -------------------------
  # Minimal plotting frames
  # -------------------------
   message("Columns in dat_LISST:")
   print(names(dat_LISST))

   message("Columns in dat_ISIIS_particles:")
   print(names(dat_ISIIS_particles))

   message("Columns in dat_ISIIS_plankton:")
   print(names(dat_ISIIS_plankton))


   lisst <- dat_LISST %>%
    dplyr::transmute(
      source   = "LISST",
      group_id = group_id,
      size_num = as.numeric(size_class_mid_um),
      NB       = as.numeric(NB)
    ) %>%
    dplyr::filter(
      is.finite(size_num), size_num > 0,
      is.finite(NB), NB > 0
    )
  
  isiis_particles <- dat_ISIIS_particles %>%
    dplyr::transmute(
      source   = "ISIIS_allParticles",
      group_id = group_id,
      size_num = as.numeric(ECD_um_mid),
      NB       = as.numeric(NB)
    ) %>%
    dplyr::filter(
      is.finite(size_num), size_num > 0,
      is.finite(NB), NB > 0
    )
  
  isiis_plankton <- dat_ISIIS_plankton %>%
    dplyr::transmute(
      source   = "ISIIS_plankton",
      group_id = group_id,
      size_num = as.numeric(ECD_um_mid),
      NB       = as.numeric(NB)
    ) %>%
    dplyr::filter(
      is.finite(size_num), size_num > 0,
      is.finite(NB), NB > 0
    )
  
  if (nrow(lisst) == 0 && nrow(isiis_particles) == 0 && nrow(isiis_plankton) == 0) {
    message("Skipping empty panel: ", Cruise, " / ", Transect)
    return(NULL)
  }
  
  # -------------------------
  # Shared ordered size axis
  # -------------------------
  levels_x <- sort(unique(c(
    lisst$size_num,
    isiis_particles$size_num,
    isiis_plankton$size_num
  )))
  
  lisst <- lisst %>%
    dplyr::mutate(size_fac = factor(size_num, levels = levels_x))
  
  isiis_particles <- isiis_particles %>%
    dplyr::mutate(size_fac = factor(size_num, levels = levels_x))
  
  isiis_plankton <- isiis_plankton %>%
    dplyr::mutate(size_fac = factor(size_num, levels = levels_x))
  
  pts_df <- dplyr::bind_rows(lisst, isiis_particles, isiis_plankton) %>%
    dplyr::arrange(source, group_id, size_fac) %>%
    dplyr::mutate(
      Cruise = Cruise,
      Transect_ID = Transect
    )
  
  med_df <- pts_df %>%
    dplyr::group_by(Cruise, Transect_ID, source, size_fac) %>%
    dplyr::summarise(
      NB_med = median(NB, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(size_num = as.numeric(as.character(size_fac))) %>%
    dplyr::arrange(Cruise, Transect_ID, source, size_fac)
  
  # -------------------------
  # Figure 3 only:
  # median-equation labels + fitted lines
  # allParticles coefficients come from LISST-ISIIS merged dataset
  # plankton coefficients come from ISIIS plankton dataset
  # coefficients are in biovolume space, so convert fit x to ECD for plotting
  # -------------------------
  
  particles_fit_b_log10 <- suppressWarnings(
    median(as.numeric(dat_ISIIS_particles_LISST_ISIIS$intercept), na.rm = TRUE)
  )
  particles_fit_s <- suppressWarnings(
    median(as.numeric(dat_ISIIS_particles_LISST_ISIIS$slope), na.rm = TRUE)
  )
  
  message(
    "Panel ", Cruise, "/", Transect,
    " | merged allParticles median intercept(log10) = ", signif(particles_fit_b_log10, 4),
    " | slope = ", signif(particles_fit_s, 4)
  )
  
  plankton_fit_b_log10 <- suppressWarnings(
    median(as.numeric(dat_ISIIS_plankton$intercept), na.rm = TRUE)
  )
  plankton_fit_s <- suppressWarnings(
    median(as.numeric(dat_ISIIS_plankton$slope), na.rm = TRUE)
  )
  
  particles_fit_b_untrans <- 10^particles_fit_b_log10
  plankton_fit_b_untrans  <- 10^plankton_fit_b_log10
  
  fit_df_medonly <- tibble::tibble()
  eq_df_medonly  <- tibble::tibble()
  
  # fixed BIovolume range used for the power-law fits
  v1 <- 12.78849
  v2 <- 4.547442e13
  
  v_seq <- 10^seq(log10(v1), log10(v2), length.out = 200)
  ecd_seq <- biovol_to_ecd(v_seq)
  
  #to control equation positions
  label_y_top <- 4e9
  label_gap_factor <- 18
  label_y_bottom <- label_y_top / label_gap_factor
  
  # allParticles fit from LISST-ISIIS merged coeffs
  # allParticles fit from LISST-ISIIS merged coeffs
  if (nrow(dat_ISIIS_particles_LISST_ISIIS) > 0 &&
      is.finite(particles_fit_b_untrans) && is.finite(particles_fit_s)) {
    
    fit_df_medonly <- dplyr::bind_rows(
      fit_df_medonly,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        fit_source = "allParticles_fit",
        size_num = ecd_seq,
        NB_fit = particles_fit_b_untrans * v_seq^particles_fit_s
      )
    )
    
    eq_df_medonly <- dplyr::bind_rows(
      eq_df_medonly,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        fit_source = "allParticles_fit",
        x = max(ecd_seq, na.rm = TRUE) * 0.98,
        y = label_y_top,
        label = paste0(
          "NB==", fmt_plotmath_sci(signif(particles_fit_b_untrans, 3)),
          " * (mu*m^{3})^{", fmt_num(round(particles_fit_s, 2)), "}"
        )
      )
    )
  }
  
  # plankton fit
  # plankton fit
  if (nrow(dat_ISIIS_plankton) > 0 &&
      is.finite(plankton_fit_b_untrans) && is.finite(plankton_fit_s)) {
    
    fit_df_medonly <- dplyr::bind_rows(
      fit_df_medonly,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        fit_source = "plankton_fit",
        size_num = ecd_seq,
        NB_fit = plankton_fit_b_untrans * v_seq^plankton_fit_s
      )
    )
    
    eq_df_medonly <- dplyr::bind_rows(
      eq_df_medonly,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        fit_source = "plankton_fit",
        x = max(ecd_seq, na.rm = TRUE) * 0.98,
        y = label_y_bottom,
        label = paste0(
          "NB==", fmt_plotmath_sci(signif(plankton_fit_b_untrans, 3)),
          " * (mu*m^{3})^{", fmt_num(round(plankton_fit_s, 2)), "}"
        )
      )
    )
  }
  
  
  # -------------------------
  # Equation labels for 2-source figures
  # -------------------------
  lisst_b <- suppressWarnings(median(as.numeric(dat_LISST$intercept), na.rm = TRUE))
  lisst_s <- suppressWarnings(median(as.numeric(dat_LISST$slope),     na.rm = TRUE))
  
  particles_b <- suppressWarnings(median(as.numeric(dat_ISIIS_particles$intercept), na.rm = TRUE))
  particles_s <- suppressWarnings(median(as.numeric(dat_ISIIS_particles$slope),     na.rm = TRUE))
  
  plankton_b <- suppressWarnings(median(as.numeric(dat_ISIIS_plankton$intercept), na.rm = TRUE))
  plankton_s <- suppressWarnings(median(as.numeric(dat_ISIIS_plankton$slope),     na.rm = TRUE))
  
  eq_df <- tibble::tibble()
  
  if (nrow(pts_df %>% filter(source == "LISST")) > 0 &&
      is.finite(lisst_b) && is.finite(lisst_s)) {
    
    lisst_y <- stats::quantile(
      pts_df$NB[pts_df$source == "LISST"], 0.85, na.rm = TRUE
    )
    
    eq_df <- bind_rows(
      eq_df,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        source = "LISST",
        x = 1000,
        y = lisst_y,
        label = paste0("y==", fmt_num(lisst_b), " * x^{", fmt_num(lisst_s), "}")
      )
    )
  }
  
  if (nrow(pts_df %>% filter(source == "ISIIS_allParticles")) > 0 &&
      is.finite(particles_b) && is.finite(particles_s)) {
    
    particles_y <- stats::quantile(
      pts_df$NB[pts_df$source == "ISIIS_allParticles"], 0.20, na.rm = TRUE
    )
    
    eq_df <- bind_rows(
      eq_df,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        source = "ISIIS_allParticles",
        x = 70,
        y = particles_y,
        label = paste0("y==", fmt_num(particles_b), " * x^{", fmt_num(particles_s), "}")
      )
    )
  }
  
  if (nrow(pts_df %>% filter(source == "ISIIS_plankton")) > 0 &&
      is.finite(plankton_b) && is.finite(plankton_s)) {
    
    plankton_y <- stats::quantile(
      pts_df$NB[pts_df$source == "ISIIS_plankton"], 0.20, na.rm = TRUE
    )
    
    eq_df <- bind_rows(
      eq_df,
      tibble::tibble(
        Cruise = Cruise,
        Transect_ID = Transect,
        source = "ISIIS_plankton",
        x = 70,
        y = plankton_y,
        label = paste0("y==", fmt_num(plankton_b), " * x^{", fmt_num(plankton_s), "}")
      )
    )
  }
  
  list(
    pts_df = pts_df,
    med_df = med_df,
    eq_df  = eq_df,
    fit_df_medonly = fit_df_medonly,
    eq_df_medonly  = eq_df_medonly
  )
}

# ------------------------------------------------------------
# Loop all panels
# ------------------------------------------------------------
panel_grid <- expand.grid(
  Cruise = cruise_levels,
  Transect_ID = transect_levels,
  stringsAsFactors = FALSE
)

panel_data <- purrr::pmap(
  list(panel_grid$Cruise, panel_grid$Transect_ID),
  ~ load_panel_data(..1, ..2, size_var = opt$sizeVar)
)

pts_all <- purrr::map_dfr(panel_data, "pts_df")
med_all <- purrr::map_dfr(panel_data, "med_df")
eq_all  <- purrr::map_dfr(panel_data, "eq_df")

fit_medonly_all <- purrr::map_dfr(panel_data, "fit_df_medonly")
eq_medonly_all  <- purrr::map_dfr(panel_data, "eq_df_medonly")

fit_medonly_all <- fit_medonly_all %>%
  dplyr::mutate(
    Cruise = factor(Cruise, levels = cruise_levels),
    Transect_ID = factor(Transect_ID, levels = transect_levels),
    fit_source = factor(fit_source, levels = c("allParticles_fit", "plankton_fit"))
  )

eq_medonly_all <- eq_medonly_all %>%
  dplyr::mutate(
    Cruise = factor(Cruise, levels = cruise_levels),
    Transect_ID = factor(Transect_ID, levels = transect_levels),
    fit_source = factor(fit_source, levels = c("allParticles_fit", "plankton_fit"))
  )

# enforce facet order
pts_all <- pts_all %>%
  mutate(
    Cruise = factor(Cruise, levels = cruise_levels),
    Transect_ID = factor(Transect_ID, levels = transect_levels),
    source = factor(source, levels = c("LISST", "ISIIS_allParticles", "ISIIS_plankton"))
  )

med_all <- med_all %>%
  mutate(
    Cruise = factor(Cruise, levels = cruise_levels),
    Transect_ID = factor(Transect_ID, levels = transect_levels),
    source = factor(source, levels = c("LISST", "ISIIS_allParticles", "ISIIS_plankton"))
  )


eq_all <- eq_all %>%
  mutate(
    Cruise = factor(Cruise, levels = cruise_levels),
    Transect_ID = factor(Transect_ID, levels = transect_levels),
    source = factor(source, levels = c("LISST", "ISIIS_allParticles", "ISIIS_plankton"))
  )

make_panel_letters <- function(df, cruise_levels, transect_levels, x_val, y_val) {
  df %>%
    dplyr::distinct(Cruise, Transect_ID) %>%
    dplyr::mutate(
      Cruise = factor(as.character(Cruise), levels = cruise_levels),
      Transect_ID = factor(as.character(Transect_ID), levels = transect_levels)
    ) %>%
    dplyr::arrange(Transect_ID, Cruise) %>%
    dplyr::mutate(
      panel_letter = paste0(letters[dplyr::row_number()], "."),
      x = x_val,
      y = y_val
    )
}

panel_letter_df_medians <- make_panel_letters(
  df = med_all,
  cruise_levels = cruise_levels,
  transect_levels = transect_levels,
  x_val = 3.5,
  y_val = 1e10
)

# ------------------------------------------------------------
# Figure 1: LISST vs ISIIS_allParticles
# ------------------------------------------------------------
pts_particles <- pts_all %>%
  filter(source %in% c("LISST", "ISIIS_allParticles")) %>%
  mutate(source = factor(source, levels = c("LISST", "ISIIS_allParticles")))

panel_letter_df_particles <- make_panel_letters(
  df = pts_particles,
  cruise_levels = cruise_levels,
  transect_levels = transect_levels,
  x_val = 3.5,
  y_val = 1e11
)

med_particles <- med_all %>%
  filter(source %in% c("LISST", "ISIIS_allParticles")) %>%
  mutate(source = factor(source, levels = c("LISST", "ISIIS_allParticles")))

eq_particles <- eq_all %>%
  filter(source %in% c("LISST", "ISIIS_allParticles")) %>%
  mutate(source = factor(source, levels = c("LISST", "ISIIS_allParticles")))


p_particles <- ggplot() +
  geom_point(
    data = pts_particles,
    aes(x = size_num, y = NB, color = source),
    shape = 1,
    size = 0.7,
    stroke = 0.35,
    alpha = 0.25
  ) +
  geom_line(
    data = med_particles,
    aes(x = size_num, y = NB_med, color = source),
    linewidth = 0.8,
    linetype = "solid",
    alpha = 0.95
  ) +
  geom_point(
    data = med_particles,
    aes(x = size_num, y = NB_med, color = source),
    size = 0.5,
    alpha = 0.95
  ) +
  geom_text(
    data = panel_letter_df_particles,
    aes(x = x, y = y, label = panel_letter),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 6,
    fontface = "bold",
    color = "black"
  )  +
  scale_color_manual(
    values = src_cols[c("LISST", "ISIIS_allParticles")],
    name = "Source",
    labels = c("LISST", "ISIIS all particles")
  ) +
  scale_x_log10() +
  scale_y_log10(limits = c(1e-5, 1e11)) +
  facet_grid(
    Transect_ID ~ Cruise,
    drop = FALSE,
    labeller = labeller(
      Cruise = as_labeller(cruise_label_map),
      Transect_ID = as_labeller(transect_label_map)
    )
  ) +
  labs(
    x = expression(Equivalent~spherical~diameter~(mu*m)),
    y = expression(Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1)))

if (opt$withEq && nrow(eq_particles) > 0) {
  p_particles <- p_particles +
    geom_text(
      data = eq_particles %>% filter(source == "LISST"),
      aes(x = x, y = y, label = label),
      color = src_cols["LISST"],
      parse = TRUE,
      hjust = 1,
      size = 2.4,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = eq_particles %>% filter(source == "ISIIS_allParticles"),
      aes(x = x, y = y, label = label),
      color = src_cols["ISIIS_allParticles"],
      parse = TRUE,
      hjust = 0,
      size = 2.4,
      fontface = "bold",
      inherit.aes = FALSE
    )
}

# ------------------------------------------------------------
# Figure 2: LISST vs ISIIS_plankton
# ------------------------------------------------------------
pts_plankton <- pts_all %>%
  filter(source %in% c("LISST", "ISIIS_plankton")) %>%
  mutate(source = factor(source, levels = c("LISST", "ISIIS_plankton")))

panel_letter_df_plankton <- make_panel_letters(
  df = pts_plankton,
  cruise_levels = cruise_levels,
  transect_levels = transect_levels,
  x_val = 3.5,
  y_val = 1e11
)

med_plankton <- med_all %>%
  filter(source %in% c("LISST", "ISIIS_plankton")) %>%
  mutate(source = factor(source, levels = c("LISST", "ISIIS_plankton")))

eq_plankton <- eq_all %>%
  filter(source %in% c("LISST", "ISIIS_plankton")) %>%
  mutate(source = factor(source, levels = c("LISST", "ISIIS_plankton")))

p_plankton <- ggplot() +
  geom_point(
    data = pts_plankton,
    aes(x = size_num, y = NB, color = source),
    shape = 1,
    size = 0.7,
    stroke = 0.35,
    alpha = 0.25
  ) +
  geom_line(
    data = med_plankton,
    aes(x = size_num, y = NB_med, color = source),
    linewidth = 0.8,
    linetype = "solid",
    alpha = 0.95
  ) +
  geom_point(
    data = med_plankton,
    aes(x = size_num, y = NB_med, color = source),
    size = 0.5,
    alpha = 0.95
  ) +
  geom_text(
    data = panel_letter_df_plankton,
    aes(x = x, y = y, label = panel_letter),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  scale_color_manual(
    values = src_cols[c("LISST", "ISIIS_plankton")],
    name = "Source",
    labels = c("LISST", "ISIIS plankton")
  ) +
  scale_x_log10() +
  scale_y_log10(limits = c(1e-5, 1e11)) +
  facet_grid(
    Transect_ID ~ Cruise,
    drop = FALSE,
    labeller = labeller(
      Cruise = as_labeller(cruise_label_map),
      Transect_ID = as_labeller(transect_label_map)
    )
  ) +
  labs(
    x = expression(Equivalent~spherical~diameter~(mu*m)),
    y = expression(Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1)))

if (opt$withEq && nrow(eq_plankton) > 0) {
  p_plankton <- p_plankton +
    geom_text(
      data = eq_plankton %>% filter(source == "LISST"),
      aes(x = x, y = y, label = label),
      color = src_cols["LISST"],
      parse = TRUE,
      hjust = 1,
      size = 2.4,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = eq_plankton %>% filter(source == "ISIIS_plankton"),
      aes(x = x, y = y, label = label),
      color = src_cols["ISIIS_plankton"],
      parse = TRUE,
      hjust = 0,
      size = 2.4,
      fontface = "bold",
      inherit.aes = FALSE
    )
}

# ------------------------------------------------------------
# Figure 3: median-only comparison with all three sources
# ------------------------------------------------------------
p_medians_all <- ggplot() +
  
  # fitted median equation: allParticles from LISST-ISIIS merged coeffs
  geom_line(
    data = fit_medonly_all %>% dplyr::filter(fit_source == "allParticles_fit"),
    aes(x = size_num, y = NB_fit),
    color = "black",
    linewidth = 0.9,
    linetype = "dashed",
    alpha = 0.95
  ) +
  
  # fitted median equation: plankton
  geom_line(
    data = fit_medonly_all %>% dplyr::filter(fit_source == "plankton_fit"),
    aes(x = size_num, y = NB_fit),
    color = src_cols["ISIIS_plankton"],
    linewidth = 0.9,
    linetype = "dashed",
    alpha = 0.95
  ) +
  
  # existing median curves for the three sources
  geom_line(
    data = med_all,
    aes(x = size_num, y = NB_med, color = source),
    linewidth = 1.0,
    linetype = "solid",
    alpha = 0.95
  ) +
  geom_point(
    data = med_all,
    aes(x = size_num, y = NB_med, color = source),
    size = 0.65,
    alpha = 0.95
  ) +
  
  # equation text: allParticles top right
  geom_text(
    data = eq_medonly_all %>% dplyr::filter(fit_source == "allParticles_fit"),
    aes(x = x, y = y, label = label),
    color = "black",
    parse = TRUE,
    hjust = 1,
    vjust = 1,
    size = 4.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  # equation text: plankton just below
  geom_text(
    data = eq_medonly_all %>% dplyr::filter(fit_source == "plankton_fit"),
    aes(x = x, y = y, label = label),
    color = src_cols["ISIIS_plankton"],
    parse = TRUE,
    hjust = 1,
    vjust = 1,
    size = 4.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = panel_letter_df_medians,
    aes(x = x, y = y, label = panel_letter),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_color_manual(
    values = src_cols,
    name = "Source",
    labels = c("LISST", "ISIIS all particles", "ISIIS plankton")
  ) +
  scale_x_log10() +
  scale_y_log10(limits = c(1e-2, 1e10)) +
  facet_grid(
    Transect_ID ~ Cruise,
    drop = FALSE,
    labeller = labeller(
      Cruise = as_labeller(cruise_label_map),
      Transect_ID = as_labeller(transect_label_map)
    )
  ) +
  labs(
    x = expression(Equivalent~spherical~diameter~(mu*m)),
    y = expression(Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    
    legend.text = element_text(size = 14),
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(0.9, "cm"),
    
    legend.spacing.x = unit(0.8, "cm"),
    legend.box = "horizontal",
    
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 2.2,
        size = 3,
        alpha = 1
      )
    )
  )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave(
  filename = file.path(
    output_dir,
    paste0("multi_panel_LISST_vs_ISIIS_allParticles", if (opt$withEq) "_withEq" else "", ".png")
  ),
  plot = p_particles,
  width = 14,
  height = 16,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    paste0("multi_panel_LISST_vs_ISIIS_plankton", if (opt$withEq) "_withEq" else "", ".png")
  ),
  plot = p_plankton,
  width = 14,
  height = 16,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "multi_panel_medians_LISST_ISIIS_allParticles_ISIIS_plankton.png"
  ),
  plot = p_medians_all,
  width = 14,
  height = 16,
  units = "in",
  dpi = 300,
  bg = "white"
)

