#script to generate size classes based on minimum sizes from training library
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
base_dir<- getwd()

# --- 0) Paths & identifiers (ADAPT: uses your base_dir and w23_NH) ----
# main directory for outputs (pick your SPECTRA workflow root)

source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


pixel_size<- 55 # for ISIIS-3
# size classes will be built on contour:

df_train <- suppressMessages(
  readr::read_csv(file.path(base_dir, "ancillary", "training_lib_sizes.csv"),
                  show_col_types = FALSE,
                  col_select = -1))

df_train$biovol<-(4/3)* pi * ((((df_train$contour*(pixel_size**2))/pi)**0.5) **3)

df_train$ECD<-2*sqrt((df_train$contour*(pixel_size**2))/pi)

cols <- c("biovol", "ECD", "contour")


df_ranges <- df_train |>
  dplyr::select(dplyr::all_of(cols)) |>
  tidyr::pivot_longer(dplyr::everything(), names_to = "metric", values_to = "value") |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(factor(metric, levels = cols))

df_ranges

# size bin construction from Kiko et al. (2022) and Dugenne et al. (2024): logarithmically spaced
#using a base of 2 and an increment of 1=3 so that a doubling in equivalent circular diameter (ECD) is observed in
#every third bin (equivalent to a doubling in biovolume observed in every bin), with a range between minimum size on training library to 50 000 000 µm

# Inputs from training library minimum sizes
ECD_min <- df_ranges %>% dplyr::filter(metric == "ECD") %>% dplyr::pull(min) %>% as.numeric()
ECD_max_target <- 5e7                          # 50,000,000 µm (blue whale size)
ratio <- 2^(1/3)                               # 3 bins per doubling in ECD

make_plankton_bins <- function(ECD_min, ECD_max = 5e7, ratio = 2^(1/3)) {
  stopifnot(is.finite(ECD_min), ECD_min > 0, ECD_max > ECD_min, ratio > 1)
  
  # Number of steps needed and edges in ECD
  k <- ceiling(log(ECD_max / ECD_min, base = ratio))
  ecd_edges <- ECD_min * ratio^(0:k)
  if (tail(ecd_edges, 1) < ECD_max) ecd_edges <- c(ecd_edges, tail(ecd_edges, 1) * ratio)
  
  # Convert edges to biovolume (µm^3): V = (pi/6) * ECD^3
  bv_edges <- (pi / 6) * (ecd_edges^3)
  
  tibble(
    bin_id        = seq_len(length(ecd_edges) - 1),
    # ECD bins (µm)
    ECD_lower     = ecd_edges[-length(ecd_edges)],
    ECD_upper     = ecd_edges[-1],
    ECD_mid       = sqrt(ECD_lower * ECD_upper),
    ECD_width     = ECD_upper - ECD_lower,
    # Biovolume bins (µm^3)
    biovol_lower  = bv_edges[-length(bv_edges)],
    biovol_upper  = bv_edges[-1],
    biovol_um3_mid    = sqrt(biovol_lower * biovol_upper),
    biovol_um3_width  = biovol_upper - biovol_lower
  )
}

bins <- make_plankton_bins(ECD_min = ECD_min, ECD_max = ECD_max_target, ratio = ratio)
bins

write.csv(bins, paste(base_dir, 'ancillary', 'custom_ISIIS-3_size_classes.csv', sep='/'))

