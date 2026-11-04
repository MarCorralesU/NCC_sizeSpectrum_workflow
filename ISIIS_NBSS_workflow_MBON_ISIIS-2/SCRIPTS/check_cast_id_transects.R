#!/usr/bin/env Rscript
setwd('/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_MBON_ISIIS-2')
library(optparse)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)

base_dir<-getwd()

option_list <- list(
  make_option(c("-c","--cruise"), type = "character"),
  make_option(c("-t","--transect"), type = "character"),
  make_option(c("-o","--tow"), type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))


binned_file <- file.path(
  base_dir,
  "data",
  "binned_abundance",
  paste0(opt$transect, "_", opt$cruise, "_", opt$tow, "_binned_abundance_chunk_1.csv")
)

size_dir <- file.path(
  base_dir,
  "data",
  "size_spectra_plankton_5-10m-bins_thresholded",
  "contour",
  opt$cruise,
  opt$transect
)

size_files <- list.files(
  size_dir,
  pattern = paste0("^", opt$cruise, "-", opt$transect, "-", opt$tow,
                   "_cast_[0-9]+-ISIIS_NBSS_by_community\\.csv$"),
  full.names = TRUE
)

message("Reading binned abundance:")
message(binned_file)

dat_binned <- read_csv(binned_file, show_col_types = FALSE) %>%
  transmute(
    Longitude = mean_Longitude,
    depth_mid = depth_mid,
    time = mean_time_Pacific,
    source = "Binned abundance"
  )

message("Reading size spectra files:")
message(length(size_files), " files found")

dat_size <- map_dfr(size_files, ~ read_csv(.x, show_col_types = FALSE)) %>%
  transmute(
    Longitude = Longitude,
    depth_mid = depth_mid,
    time = time_Pacific,
    source = "Size spectra"
  )

dat <- bind_rows(dat_binned, dat_size) %>%
  filter(!is.na(Longitude), !is.na(depth_mid))

dat <- dat %>%
  arrange(source, time)

p <- ggplot() +
  
  ## Binned abundance: thick, transparent
  geom_path(
    data = filter(dat, source == "Binned abundance"),
    aes(x = Longitude, y = depth_mid),
    colour = "#D55E00",
    linewidth = 2,
    alpha = 0.35
  ) +
  
  ## Size spectra: thin, solid
  geom_path(
    data = filter(dat, source == "Size spectra"),
    aes(x = Longitude, y = depth_mid),
    colour = "#0072B2",
    linewidth = 0.7
  ) +
  
  scale_y_reverse() +
  labs(
    x = "Longitude",
    y = "Depth midpoint (m)"
  ) +
  theme_bw()

out_dir <- file.path(base_dir, "data", "towYo_diagnostic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(
  out_dir,
  paste0(opt$cruise, "_", opt$transect, "_", opt$tow,  "_towYo_diagnostic.png")
)

ggsave(outfile, p, width = 12, height = 8, dpi = 300)

message("Saved to:")
message(outfile)


ggsave(outfile, p, width = 10, height = 8, dpi = 300)

message("Saved to:")
message(outfile)
