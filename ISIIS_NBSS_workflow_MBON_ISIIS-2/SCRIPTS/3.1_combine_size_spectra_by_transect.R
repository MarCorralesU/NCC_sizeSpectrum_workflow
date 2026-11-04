#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)
library(purrr)

base_dir <- "/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON"
spectra_dir <- file.path(base_dir, "data", "size_spectra")
output_dir <- spectra_dir  # Output to same dir

# List all CSV files
csv_files <- list.files(spectra_dir, pattern = "\\.csv$", full.names = TRUE)

# Extract grouping keys
meta_df <- tibble(
  file = csv_files,
  filename = basename(csv_files)
) %>%
  mutate(
    transect_name = str_extract(filename, "^[^_]+"),
    file_type = case_when(
      str_detect(filename, "by_taxa") ~ "by_taxa",
      str_detect(filename, "by_community_no_threshold") ~ "by_community_no_threshold",
      str_detect(filename, "by_community") ~ "by_community",
      TRUE ~ "unknown"
    )
  ) %>%
  filter(file_type != "unknown")

# Combine per (transect_name, file_type)
meta_df %>%
  group_by(transect_name, file_type) %>%
  group_split() %>%
  walk(function(group_df) {
    combined <- group_df$file %>%
      map_dfr(read_csv, show_col_types = FALSE)
    
    out_name <- paste0("COMBINED_", group_df$transect_name[1], "_", group_df$file_type[1], ".csv")
    out_path <- file.path(output_dir, out_name)
    
    write_csv(combined, out_path)
    message("✅ Wrote: ", out_path)
  })

message("🎉 All done combining size spectra.")

