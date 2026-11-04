library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("-c", "--CRUISE"), type = "character", help = "Cruise to get mean dateTime and location(Lat/Lon) ", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$CRUISE)) {
  stop("Error: --CRUISE is required.")
}

base_dir <- getwd()
input_dir <- file.path(base_dir, "data", "binned_abundance")
output_dir <- file.path(input_dir, "date_plus_location")

# Make sure output dir exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# CRUISE identifier (set this)
CRUISE <- opt$CRUISE   # <-- change to match your cruise string

# Get matching files
files <- list.files(
  input_dir, 
  pattern = paste0(CRUISE, ".*\\.csv$"), 
  full.names = TRUE
)

# Read, select cols, combine
df_all <- files %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE) %>%
             select(mean_time_Pacific, mean_Latitude, mean_Longitude))

# Output file path
out_file <- file.path(output_dir, paste(CRUISE, "dateTime_locations.csv", sep = "_"))

# Write combined csv
write_csv(df_all, out_file)

message("Saved combined file: ", out_file)
