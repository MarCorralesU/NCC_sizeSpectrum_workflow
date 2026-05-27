library(tidyverse)
library(data.table)
library(lubridate)
library(stringr)
library(fs)
library(optparse)
library(parsedate)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-c", "--seasonYear"), type = "character", help = "string that identifies a cruise [w22/s22/w23/s23]", metavar = "character"),
  make_option(c("-t", "--Transect_ID"), type = "character", help = "string that identifies a transect [GH/OB/CR/CM/NH/HH/RR]", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$seasonYear) || is.null(opt$Transect_ID)) {
  stop("Error: Both --seasonYear and --Transect_ID are required.")
}

#extract arguments
seasonYear <- opt$seasonYear
Transect_ID <- opt$Transect_ID


# === substring designation to search for files  ===
combined_dir <- "data/combined_data/"
substrings <- c(seasonYear, Transect_ID)  # dataset identifier substrings
camtop_tag <- "CamTop"
cambottom_tag <- "CamBottom"

# === Step 1: Filter files ===
all_files <- list.files(paste(base_dir, combined_dir, sep='/'), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

# Only files that match all substrings
filtered_files <- all_files %>% keep(~ all(str_detect(., substrings)))
camtop_files <- filtered_files[str_detect(filtered_files, camtop_tag)]
cambottom_files <- filtered_files[str_detect(filtered_files, cambottom_tag)]


valid_pattern <- paste0("^", seasonYear, "_", Transect_ID, ".*")

camtop_files <- camtop_files %>%
  keep(~ str_detect(basename(dirname(.)), valid_pattern))

cambottom_files <- cambottom_files %>%
  keep(~ str_detect(basename(dirname(.)), valid_pattern))


# === Step 2: Extract ALL timestamps from files ===
# Summarize time info per file
summarize_times <- function(file_path) {
  message("Reading timestamps from: ", basename(file_path))
  df <- fread(file_path, colClasses = list(character = "time"), showProgress = TRUE)
  if (any(grepl("^0022", df$time))) {
     df$time <- gsub("^0022", "2022", df$time)
  }
 #strip timezone offset, and set timezone
 df$time <- sub("([-+]\\d{2}:\\d{2})$", "", df$time)
 #retain time info and format
 df$time <- as.POSIXct(df$time, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Los_Angeles") 

  if (length(unique(df$time)) > 10) {
    message("Timestamps for this file (first 10 shown): ", 
            paste(head(unique(df$time), 10), collapse = ", "), " ...")
  } else {
    message("Timestamps for this file: ", paste(unique(df$time), collapse = ", "))
  }

  tibble(
    source_file = file_path,
    min_time = min(df$time, na.rm = TRUE),
    max_time = max(df$time, na.rm = TRUE),
    median_time = median(df$time, na.rm = TRUE)
  )
}

cat("summarizing CamTop timestamps...\n")
camtop_times <- map_dfr(camtop_files, summarize_times)

cat("Reading CamBottom timestamps...\n")
cambottom_times <- map_dfr(cambottom_files, summarize_times)



# === Step 3: Match timestamps (overlapping) ===
# We group all CamTop timestamps and match against all CamBottom timestamps
overlap_pairs <- expand_grid( camtop = camtop_times, cambottom = cambottom_times) %>%
  filter( camtop$min_time <= cambottom$max_time & camtop$max_time >= cambottom$min_time) %>%
  transmute(
    camtop_file = camtop$source_file,
    cambottom_file = cambottom$source_file,
    overlap_start = pmax(camtop$min_time, cambottom$min_time),
    overlap_end = pmin(camtop$max_time, cambottom$max_time)
  )

# Fix 'Camera_Top*' → 'CamTop' and 'Camera_Bottom*' → 'CamBottom'
clean_camera_labels <- function(df) {
  df %>%
    mutate(
      tar_id  = str_replace_all(tar_id, "Camera_Top", "CamTop"),
      tar_id  = str_replace_all(tar_id, "Camera_Bottom", "CamBottom"),
      image   = str_replace_all(image, "Camera_Top", "CamTop"),
      image   = str_replace_all(image, "Camera_Bottom", "CamBottom"),
      img   = str_replace_all(img, "Camera_Top", "CamTop"),
      img   = str_replace_all(img, "Camera_Bottom", "CamBottom")
    )
}

join_overlap_pair <- function(camtop_file, cambottom_file, overlap_start, overlap_end) {
  message("🔄 Processing:\n  CamTop:    ", basename(camtop_file), 
          "\n  CamBottom: ", basename(cambottom_file), 
          "\n  Range:     ", overlap_start, " → ", overlap_end)

  # Read and filter CamTop
  camtop <- fread(camtop_file, showProgress = FALSE, colClasses = list(character = "time"))
  if (any(grepl("^0022", camtop$time))) {
     camtop$time <- gsub("^0022", "2022", camtop$time)
  }
  #strip timezone offset, and set timezone
  camtop$time <- sub("([-+]\\d{2}:\\d{2})$", "", camtop$time)
  #retain time info and format
  camtop$time <- as.POSIXct(camtop$time, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Los_Angeles") 

  camtop <- camtop %>% filter(time >= overlap_start & time <= overlap_end)

  # Read and filter CamBottom
  cambottom <- fread(cambottom_file, showProgress = FALSE, colClasses = list(character = "time"))
  if (any(grepl("^0022", cambottom$time))) {
     cambottom$time <- gsub("^0022", "2022", cambottom$time)
  }
  #strip timezone offset, and set timezone
  cambottom$time <- sub("([-+]\\d{2}:\\d{2})$", "", cambottom$time)
  #retain time info and format
  cambottom$time <- as.POSIXct(cambottom$time, format = "%Y-%m-%d %H:%M:%OS", tz = "America/Los_Angeles") 
  
  cambottom <- cambottom %>% filter(time >= overlap_start & time <= overlap_end)

  # Combine and clean
  combined <- bind_rows(camtop, cambottom) %>%
    arrange(time) %>%
    clean_camera_labels()

  return(combined)
}


# Output directory
output_dir <- file.path(base_dir, "data", "combined_data_2Cams", seasonYear, Transect_ID)
if (!dir_exists(output_dir)) {
  dir_create(output_dir, recurse = TRUE)
  message("📂 Created directory: ", output_dir)
}

# Time formatting helper
#format_timestamp <- function(dt) {
  #format(dt, "%Y%m%dT%H%M%S", tz = "UTC")
#}

# Process row-by-row
for (i in seq_len(nrow(overlap_pairs))) {
  row <- overlap_pairs[i, ]

  result <- join_overlap_pair(
    camtop_file = row$camtop_file,
    cambottom_file = row$cambottom_file,
    overlap_start = row$overlap_start,
    overlap_end = row$overlap_end
  )

  start_str <- gsub(" ", "_", row$overlap_start) 
  end_str   <- gsub(" ", "_", row$overlap_end)

  outfile <- file.path(
    output_dir,
    paste0(seasonYear, "_", Transect_ID, "_", "joined_", start_str, "_to_", end_str, ".csv")
  )

  fwrite(result, outfile)
  message("✅ Wrote: ", outfile)
}









