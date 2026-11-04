
library(tidyverse)
library(scales)
library(fs)
library(stringr)
library(chunked)
library(optparse)
library(LaF)

# 🧠 Limit multithreading to avoid crashes/memory bloat


base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-m", "--mergedDir"), type = "character", help = "Path to ISIIS merged CSV", metavar = "character"),
  make_option(c("-l", "--listFile"), type = "character", default = "NO", help = "Generate .txt file with CSV paths [YES/NO]", metavar = "character"),
  make_option(c("-c", "--CSVpath"), type = "character", default = "CSV", help = "path to the CSVs inside the project", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$mergedDir)) {
  stop("Error: --mergedDir is required.")
}

merged_dir <- opt$mergedDir # #use this for testing merged_dir<- "/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/data/merged_files/105_MBON_merge.csv" 
file_name <- tools::file_path_sans_ext(basename(merged_dir))
message("merged file to be combined: ", toString(file_name))
make_csv_text <- toupper(opt$listFile)  # Make it case-insensitive

if (make_csv_text == "YES") {
  csv_dir <- file.path(base_dir, opt$CSVpath)
  output_file <- file.path(base_dir, "ancillary", "size_csv_paths_MBON.txt")
  
  # Make sure the output directory exists
  dir_create(dirname(output_file))
  
  # Bash one-liner: find all .csv files and write full paths
  cmd <- paste0("find ", shQuote(csv_dir), " -type f -name '*.csv' > ", shQuote(output_file))
  system(cmd)
  message("✅ Wrote size_csv_paths_MBON.txt and exiting as requested.")
  quit(save = "no")
}

#old way of doing it with LaF, but produces memory crash at C level
#get column info
# load data


column_names <- names(read_csv(merged_dir, n_max = 0))
#load list of size files
size_csv_paths <- readLines(file.path(base_dir, "ancillary", "size_csv_paths_MBON.txt"))

chunk_size <- 1e6

for (i in 1:1000){
    message("\n--- Iteration ", i, " ---")
    skip_rows <- chunk_size * (i - 1)
    cols_to_drop <- c("HD_ID", "Vol.imaged", "Altitude", "Pitch", "Roll", "Heading", "Trans.Velocity", "Long..Velocity", "Normal.Velocity", "Wing.Fin.Angle", "Roll.Fin.Angle", "Transect_ID")
    result <- try({
        binned <- read_csv(  merged_dir, skip = skip_rows, n_max = chunk_size, col_names = column_names, col_types = cols(.default = "c"))
        # If first iteration, discard the first row (header)
        if (i == 1) {
            binned <- binned[-1, , drop = FALSE]
        }
        if (nrow(binned) == 0) {
            message(" ⚠️  No more data after chunk ", i - 1, ". Exiting loop.")
            break}
        binned <- binned %>% select(-any_of(cols_to_drop))
        binned$tar_id <- str_extract(binned$img, "^[^_]+")

        #Get unique tar_ids
        tar_ids <- unique(binned$tar_id)
        message("Tar ids in this chunk: ", toString(tar_ids))
        if (length(tar_ids) == 0) {
            message(" ⚠️  No tar_ids found in this chunk, skipping iteration.")
            next}
        #Filter paths
        filtered_paths <- size_csv_paths[str_detect(size_csv_paths, paste(tar_ids, collapse = "|"))]
        if (length(filtered_paths) == 0) {
            message(" ⚠️  No matching .csv files found in size_csv_paths for the current binned$tar_id values, skipping iteration.")
            next}
        message("✅ Matching paths found. Proceeding to read .csv files...")
        message("Corresponding paths are:\n", paste(filtered_paths, collapse = "\n"))
 
        #Read all CSVs and combine
        size_data <- bind_rows(lapply(filtered_paths, function(path) {read_csv(path, col_types = cols(.default = "c"))} )) 
        binned <- left_join(binned, size_data, by = c("img" = "crop")) %>% filter(!is.na(contour))
        
        # Create output directory if it doesn't exist
        output_dir <- file.path(base_dir, "data", "combined_data", file_name)
        if (!dir_exists(output_dir)) {
            dir_create(output_dir)
            message("📁 Created directory: ", output_dir)
        }

        #save in chunks
        output_path_binned <- file.path(output_dir, paste0(file_name, "-combined_data_", i, ".csv"))
        write_csv(binned, output_path_binned)
        message("✅ Saved to: ", output_path_binned)

        #delete to avoid memory overload
        rm(binned)
        rm(size_data)
    }, silent=TRUE)
    if (inherits(result, "try-error")) {
        message("⚠️ [Iteration ", i, "] failed due to internal error.")
    next
    }

}


#old stuff when i tried to create one file per merged

