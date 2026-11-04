library(dplyr)
library(magrittr)
library(tidyr)
library(NLP)
library(here)
library(hms)
library(lubridate)
library(tidyverse)
library(scales)
library(fs)
library(stringr)
library(optparse)
library(vroom)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-d", "--combinedDir"), type = "character", help = "Path to combined ISIIS data of the merged plus size data ", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$combinedDir)) {
  stop("Error: --combinedDir is required.")
}

combined_dir <- opt$combinedDir #use this for testing combined_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/combined_data_2Cams/s23/NH"


df_list <- list.files(file.path(combined_dir), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)
# Extract Transect (last folder)
Transect <- basename(combined_dir)
# Extract Cruise (the folder above)
Cruise <- basename(dirname(combined_dir))

n_frames_avi <- data.frame()

for (d in df_list){
    n_frames_tar_id <- read_csv(d) %>%
                            mutate(frame_id = str_remove(img, "_crop_.*$")) %>% 
                            group_by(tar_id, frame_id) %>%
                            summarise(mean_Latitude = mean(Latitude, na.rm = TRUE),
                                      mean_Longitude = mean(Longitude, na.rm = TRUE), .groups = "drop")

    n_frames_avi <- rbind(n_frames_avi, n_frames_tar_id)

}

#routine to get frame count per avi

n_frames_avi<- n_frames_avi %>%
                   group_by(tar_id) %>%
                   summarise(n_frames = n_distinct(frame_id),
                             mean_Latitude = mean(mean_Latitude, na.rm = TRUE),
                             mean_Longitude = mean(mean_Longitude, na.rm = TRUE), .groups = "drop")
frame_n_avi_path <- file.path(base_dir, paste(Cruise, Transect, 'frame_count_per_avi.csv', sep='_'))
write_csv(n_frames_avi, frame_n_avi_path)
