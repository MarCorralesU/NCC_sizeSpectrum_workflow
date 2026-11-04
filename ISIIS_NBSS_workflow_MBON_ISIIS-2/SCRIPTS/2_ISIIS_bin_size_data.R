
#Goal: first step toward estimating size spectra from ISIIS, This will be based on Mo's and Elena's previous scripts
library(data.table)
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
library(LaF)
library(optparse)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-c", "--combinedDir"), type = "character", help = "Path to combined ISIIS data of the merged plus size data ", metavar = "character"),
  make_option(c("-i", "--ISIISver"), type = "integer", default = 2, help = "ISIIS version, either 2 or 3", metavar = "integer"),
  make_option(c("-v", "--dataVersion"), type = "character", default = "MEZCAL", help = "the version of the data, either MEZCAL (use this for MBON) or SPECTRA", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-t", "--transectLats"), type = "character", default = "ancillary/transect_lats.csv", help = "Path to CSV with reference transect latitudes", metavar = "character"),
  make_option(c("-b", "--shelfBreak"), type = "character", default = "ancillary/shelf_break_200.csv", help = "Path to CSV with shelf break longitudes per transect", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$combinedDir)) {
  stop("Error: --combinedDir is required.")
}

#extract arguments
combined_dir <- opt$combinedDir # #use this for testing combined_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/data/combined_data/101_MBON_merge/101_MBON_merge-combined_data_1.csv" 
file_name <- gsub("-combined_data", "", tools::file_path_sans_ext(basename(combined_dir)))
ISIISver<- opt$ISIISver
project<- toupper(opt$dataVersion)  # Make it case-insensitive
size_var<- opt$sizeVar
transect_lat_path <- file.path(base_dir, opt$transectLats)
shelf_break_path <- file.path(base_dir, opt$shelfBreak)



#generate binned files
binned<- read_csv(combined_dir)  %>% filter_low_class_scores(., project=project) %>% convert_latlong_dms() %>% assign_transect_name(transect_lat_path) %>% assign_shelf_location(shelf_break_path) %>% filter(!grepl("artifact", predicted_taxon)) %>% assign_season(.)
binned<-binned %>% depth_bin_calc(., depth_bin_width=1) 
binned<-binned %>% coord_bin_calc(.,lat_var='Lat',lon_var='Long', coord_bin_width=0.25) 
if (nrow(binned) == 0) {
  message("🛑 coord_bin_calc() returned an empty dataframe (likely all depths <3 m). Stopping process.")
  quit(save = "no", status = 1, runLast = FALSE)
}
binned<-binned %>% size_binning(., size_var=size_var, depth_bin_width=5, ISIIS_ver=ISIISver)
binned<-binned %>% select(c(Transect_name, shelf_loc_200,  time, season, img, Temp, Salinity, Fluoro, Oxygen, Irrandiance,  predicted_taxon, depth_mid, lat_mid, lon_mid, biovol, biovol_binned, biovol_um3_mid, biovol_um3_width)) #%>%

binned$grouping_col<- paste(as.character(binned$Transect_name), as.character(binned$shelf_loc_200),  as.character(binned$season), as.character(binned$depth_mid), as.character(binned$lat_mid), as.character(binned$lon_mid), sep='_')

output_dir <- file.path(base_dir, "data", "binned_data", size_var)
if (!dir_exists(output_dir)) {
  dir_create(output_dir, recurse = TRUE)
  message("📁 Created directory: ", output_dir)
}
  
for (g in unique(binned$grouping_col)){
     df_subset<-binned[binned$grouping_col==g,]
     output_path<- file.path(output_dir, paste(file_name, g, size_var,'binned.csv', sep='_'))
     write_csv(df_subset, output_path)
     message("✅ Saved to: ", output_path)}  
