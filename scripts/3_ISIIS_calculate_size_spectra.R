
library(tidyverse)
library(scales)
library(fs)
library(stringr)
library(chunked)
library(optparse)
library(data.table)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-b", "--binnedDir"), type = "character", help = "Path to binned data ", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-i", "--ISIISver"), type = "integer", default = 2, help = "ISIIS version, either 2 or 3", metavar = "integer"),
  make_option(c("-g", "--groupingVars"), type = "character", default = "Transect_name,shelf_loc_200,season,depth_mid,lat_mid,lon_mid", help = "Comma-separated grouping variables", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temp,Salinity,Fluoro,Oxygen,Irrandiance", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "character", default = "TRUE", help = "whether to threshold the size spectra TRUE/FALSE", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$binnedDir)) {
  stop("Error: binnedDir is required.")
}

#extract arguments
binned_dir <- opt$binnedDir # #use this for testing binned_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_MBON/data/binned_data"
size_var<- opt$sizeVar #size_var<-'contour'
ISIIS_ver <- opt$ISIISver#ISIIS_ver<-2
grouping_vars <- strsplit(opt$groupingVars, ",")[[1]] #grouping_vars <- strsplit(c("Transect_name,shelf_loc_200,season,depth_mid,lat_mid,lon_mid"), ",")[[1]]
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit(c("Temp,Salinity,Fluoro,Oxygen,Irrandiance"), ",")[[1]]
threshold <- tolower(opt$threshold) %in% c("true", "t", "1")

#create binned data if it doesn't exist
if (!dir.exists(file.path(base_dir, "data", "size_spectra"))) {
    dir.create(file.path(base_dir, "data", "size_spectra"), recursive = TRUE)
}

#get data grouping identifiers
df_list <- list.files(file.path(binned_dir,size_var), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)
pattern <- paste0("merge_\\d+_(.*?)_", size_var, "_binned\\.csv")
group_id <- str_match(basename(df_list), pattern)[,2]


df_all<- data.frame()

if (threshold==TRUE){
    df_noThres_all <- data.frame()
}

df_all_taxa<- data.frame()


for (g in unique(group_id)){
    message("⏳ Processing group: ", g)
   list_group <- df_list[grep(g, df_list)]
   #print(list_group)
   df_raw <- list_group %>% lapply(fread,  data.table=FALSE) %>% bind_rows %>% as.data.frame()

   if (threshold==FALSE){
       df_noThres<- df_raw %>% 
           group_abundance_func(., grouping_vars = grouping_vars,
                                env_vars = env_vars,
                                group_by_taxa =FALSE, ISIIS_ver=ISIIS_ver, thresholding = FALSE)
   }

   df_no_taxa<- df_raw %>% 
        group_abundance_func(., grouping_vars = grouping_vars,
                             env_vars = env_vars,
                             group_by_taxa =FALSE, ISIIS_ver=ISIIS_ver) %>%
        lin_fit_func(., filter_min_size_class=TRUE, n_min_size_class=3)

   df_taxa <- df_raw %>% 
              group_abundance_func(., grouping_vars = grouping_vars,
                                   env_vars = env_vars, 
                                   group_by_taxa =TRUE, ISIIS_ver=ISIIS_ver) %>% 
              lin_fit_func(., filter_min_size_class=TRUE, n_min_size_class=3)
  # Combine results
  df_all <- if (exists("df_all")) bind_rows(df_all, df_no_taxa) else df_no_taxa
  if (threshold==FALSE){
      df_noThres_all <- if (exists("df_noThres_all")) bind_rows(df_noThres_all, df_noThres) else df_noThres
  }
  df_all_taxa <- if (exists("df_all_taxa")) bind_rows(df_all_taxa, df_taxa) else df_taxa

  # Free memory
  if (threshold==FALSE){
      rm(df_raw, df_noThres, df_no_taxa, df_taxa)
  }else{
      rm(df_raw, df_no_taxa, df_taxa)
  }
  gc(verbose = FALSE)
}

#assign group numbers
df_all<- df_all %>% group_id_func(., grouping_vars=c(grouping_vars, env_vars), group_by_taxa =FALSE)# this function is necesary because of how the groups are being separated in the binning step
if (threshold==FALSE){
    df_noThres_all <- df_noThres_all  %>% group_id_func(., grouping_vars=c(grouping_vars, env_vars), group_by_taxa =FALSE)
}
df_all_taxa<- df_all_taxa %>% group_id_func(., grouping_vars=c(grouping_vars, env_vars), group_by_taxa =TRUE)


write.csv(df_all, file.path(base_dir, "data", "size_spectra",'ISIIS_NBSS_by_community.csv'), row.names=FALSE)

if (threshold==FALSE){
    write.csv(df_noThres_all, file.path(base_dir, "data", "size_spectra",'ISIIS_NBSS_by_community_no_threshold.csv'), row.names=FALSE)
}
write.csv(df_all_taxa, file.path(base_dir, "data", "size_spectra", 'ISIIS_NBSS_by_taxa.csv'), row.names=FALSE)

