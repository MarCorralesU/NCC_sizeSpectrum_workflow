library(dplyr)
library(stringr)
library(optparse)

# Define options
option_list <- list(
  make_option(c("-d", "--dirs"), type = "character", action = "store", help = "Comma-separated list of directories where to look for tar files"),
  make_option(c("-s", "--split"), type = "character", default = "noSplit", help = "Whether to split tarfile list (split/noSplit)"),
  make_option(c("-n", "--nfiles"), type = "integer", default = 0, help = "Number of files per split (0 means no split)"),
  make_option(c("-p", "--pipeline"), type = "character", default = "SPECTRA", help = "Pipeline name (e.g., SPECTRA)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

base_dir <- getwd()
dirs <- str_split(opt$dirs, pattern = ",")[[1]]
split_YN <- opt$split
files_per_split <- opt$nfiles
pipeline <- opt$pipeline

# Extract tar paths
all_tars <- list()
for (d in dirs) {
  files <- list.files(path = d, full.names = TRUE, pattern = "*.tar.gz", recursive = TRUE)
  all_tars <- c(all_tars, files)
}
all_tars <- Filter(function(x) !any(grepl("testdata", x)), all_tars)

if (pipeline == "SPECTRA") {
  columns <- c('filePath', 'Cruise', 'year', 'season', 'transect')
  df <- data.frame(matrix(nrow = length(all_tars), ncol = length(columns)))
  colnames(df) <- columns
  df$filePath <- all_tars

  df <- df %>% 
    mutate(Cruise = case_when(
      grepl("w22", filePath) ~ "2022_winter",
      grepl("s22", filePath) ~ "2022_summer",
      grepl("w23", filePath) ~ "2023_winter",
      grepl("s23", filePath) ~ "2023_summer"
    )) %>% 
    mutate(transect = case_when(
      grepl("GH", filePath) ~ "GraysHarbor",
      grepl("CR", filePath) ~ "ColumbiaRiver",
      grepl("CM", filePath) ~ "CapeMeares",
      grepl("NH", filePath) ~ "NewportLine",
      grepl("HH", filePath) ~ "HecetaHead",
      grepl("RR", filePath) ~ "RogueRiver",
      grepl("OB", filePath) ~ "OceanBeaches"
    ))

  df$Cruise <- factor(df$Cruise, levels = c("2022_winter", "2022_summer", "2023_winter", "2023_summer"))
  df$transect <- factor(df$transect, levels = c("NewportLine", "GraysHarbor", "RogueRiver", "ColumbiaRiver", "HecetaHead", "CapeMeares", "OceanBeaches"))
  df <- df[order(df$transect, df$Cruise), ]
  df[c('year', 'season')] <- str_split_fixed(df$Cruise, '_', 2)
  df <- apply(df, 2, as.character)
  write.csv(df, file.path(base_dir, 'ancillary/tar_file_list_ordered.csv'), row.names = FALSE)
}

# Create inputTars.txt
fileTars <- file(file.path(base_dir, "inputTars.txt"))
writeLines(all_tars, fileTars)
close(fileTars)

# Create symlinks in input_tars dir
symlink_dir <- file.path(base_dir, "input_tars")
if (!dir.exists(symlink_dir)) dir.create(symlink_dir)

unique_filenames <- basename(all_tars)
for (i in seq_along(all_tars)) {
  target <- all_tars[i]
  link_name <- file.path(symlink_dir, unique_filenames[i])
  if (!file.exists(link_name)) file.symlink(target, link_name)
}

# Optionally split inputTars.txt
if (split_YN == 'split' && files_per_split > 0) {
  splitIndices <- seq(1, length(all_tars), files_per_split)
  for (i in seq_along(splitIndices)) {
    minIndex <- splitIndices[i]
    maxIndex <- min(minIndex + files_per_split - 1, length(all_tars))
    DF <- all_tars[minIndex:maxIndex]
    out_path <- file.path(base_dir, sprintf("ancillary/tar_lists/inputTars_%d.txt", i))
    writeLines(as.character(DF), file(out_path))
  }
}
