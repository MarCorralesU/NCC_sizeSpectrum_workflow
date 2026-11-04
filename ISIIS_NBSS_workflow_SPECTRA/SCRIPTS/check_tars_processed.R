#!/jet/home/corrales/R
library(dplyr)

base_dir = getwd()

files_tars<- gsub('.tar.gz', '', basename(readLines(paste(base_dir, 'ancillary/tar_lists/inputTars_7.txt', sep='/'))))

files_csv <- gsub('.csv', '', basename(list.files(path=paste(base_dir, 'CSV', sep='/'),full.names= TRUE, pattern = "*.csv", recursive=TRUE)))

no_csv<- setdiff(files_tars, files_csv)

print(no_csv)
