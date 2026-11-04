#!/jet/home/corrales/R
library(dplyr)
library(stringr)

base_dir = getwd()
args <- commandArgs(trailingOnly = TRUE)
print(args)
#stop the script if no command line argument
if(length(args)==0){
  print("this function requires three arguments:")
  print("1) one or more directories where to look for tar files!")
  print("2) whether you need to generate the split tarfile list for computing (split/noSplit)")
  print("3) the number of tar files that will be listed in the split text files")
  stop("Requires command line argument.")
}
# get arguments
dirs<- grep('/', args, value=TRUE)
split_YN<- head(tail(args,2),1)
files_per_split<-as.numeric(tail(args,1))

#extract tar paths
all_tars<-list()
for (d in dirs){
    files <- list.files(path= d,full.names= TRUE, pattern = "*.tar.gz", recursive=TRUE)
    all_tars <- c(all_tars, files)
}
all_tars <- Filter(function(x) !any(grepl("testdata", x)), all_tars)

# create dataframe for tracking processing progress and to sort the tar paths in processing priority
columns<-c('filePath', 'Cruise', 'year', 'season', 'transect')
df<-data.frame(matrix(nrow = length(all_tars), ncol = length(columns)))
colnames(df) = columns
df$filePath<-all_tars

df<-df %>%
  mutate(Cruise = case_when(
    grepl(pattern = "w22", x = filePath) ~ "2022_winter",
    grepl(pattern = "s22", x = filePath) ~ "2022_summer",
    grepl(pattern = "w23", x = filePath) ~ "2023_winter",
    grepl(pattern = "s23", x = filePath) ~ "2023_summer"
    
  )) %>%
mutate(transect = case_when(
  grepl(pattern = "GH", x = filePath) ~ "GraysHarbor",
  grepl(pattern = "CR", x = filePath) ~ "ColumbiaRiver",
  grepl(pattern = "CM", x = filePath) ~ "CapeMeares",
  grepl(pattern = "NH", x = filePath) ~ "NewportLine",
  grepl(pattern = "HH", x = filePath) ~ "HecetaHead",
  grepl(pattern = "RR", x = filePath) ~ "RogueRiver",
  grepl(pattern = "OB", x = filePath) ~ "OceanBeaches",
  
))



df$Cruise <- factor(df$Cruise, levels = c("2022_winter", "2022_summer", "2023_winter", "2023_summer"))
df$transect <- factor(df$transect, levels = c("NewportLine","GraysHarbor","RogueRiver","ColumbiaRiver","HecetaHead","CapeMeares","OceanBeaches"))

df<-df[order(df$transect,df$Cruise),]

df[c('year', 'season')] <- str_split_fixed(df$Cruise, '_', 2)
df <- apply(df,2,as.character)

write.csv(df, paste(base_dir, 'ancillary/tar_file_list_ordered.csv', sep='/'),row.names=FALSE)


# extract the ordered tar path list and save as text
all_tars<-as.list(df[,'filePath'])
fileTars<-file(paste(base_dir, "inputTars.txt", sep='/'))
writeLines(as.character(all_tars), fileTars)
close(fileTars)

if (split_YN=='split'){
splitIndices = seq(1,length(all_tars),files_per_split)

#for each index, split dataset and export
lapply(splitIndices,function(x)  { 
    minIndex = x; 
    maxIndex = min(x+files_per_split-1,length(all_tars)); 
    DF = all_tars[minIndex:maxIndex]
    subset_tar<- file(paste0(base_dir,"/ancillary/tar_lists/inputTars_",1+as.integer(x/files_per_split),".txt")) 
    writeLines(as.character(DF), subset_tar)
})
}
