setwd('/Users/user/SPECTRA_proj/')
library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)
library(NLP)
library(here)

base_dir<-getwd()

# assign season label
assign_season <- function(df_clean) {
  # derive helpers
  MONTH <- lubridate::month(df_clean$dateTime)
  DAY   <- lubridate::mday(df_clean$dateTime)
  DPM   <- DAY + MONTH * 100L
  
  # base R season classification
  season <- ifelse(DPM >= 320  & DPM <= 620,  "Spring",
            ifelse(DPM >= 621  & DPM <= 922,  "Summer",
            ifelse(DPM >= 923  & DPM <= 1221, "Fall", "Winter")))
  
  df_clean$season <- season
  df_clean
}

format_LISST_files_processed<- function(filename) {
  # load the middle of size classes.
  mid_size_class <- as.data.frame(transpose(list(read.csv(file.path(base_dir, "ancillary", "LISST_size_classes_um.csv"), sep = ',')$size_class_mid_um)))
  df_bin_width<-read.csv(file.path(base_dir,"ancillary", "LISST_size_classes_um.csv"), sep = ',')%>% mutate(size_class_mid_um = factor(size_class_mid_um))
  # refer to LISST documentation to know which columns refer to which data type. These are the target data:
  # 1:36: volume concentration of size classes 
  # 41: Depth (m of seawater)
  # 42: Temperature (C)
  # 43: Year
  # 44: Month
  # 45: Day
  # 46: hour
  # 47: minute
  # 48: second
  # separate data sources, modify and merge
  df <- read.csv(filename, sep = ',', header=F)
  # partricle concentrations:
  conc <- df[,1:36] %>% set_colnames(mid_size_class[1,])
  depth <- data.frame(df[,41]) %>% set_colnames('depth_m')
  temperature <- data.frame(df[,42]) %>% set_colnames('temperature_c')
  date <- apply(df[,43:45], 1:2, function(x) sprintf("%02s", x)) %>% 
    data.frame() %>%
    unite('date', colnames(df[,43:45])[1:3], sep='-')
  time <- apply(df[,46:48], 1:2, function(x) sprintf("%02s", x)) %>%   
    data.frame() %>%
    unite('time', colnames(df[,46:48])[1:3], sep=':') 
  dateTime <-cbind(date, time) %>% 
    unite( col='dateTime' ,c('date', 'time'), sep=' ') %>% 
    lapply( function(x) strptime(x, format = '%Y-%m-%d %H:%M:%S')) %>%
    as.data.frame()
  df_clean <-melt(as.data.table(cbind(dateTime, depth, temperature, conc)), id.vars=c('dateTime','depth_m', 'temperature_c'), variable.name='size_class_mid_um', value.name = 'concentration_uL_L') %>% 
    dplyr::arrange(dateTime) %>% assign_season()
  df_clean<-merge(df_clean, df_bin_width, by = "size_class_mid_um", all.x = TRUE)
  df_clean$size_class_mid_ul<-(((as.numeric(paste(df_clean$size_class_mid_um))/2)^3)*1.333*3.14)*1e-9
  df_clean$bin_width_ul<-(((df_clean$bin_width/2)^3)*1.333*3.14)*1e-9
  df_clean <- df_clean %>% filter(depth_m> 0) %>% filter(concentration_uL_L> 0)#[df_clean$depth_m > 0, ]
  
  rm(conc); rm(depth); rm(temperature); rm(date); rm(time); rm(dateTime); rm(mid_size_class); rm(df_bin_width)
  return(df_clean)
}

format_LISST_files_custom_processed<- function(filename) {
  # load the middle of size classes.
  mid_size_class <- as.data.frame(transpose(list(read.csv(file.path(base_dir, "ancillary", "LISST_size_classes_um.csv"), sep = ',')$size_class_mid_um)))
  df_bin_width<-read.csv(file.path(base_dir, "ancillary", "LISST_size_classes_um.csv"), sep = ',')%>% mutate(size_class_mid_um = factor(size_class_mid_um))
  # refer to LISST documentation on custom MATLAB processing and custom script process_LISST_RBN.m to know which columns refer to which data type. These are the target data:
  # 1:36: volume concentration of size classes 
  # 37: Depth (m of seawater)
  # 38: Temperature (C)
  # 39: Year
  # 40: Month
  # 41: Day
  # 42: hour
  # 43: minute
  # 44: second
  # separate data sources, modify and merge
  df <- read.csv(filename, sep = ',', header = TRUE, check.names = FALSE)
  # partricle concentrations:

  df_clean <-melt(df, id.vars=c('dateTime','Depth', 'Temperature'), variable.name='size_class_mid_um', value.name = 'concentration_uL_L', variable.factor = FALSE ) %>%  dplyr::arrange(dateTime)
  df_clean <- assign_season(df_clean)
  df_clean$concentration_um3_L<- df_clean$concentration_uL_L*1e9
  
  df_clean$size_class_mid_um <- as.numeric(sub("^X", "", df_clean$size_class_mid_um))
  df_clean$size_class_mid_um <- round(df_clean$size_class_mid_um, 6)
  
  df_bin_width$size_class_mid_um <- round(as.numeric(as.character(df_bin_width$size_class_mid_um)), 6)
  
  df_clean<-merge(df_clean, df_bin_width, by = "size_class_mid_um", all.x = TRUE)
  
  df_clean$size_class_mid_um3<-(((as.numeric(paste(df_clean$size_class_mid_um))/2)^3)*1.333*pi)
  
  df_clean$lower_um3<-(((as.numeric(paste(df_clean$lower))/2)^3)*1.333*pi)
  df_clean$upper_um3<-(((as.numeric(paste(df_clean$upper))/2)^3)*1.333*pi)
  df_clean$bin_width_um3<-df_clean$upper_um3 - df_clean$lower_um3
  
  df_clean <- df_clean %>% filter(Depth > 0) %>% filter(concentration_uL_L> 0)#[df_clean$depth_m > 0, ]
  
  return(df_clean)
}

#gets midpoints of bins obtained by the cut() function
midpoints <- function(x, dp=1){
  lower <- as.numeric(gsub(',.*','', gsub('\\(|\\[|\\)|\\]','', x)))# two regex searches: first gsub searches for ( or [ or ) or ) and removes it. Operators: |=or, \=silence (escape)
  upper <- as.numeric(gsub('.*,','', gsub('\\(|\\[|\\)|\\]','', x))) #the other regex search is to find any characters after or before the comma, and remove them 
  return(round(lower+(upper-lower)/2, dp))
}

depth_bin_calc <-function(df, depth_var='Depth', depth_bin_width=5){
  # assing depth bins to data.  initial depth binning will be each 5 m
  df<-df %>%mutate(Depth= as.numeric(.[[depth_var]])) %>% mutate(Depth = replace(Depth, Depth>100, 100)) %>% mutate(depth_binned = cut(.[[depth_var]], breaks=seq(0, 100, by = depth_bin_width)))
  df$depth_mid<-midpoints(df$depth_binned)
  return(df)
}

latlon_depth_binning <- function(df,bin_size){
  #test spatial binning (in lat/lon) will be one each kilometer. If ISIIS is towed at 2.5m/s, it takes 400 seconds(~7 minutes) to travel 1000 m
  time_travelled_min <-(((bin_size)*1000)/2.5)/60 # bin_size is in Km, will be transformed to meters to then use ship speed to get seconds travelled to get the desired bin size and then transformed to minutes
  str_time_travelled <- paste(as.character(round(time_travelled_min, 0)), 'min', sep = ' ')
  df$Year_Month_day_hr <- format(df$dateTime, "%Y-%m-%d %H:%M")
  
  df <- df %>%
    mutate(dateTime = ymd_hms(dateTime)) %>%
    group_by(size_class_mid_um, size_class_mid_ul,bin_width_ul,depth_mid, season, dateTime = floor_date(dateTime, str_time_travelled)) %>% summarise(concentration_uL_L = sum(concentration_uL_L)) %>% 
    dplyr::arrange(dateTime) %>% mutate(dateTime = as.factor(dateTime))%>%mutate(depth_mid = as.factor(depth_mid)) %>%mutate(size_class_mid_um = as.numeric(as.character(size_class_mid_um))) 
  
  return(df)
}

particle_count <- function(df){
  # approximate particle count assuming sphericity
  df$particle_count_L <- df$concentration_uL_L /
    ((((df$size_class_mid_um/2)^3) * 1.333 * 3.14) * 1e-9)
  
  # round to nearest and coerce to integer
  df$particle_count_L <- as.integer(round(df$particle_count_L))
  
  return(df)
}



