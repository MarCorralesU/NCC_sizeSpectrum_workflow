# Goal: to produce a unified dataset with all the LISST data, with accompanying metadata (date, time, depth, temperature, size class, abundance).

setwd('/Users/user/SPECTRA_proj/')
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
library(purrr)
library(dplyr)
library(fuzzyjoin)
library(vctrs)

# annoying version of base R cast_id function
assign_cast_id_base <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n <= 1L) return(rep.int(1L, n))
  cast  <- integer(n); cast[1L] <- 1L
  lastd <- 0L
  for (i in 2L:n) {
    d <- sign(x[i] - x[i - 1L])   # -1 down, 0 flat, +1 up
    if (d == 0L) d <- lastd       # carry direction across flats
    # start a new cast when direction flips (ignore first transition)
    cast[i] <- cast[i - 1L] + as.integer(i > 2L && d != lastd)
    lastd <- d
  }
  cast
}


base_dir<-getwd()
mid_size_class <- read.csv(
  file.path(base_dir, "ancillary", "LISST_size_classes_um.csv"),
  stringsAsFactors = FALSE
)[["size_class_mid_um"]]


source(file.path(base_dir, 'scripts','functions', 'LISST_funcs.R'))

funcs_dir <- file.path('~/SCRIPTS_backups', "ISIIS_NBSS_workflow_SPECTRA")
source(file.path(funcs_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


df_list<-Sys.glob(file.path(base_dir, "data", "LISST",  '**', 'processed_matlab', '*.csv') )
file_names <- basename(df_list)
df_list <- df_list[!duplicated(file_names)]
df_list <- df_list[!grepl("_rs", df_list)]

df_list_W22 <- df_list[grepl("W22", df_list)]


df_full<-data.frame()



df_full <- map_dfr(df_list, ~ format_LISST_files_custom_processed(.x),
                   .id = "source_file") 
df_full$dateTime<- as.POSIXct(df_full$dateTime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")


df_full$yearMonth<- format(df_full$dateTime, "%Y-%m")


write.csv(df_full, file.path(base_dir, "data", "LISST_raw_all.csv"), row.names = FALSE)

df_full <- suppressMessages(
  readr::read_csv(file.path(base_dir, "data", "LISST_raw_all.csv"),
                  show_col_types = FALSE,
                  col_select = -1))


for (m in unique(df_full$yearMonth)) {
  sel <- df_full$yearMonth == m
  message("Month: ", m)
  message("  start time was: ", min(df_full$dateTime[sel], na.rm = TRUE))
  message("  end time was: ",   max(df_full$dateTime[sel], na.rm = TRUE))
}

# Changes to match ISIIS size spectra processing: data will be grouped in 2 m depth bins AND BY SIZE CLASS, just as the community NBSS, get average time, and sum of relevant
# IMPOSSIBLE AT THIS STAGE, CAST ID IS ALSO NEEDED
  
#subset w23  NH data only:

#subset s22  NH data only:
start <- as.POSIXct("2022-07-24 18:30:00", tz = "UTC")
end   <- as.POSIXct("2022-07-26 02:50:00", tz = "UTC")  # end-exclusive is safer

df_w23_NH_base <-  df_full %>%
              arrange(dateTime) %>%
              filter(dateTime >= start, dateTime< end)  %>% 
              ungroup() 

# FUZZY JOIN OF ENVIRO PARAMETERS
env_data_list<- list.files('/Users/user/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA/data/merge_environmental', full.names = TRUE, recursive  = TRUE)

env_data_all<- map_dfr( env_data_list, ~ readr::read_csv(.x, show_col_types = FALSE))

env_data_all <-  env_data_all %>%
  mutate(Time = as.POSIXct(Time)) %>%
  arrange(Time) %>%
  filter(Time >= start, Time< end)  %>% 
  ungroup() 



wanted <- c(
  "Conductivity","Pressure","Salinity","Sound Velocity","Altitude",
  "Raw Oxygen [V]","Raw PAR [V]","Raw pH [V]",
  "Horizontal Speed","Vertical Speed","Chlorophyll Wavelength",
  "Chlorophyll","Latitude","Longitude"
)

# env: keep only what we need, align types, suffix to avoid name clashes
env_small <- as.data.table(env_data_all)[, c("Time", wanted), with = FALSE]
env_small[, Time := as.POSIXct(Time, tz = "UTC")]
setkey(env_small, Time)

# df: make sure time is POSIXct
DT <- as.data.table(df_w23_NH_base)
DT[, dateTime := as.POSIXct(dateTime, tz = "UTC")]

# exact/nearest join (no tolerance needed since dt_abs == 0)
df_w23_NH <- env_small[DT, on = .(Time = dateTime), roll = "nearest"]

# optional sanity: should be same number of rows
stopifnot(nrow(df_w23_NH) == nrow(df_w23_NH_base))

##apply correction factors:
df_w23_NH <- df_w23_NH %>% apply_calibration_factors(.) %>%
  dplyr::rename(time_Pacific = Time) %>% # this is needed for the interplation
  rename_with(~ str_replace_all(., " ", "_")) %>% 
  interpolate_missing_values(., c("Temperature", "Conductivity", "Pressure", "Depth", "Salinity", "Altitude", "Sound_Velocity", "Horizontal_Speed", "Vertical_Speed", 
                                  "Chlorophyll_Wavelength", "Oxygen_ml_l", "pH", "PAR", "Latitude", "Longitude", "chl_a_ug_l")) 



## ---SIZE SPECTRUM CALCULATIONS, alternative is to threshold, and size spectra for the FINEST scale of LISSt data as possible:


#assign groups to data in order to threshold:
df_w23_NH<- df_w23_NH%>% particle_count(.)
df_w23_NH$count_uncertainty<-dpois(x=df_w23_NH$particle_count_L, lambda=df_w23_NH$particle_count_L)

df_w23_NH$concentration_um3_m3 <-(df_w23_NH$concentration_um3_L)/0.001 # concentration in cubic meters

df_w23_NH$NB<-df_w23_NH$concentration_um3_m3/df_w23_NH$bin_width_um3
setDT(df_w23_NH)[, group_id := .GRP, .(time_Pacific, Depth)]
df_w23_NH$group_id<-as.factor(df_w23_NH$group_id)
setcolorder(df_w23_NH, c("group_id", setdiff(names(df_w23_NH), "group_id")))





#qverage NBSS to extract mode:
avg_NB_by_size <- df_w23_NH %>%
  dplyr::group_by(size_class_mid_um) %>%
  dplyr::summarise(NB_mean = mean(NB, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(size_class_mid_um)

peak_idx <- which.max(avg_NB_by_size$NB_mean)
peak_size_class <- avg_NB_by_size$size_class_mid_um[peak_idx]
peak_NB_mean    <- avg_NB_by_size$NB_mean[peak_idx]

#THRESHOLD
df_w23_NH_thres <- df_w23_NH %>%
  dplyr::filter(!(!is.na(count_uncertainty) & count_uncertainty > 0.3))

#now, detect empty size classes
mid_size_class <- sort(unique(as.numeric(mid_size_class)))
mid_size_class <- round(as.numeric(mid_size_class), 3)

df_w23_NH_thres <- df_w23_NH_thres %>%
  dplyr::mutate(size_class_mid_um = round(as.numeric(size_class_mid_um), 3))

df_w23_NH_filled <- df_w23_NH_thres %>%
  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
  group_by(group_id) %>%
  complete(size_class_mid_um = mid_size_class) %>%   # add missing sizes per group
  arrange(group_id, size_class_mid_um) %>%
  ungroup()

df_w23_NH_filled <- df_w23_NH_filled %>%
  dplyr::relocate(NB, .before = 3)   

#cut by maximum
df_w23_NH_threshold<- df_w23_NH_filled %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(
    .maxNB   = max(NB, na.rm = TRUE),
    .cutoff  = max(size_class_mid_um[NB == .maxNB], na.rm = TRUE)  # largest size among ties
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(size_class_mid_um >= .cutoff) %>%
  dplyr::select(-.maxNB, -.cutoff)

#filter data after three consecutive size classes
df_w23_NH_threshold <- df_w23_NH_threshold %>%
  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
  arrange(group_id, size_class_mid_um) %>%
  group_by(group_id) %>%
  group_modify(function(.x, .g) {
    na_vec <- is.na(.x$NB)
    r      <- rle(na_vec)
    pos    <- which(r$values & r$lengths >= 3)[1]        # first 3-NA run
    keep_n <- if (is.na(pos)) nrow(.x) else sum(r$lengths[seq_len(pos - 1)])
    if (keep_n == 0) .x[0, , drop = FALSE] else .x[seq_len(keep_n), , drop = FALSE]
  }) %>%
  ungroup()


## -- AFTER NORMALIZED BIOVOLUME CALCULATIONS PER FINEST SAMLPE AND THRESHOLDING, averaging NB by depth_mid and cast id


# Normalized biovolume and size spectra will be averaged by depth_bin, so apply depth binning, despiking, and averaging:
df_w23_NH_threshold <- df_w23_NH_threshold %>% depth_bin_calc(., depth_bin_width=2) %>% dplyr::relocate(depth_mid, .before = 1) %>% dplyr::relocate(Depth, .before = 2)  



df_w23_NH_cast_id_key <- df_w23_NH_threshold %>%
  dplyr::select(time_Pacific, depth_mid) %>%
  drop_na(time_Pacific, depth_mid) %>%
  distinct() %>%                  # unique combos of the three cols
  arrange(time_Pacific)%>%
  repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") 

df_w23_NH_cast_id_key$cast_id <- assign_cast_id_base(df_w23_NH_cast_id_key$depth_mid)
df_w23_NH_cast_id_key$depth_mid <- NULL

df_w23_NH_threshold <- df_w23_NH_threshold %>%
  left_join(df_w23_NH_cast_id_key, by = c("time_Pacific")) %>% dplyr::relocate(cast_id, .before = 1)


setDT(df_w23_NH_threshold)[, group_id := .GRP, .(cast_id, depth_mid)]
df_w23_NH_threshold$group_id<-as.factor(df_w23_NH_threshold$group_id) 
df_w23_NH_threshold <- df_w23_NH_threshold %>% dplyr::relocate(group_id, .before = 1)



# keys
grp_cols <- c("group_id","cast_id","depth_mid","size_class_mid_um", "bin_width", "size_class_mid_um3", "bin_width_um3", "season")

# numeric fields to average (excluding time_Pacific)
avg_cols <- c("Depth", "NB", "concentration_uL_L","concentration_um3_m3", "particle_count_L",
  "Conductivity","Pressure","Salinity","Sound_Velocity","Altitude",
  "Horizontal_Speed","Vertical_Speed","Chlorophyll_Wavelength",
  "Oxygen_ml_l","pH","PAR","Latitude","Longitude",
  "Temperature")

max_cols <-c("n_size_classes")

df_w23_NH_grouped <- df_w23_NH_threshold %>%
  dplyr::group_by(dplyr::across(all_of(grp_cols))) %>%
  dplyr::summarise(
    # keep this exactly as you wanted
    mean_time_Pacific = mean(time_Pacific, na.rm = TRUE),
    
    # means for avg_only
    dplyr::across(dplyr::any_of(avg_cols), ~ mean(as.numeric(.), na.rm = TRUE)),
    
    # max for max_cols
    dplyr::across(dplyr::any_of(max_cols), ~ max(as.numeric(.), na.rm = TRUE)),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_time_Pacific, cast_id) %>%
  dplyr::select(cast_id, mean_time_Pacific, depth_mid, dplyr::everything())

df_w23_NH_grouped <- df_w23_NH_grouped %>%
  dplyr::arrange(cast_id, group_id, size_class_mid_um)

# after this, a new threshloding and size spectrum  calclulation needs to be performed:

#second round of threshloding 
df_w23_NH_grouped<- df_w23_NH_grouped%>% particle_count(.)
df_w23_NH_grouped$count_uncertainty<-dpois(x=df_w23_NH_grouped$particle_count_L, lambda=df_w23_NH_grouped$particle_count_L)


#qverage NBSS to extract mode:
avg_NB_by_size <- df_w23_NH_grouped %>%
  dplyr::group_by(size_class_mid_um) %>%
  dplyr::summarise(NB_mean = mean(NB, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(size_class_mid_um)

peak_idx <- which.max(avg_NB_by_size$NB_mean)
peak_size_class <- avg_NB_by_size$size_class_mid_um[peak_idx]
peak_NB_mean    <- avg_NB_by_size$NB_mean[peak_idx]

#THRESHOLD
df_w23_NH_thres_grouped <- df_w23_NH_grouped %>%
  dplyr::filter(!(!is.na(count_uncertainty) & count_uncertainty > 0.3))


df_w23_NH_filled_grouped <- df_w23_NH_thres_grouped %>%
  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
  group_by(group_id) %>%
  complete(size_class_mid_um = mid_size_class) %>%   # add missing sizes per group
  arrange(group_id, size_class_mid_um) %>%
  ungroup()

df_w23_NH_filled_grouped <- df_w23_NH_filled_grouped %>%
  dplyr::relocate(NB, .before = 3)   

#cut by maximum
df_w23_NH_threshold_grouped<- df_w23_NH_filled_grouped %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(
    .maxNB   = max(NB, na.rm = TRUE),
    .cutoff  = max(size_class_mid_um[NB == .maxNB], na.rm = TRUE)  # largest size among ties
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(size_class_mid_um >= .cutoff) %>%
  dplyr::select(-.maxNB, -.cutoff)

#filter data after three consecutive size classes
df_w23_NH_threshold_grouped <- df_w23_NH_threshold_grouped %>%
  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
  arrange(group_id, size_class_mid_um) %>%
  group_by(group_id) %>%
  group_modify(function(.x, .g) {
    na_vec <- is.na(.x$NB)
    r      <- rle(na_vec)
    pos    <- which(r$values & r$lengths >= 3)[1]        # first 3-NA run
    keep_n <- if (is.na(pos)) nrow(.x) else sum(r$lengths[seq_len(pos - 1)])
    if (keep_n == 0) .x[0, , drop = FALSE] else .x[seq_len(keep_n), , drop = FALSE]
  }) %>%
  ungroup()


## ---- FINAL SIZE SPECTRA FIT

df_w23_NH_threshold_grouped$log_size_class_mid_um3<-log10(df_w23_NH_threshold_grouped$size_class_mid_um3)
df_w23_NH_threshold_grouped$log_NB<-log10(df_w23_NH_threshold_grouped$NB)

fitted_models = df_w23_NH_threshold_grouped %>% group_by(group_id) %>% do(model = lm(log_NB ~ log_size_class_mid_um3, data = .))
# this is the structure of how the data is obstained: fitted_models[[2]][['group_id']]$coefficients[1] for intercept, coefficients[2] for slope


for (n in 1:nrow(fitted_models)){
  g= fitted_models[[n,1]]
  df_w23_NH_threshold_grouped$intercept[df_w23_NH_threshold_grouped$group_id==g]<-fitted_models[[2]][[n]]$coefficients[1]#intercept
  df_w23_NH_threshold_grouped$slope[df_w23_NH_threshold_grouped$group_id==g]<-fitted_models[[2]][[n]]$coefficients[2]#slope
  df_w23_NH_threshold_grouped$R2[df_w23_NH_threshold_grouped$group_id==g]<-summary(fitted_models[[2]][[n]])$r.squared
  df_w23_NH_threshold_grouped$n_size_classes[df_w23_NH_threshold_grouped$group_id==g]<-nrow(fitted_models[[2]][[n]]$model)# number of size classes to estimate the NBSS. 
}

df_w23_NH_threshold_grouped <- df_w23_NH_threshold_grouped %>%
  dplyr::filter(!is.na(NB), !is.na(n_size_classes), n_size_classes >= 4)

#--- FOR EXPLORATION OF RESULTS:
# 1) minimal frame, keep only finite/positive for log scales
library(dplyr); library(ggplot2); library(scales)

size_spectra_graph <- ggplot() +
  geom_line(
    data = df_w23_NH_threshold_grouped,
    aes(x = size_class_mid_um, y = NB, group = group_id, color = season),
    alpha = 0.1, linewidth = 0.1
  ) +
  geom_point(
    data = df_w23_NH_threshold_grouped,
    aes(x = size_class_mid_um, y = NB, color = season),
    size = 0.5, alpha = 0.5, shape = 16, position = position_dodge(0.2)
  ) +
  geom_line(
    data = df_w23_NH_threshold_grouped %>%
      dplyr::mutate(size_class_mid_um = as.factor(size_class_mid_um)) %>%
      dplyr::group_by(season, size_class_mid_um) %>%
      dplyr::summarise(NB = median(NB, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(size_class_mid_um = as.numeric(as.character(size_class_mid_um))),
    aes(x = size_class_mid_um, y = NB, group = season, color = season),
    alpha = 1, linewidth = 0.5
  ) +
  scale_x_log10(breaks = c(2, 20, 200)) +
  scale_y_log10(labels = label_log(digits = 2)) +
  labs(
    y = expression(log[10](Normalized~Biovolume~(mu*l~L^{-1}~u*L^{-1}))),
    x = expression(log[10](Equivalent~spherical~diameter~(mu*m)))
  ) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    axis.text.x = element_text(family = "Times", size = 15),
    axis.title.x = element_text(family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 15),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey")
  )

size_spectra_graph



write.csv(df_w23_NH_threshold_grouped, paste(base_dir, 'data', 'LISST_NBSS_s22_NH.csv', sep='/'))









