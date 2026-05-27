# Goal: to produce a unified dataset with all the LISST data, with accompanying metadata (date, time, depth, temperature, size class, abundance).

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
library(optparse)
library(fs)
library(glue)
library(stringr)

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
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

source(file.path(base_dir, 'SCRIPTS','functions', 'LISST_funcs.R'))

option_list <- list(
  make_option(c("-d", "--depthBin"), type = "integer", default = 5, help = "depth bin width to integrate the data", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

depth_bin<- opt$depthBin

#load all LISST data
df_full <- suppressMessages(
  readr::read_csv(file.path(base_dir, "data", "LISST_data", "LISST_raw_all.csv"),
                  show_col_types = FALSE,
                  col_select = -1))

#load file with cruise and transect time intervals
df_date_intervals<- read_csv(file.path(base_dir, "ancillary", "SPECTRA-MBON_CruiseTransect_timeRanges.csv"), show_col_types = FALSE)

df_date_intervals <- df_date_intervals %>%
  mutate(subset_id = paste(Project, Cruise, Transect, sep = "_"))

df_date_intervals <- df_date_intervals %>%
  mutate(
    start = dmy_hm(start_time, tz = "UTC"),
    end   = dmy_hm(end_time,   tz = "UTC")
  )

mid_size_class <- as.data.frame(transpose(list(read.csv(file.path(base_dir, "ancillary", "LISST_size_classes_um.csv"), sep = ',')$size_class_mid_um)))

# Changes to match ISIIS size spectra processing: data will be grouped in 2 m depth bins AND BY SIZE CLASS, just as the community NBSS, get average time, and sum of relevant
# IMPOSSIBLE AT THIS STAGE, CAST ID IS ALSO NEEDED
  


for (id in unique(df_date_intervals$subset_id)){

    row <- df_date_intervals[df_date_intervals$subset_id == id, ][1, ]
    start <- row$start
    end   <- row$end

 
    df_base <-  df_full %>%
                  arrange(dateTime) %>%
                  filter(dateTime >= start, dateTime< end)  %>% 
                  ungroup() 
    
    if (nrow(df_base) == 0L) {
       message(sprintf("No LISST data for %s between %s and %s — skipping.",id,format(start, "%Y-%m-%d %H:%M:%S %Z"), format(end,   "%Y-%m-%d %H:%M:%S %Z")))
    next  # go to next id in the loop
    }

    message(glue("Subsetted LISST data for {id}, from {format(start, '%Y-%m-%d %H:%M:%S %Z')} to {format(end, '%Y-%m-%d %H:%M:%S %Z')}")) 


    # JOIN OF ENVIRO PARAMETERS
    env_data_list<- list.files(file.path(base_dir,'data','merge_environmental'), full.names = TRUE, recursive  = TRUE)
    
    env_data_all<- map_dfr( env_data_list, ~ readr::read_csv(.x, show_col_types = FALSE))
    
    env_data_all <-  env_data_all %>%
      mutate(Time = as.POSIXct(Time)) %>%
      arrange(Time) %>%
      filter(Time >= start, Time< end)  %>% 
      ungroup() 
    
    if (nrow(env_data_all) == 0L) {
       message(sprintf("No environmental data for %s between %s and %s — skipping.",id, format(start, "%Y-%m-%d %H:%M:%S %Z"), format(end,   "%Y-%m-%d %H:%M:%S %Z")))
    next  # go to next id in the loop
    } 
    
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
    DT <- as.data.table(df_base)
    DT[, dateTime := as.POSIXct(dateTime, tz = "UTC")]
    
    # exact/nearest join (no tolerance needed since dt_abs == 0)
    df_ <- env_small[DT, on = .(Time = dateTime), roll = "nearest"]
    
    # optional sanity: should be same number of rows
    stopifnot(nrow(df_) == nrow(df_base))
    
    ##apply correction factors:
    df_ <- df_ %>% apply_calibration_factors(.) %>%
      dplyr::rename(time_Pacific = Time) %>% # this is needed for the interplation
      rename_with(~ str_replace_all(., " ", "_")) %>% 
      interpolate_missing_values(., c("Temperature", "Conductivity", "Pressure", "Depth", "Salinity", "Altitude", "Sound_Velocity", "Horizontal_Speed", "Vertical_Speed", 
                                      "Chlorophyll_Wavelength", "Oxygen_ml_l", "pH", "PAR", "Latitude", "Longitude", "chl_a_ug_l")) 
    
    
    
    ## ---SIZE SPECTRUM CALCULATIONS, alternative is to threshold, and size spectra for the FINEST scale of LISSt data as possible:
    
    
    #assign groups to data in order to threshold:
    df_<- df_%>% particle_count(.)
    df_$count_uncertainty<-dpois(x=df_$particle_count_L, lambda=df_$particle_count_L)
    
    df_$concentration_um3_m3 <-(df_$concentration_um3_L)/0.001 # concentration in cubic meters
    
    df_$NB<-df_$concentration_um3_m3/df_$bin_width_um3
    setDT(df_)[, group_id := .GRP, .(time_Pacific, Depth)]
    df_$group_id<-as.factor(df_$group_id)
    setcolorder(df_, c("group_id", setdiff(names(df_), "group_id")))
    
    #qverage NBSS to extract mode:
    avg_NB_by_size <- df_ %>%
      dplyr::group_by(size_class_mid_um) %>%
      dplyr::summarise(NB_mean = mean(NB, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(size_class_mid_um)
    
    peak_idx <- which.max(avg_NB_by_size$NB_mean)
    peak_size_class <- avg_NB_by_size$size_class_mid_um[peak_idx]
    peak_NB_mean    <- avg_NB_by_size$NB_mean[peak_idx]
    
   # #THRESHOLD
    df_thres <- df_ %>%
      dplyr::filter(!(!is.na(count_uncertainty) & count_uncertainty > 0.3))
   
    ##now, detect empty size classes
    #mid_size_class <- sort(unique(as.numeric(mid_size_class)))
    #mid_size_class <- round(as.numeric(mid_size_class), 3)
    #
    #df_thres <- df_thres %>%
    #  dplyr::mutate(size_class_mid_um = round(as.numeric(size_class_mid_um), 3))
    #
    #df_filled <- df_thres %>%
    #  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
    #  group_by(group_id) %>%
    #  complete(size_class_mid_um = mid_size_class) %>%   # add missing sizes per group
    #  arrange(group_id, size_class_mid_um) %>%
    #  ungroup()
    #
    #df_filled <- df_filled %>%
    #  dplyr::relocate(NB, .before = 3)   
    #
    ##cut by maximum
    #df_threshold<- df_filled %>%
    #  dplyr::group_by(group_id) %>%
    #  dplyr::mutate(
    #    .maxNB   = max(NB, na.rm = TRUE),
    #    .cutoff  = max(size_class_mid_um[NB == .maxNB], na.rm = TRUE)  # largest size among ties
    #  ) %>%
    #  dplyr::ungroup() %>%
    #  dplyr::filter(size_class_mid_um >= .cutoff) %>%
    #  dplyr::select(-.maxNB, -.cutoff)
    
    #filter data after three consecutive size classes, removed since we were substantially reducing the amount of data
    #df_threshold <- df_threshold %>%
    #  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
    #  arrange(group_id, size_class_mid_um) %>%
    #  group_by(group_id) %>%
    #  group_modify(function(.x, .g) {
    #    na_vec <- is.na(.x$NB)
    #    r      <- rle(na_vec)
    #    pos    <- which(r$values & r$lengths >= 3)[1]        # first 3-NA run
    #    keep_n <- if (is.na(pos)) nrow(.x) else sum(r$lengths[seq_len(pos - 1)])
    #    if (keep_n == 0) .x[0, , drop = FALSE] else .x[seq_len(keep_n), , drop = FALSE]
    #  }) %>%
    #  ungroup()
    
    
    ## -- AFTER NORMALIZED BIOVOLUME CALCULATIONS PER FINEST SAMLPE AND THRESHOLDING, averaging NB by depth_mid and cast id
    
    # Normalized biovolume and size spectra will be averaged by depth_bin, so apply depth binning, despiking, and averaging:
    df_threshold <- df_thres %>% depth_bin_calc(., depth_bin_width=depth_bin) %>% dplyr::relocate(depth_mid, .before = 1) %>% dplyr::relocate(Depth, .before = 2)  
    
    
    
    df_cast_id_key <- df_threshold %>%
      dplyr::select(time_Pacific, depth_mid) %>%
      drop_na(time_Pacific, depth_mid) %>%
      distinct() %>%                  # unique combos of the three cols
      arrange(time_Pacific)%>%
      repeat_despike_depth(depth_bin_col = "depth_mid", time_col = "time_Pacific") 
    
    df_cast_id_key$cast_id <- assign_cast_id_base(df_cast_id_key$depth_mid)
    df_cast_id_key$depth_mid <- NULL
    
    df_threshold <- df_threshold %>%
      left_join(df_cast_id_key, by = c("time_Pacific")) %>% dplyr::relocate(cast_id, .before = 1)
    
    
    setDT(df_threshold)[, group_id := .GRP, .(cast_id, depth_mid)]
    df_threshold$group_id<-as.factor(df_threshold$group_id) 
    df_threshold <- df_threshold %>% dplyr::relocate(group_id, .before = 1)
    
    
    
    # keys
    grp_cols <- c("group_id","cast_id","depth_mid","size_class_mid_um", "bin_width", "size_class_mid_um3", "bin_width_um3", "season")
    
    # numeric fields to average (excluding time_Pacific)
    avg_cols <- c("Depth", "NB", "concentration_uL_L","concentration_um3_m3", "particle_count_L",
      "Conductivity","Pressure","Salinity","Sound_Velocity","Altitude",
      "Horizontal_Speed","Vertical_Speed","chl_a_ug_l",
      "Oxygen_ml_l","pH","PAR",
      "Temperature")
    median_cols<- c("Latitude","Longitude")
    
    max_cols <-c("n_size_classes")
    
    df_grouped <- df_threshold %>%
      dplyr::group_by(dplyr::across(all_of(grp_cols))) %>%
      dplyr::summarise(
        # keep this exactly as you wanted
        mean_time_Pacific =  time_Pacific[ceiling(n()/2)],#mean(time_Pacific, na.rm = TRUE),
        
        # means for avg_only
        dplyr::across(dplyr::any_of(avg_cols), ~ mean(as.numeric(.), na.rm = TRUE)),
       
        #median for latitude and longitude 
        dplyr::across(dplyr::any_of(median_cols), ~ median(as.numeric(.), na.rm = TRUE)),
        
        # max for max_cols
        dplyr::across(dplyr::any_of(max_cols), ~ max(as.numeric(.), na.rm = TRUE)),
        
        .groups = "drop"
      ) %>%
      dplyr::arrange(mean_time_Pacific, cast_id) %>%
      dplyr::select(cast_id, mean_time_Pacific, depth_mid, dplyr::everything())
    
    df_grouped <- df_grouped %>%
      dplyr::arrange(cast_id, group_id, size_class_mid_um)
    
    # after this, a new threshloding and size spectrum  calclulation needs to be performed:
    
    #second round of threshloding 
    df_grouped<- df_grouped%>% particle_count(.)
    df_grouped$count_uncertainty<-dpois(x=df_grouped$particle_count_L, lambda=df_grouped$particle_count_L)
    
    #
    ##qverage NBSS to extract mode:
    avg_NB_by_size <- df_grouped %>%
      dplyr::group_by(size_class_mid_um) %>%
      dplyr::summarise(NB_mean = mean(NB, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(size_class_mid_um)
    
    peak_idx <- which.max(avg_NB_by_size$NB_mean)
    peak_size_class <- avg_NB_by_size$size_class_mid_um[peak_idx]
    peak_NB_mean    <- avg_NB_by_size$NB_mean[peak_idx]
    
    ##THRESHOLD
    df_thres_grouped <- df_grouped %>%
      dplyr::filter(!(!is.na(count_uncertainty) & count_uncertainty > 0.3))
    
    
    #df_filled_grouped <- df_thres_grouped %>%
    #  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
    #  group_by(group_id) %>%
    #  complete(size_class_mid_um = mid_size_class) %>%   # add missing sizes per group
    #  arrange(group_id, size_class_mid_um) %>%
    #  ungroup()
    #
    #df_filled_grouped <- df_filled_grouped %>%
    #  dplyr::relocate(NB, .before = 3)   
    ##
    #cut by maximum
    #df_threshold_grouped<- df_filled_grouped %>%
    #  dplyr::group_by(group_id) %>%
    #  dplyr::mutate(
    #    .maxNB   = max(NB, na.rm = TRUE),
    #    .cutoff  = max(size_class_mid_um[NB == .maxNB], na.rm = TRUE)  # largest size among ties
    #  ) %>%
    #  dplyr::ungroup() %>%
    #  dplyr::filter(size_class_mid_um >= .cutoff) %>%
    #  dplyr::select(-.maxNB, -.cutoff)
   
    # based on Reynolds et al. (2010) particle sizes <3.2 um  are not reliable
    #Reynolds, R. A., D. Stramski, V. M. Wright, and S. B.Wozniak. 2010. Measurements and characterization ofparticle size distributions in coastal waters. J. Geophys.Res. 115: C08024. doi:10.1029/2009JC005930
# for our size classes, this applies where ECD_um_mid < 3.66 um, since that size class contains particles from 3.3 to 3.9 um.
    df_threshold_grouped <- df_thres_grouped[which(df_thres_grouped$size_class_mid_um >3.2), ]
 
    #filter data after three consecutive size classes, removed since we were reducind significantly the number of data points
    #df_threshold_grouped <- df_threshold_grouped %>%
    #  mutate(size_class_mid_um = as.numeric(size_class_mid_um)) %>%
    #  arrange(group_id, size_class_mid_um) %>%
    #  group_by(group_id) %>%
    #  group_modify(function(.x, .g) {
    #    na_vec <- is.na(.x$NB)
    #    r      <- rle(na_vec)
    #    pos    <- which(r$values & r$lengths >= 3)[1]        # first 3-NA run
    #    keep_n <- if (is.na(pos)) nrow(.x) else sum(r$lengths[seq_len(pos - 1)])
    #    if (keep_n == 0) .x[0, , drop = FALSE] else .x[seq_len(keep_n), , drop = FALSE]
    #  }) %>%
    #  ungroup()
    
    ## ---- FINAL SIZE SPECTRA FIT
    
    df_threshold_grouped$log_size_class_mid_um3<-log10(df_threshold_grouped$size_class_mid_um3)
    df_threshold_grouped$log_NB<-log10(df_threshold_grouped$NB)
    
    fitted_models = df_threshold_grouped %>% group_by(group_id) %>% do(model = lm(log_NB ~ log_size_class_mid_um3, data = .))
    # this is the structure of how the data is obstained: fitted_models[[2]][['group_id']]$coefficients[1] for intercept, coefficients[2] for slope
    #rownames(fitted_models) <-fitted_models$group_id
    
    
    for (n in 1:nrow(fitted_models)){
      g= fitted_models[[n,1]]
      #print(as.numeric(g))
      #print(fitted_models[[2]][[as.numeric(g)]]$coefficients[[1]])
      #print(fitted_models[[2]][[as.numeric(g)]]$coefficients[[2]])
      df_threshold_grouped$intercept[df_threshold_grouped$group_id==g]<-fitted_models[[2]][[n]]$coefficients[1]#intercept
      df_threshold_grouped$slope[df_threshold_grouped$group_id==g]<-fitted_models[[2]][[n]]$coefficients[2]#slope
      df_threshold_grouped$R2[df_threshold_grouped$group_id==g]<-summary(fitted_models[[2]][[n]])$r.squared
      df_threshold_grouped$n_size_classes[df_threshold_grouped$group_id==g]<-nrow(fitted_models[[2]][[n]]$model)# number of size classes to estimate the NBSS. 
    }
    
    df_threshold_grouped <- df_threshold_grouped %>%
      dplyr::filter(!is.na(NB), !is.na(n_size_classes), n_size_classes >= 4)
    
    #--- FOR EXPLORATION OF RESULTS:
    # 1) minimal frame, keep only finite/positive for log scales
    library(dplyr); library(ggplot2); library(scales)
    
    size_spectra_graph <- ggplot() +
      geom_line(
        data = df_threshold_grouped,
        aes(x = size_class_mid_um, y = NB, group = group_id, color = season),
        alpha = 0.1, linewidth = 0.1
      ) +
      geom_point(
        data = df_threshold_grouped,
        aes(x = size_class_mid_um, y = NB, color = season),
        size = 0.5, alpha = 0.5, shape = 16, position = position_dodge(0.2)
      ) +
      geom_line(
        data = df_threshold_grouped %>%
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
        y = expression(log[10](Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))),
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
    
    out_path_fig <- file.path(base_dir, "data", "LISST_data", "diagnostic_figs", paste(depth_bin, 'm_depthBin', sep=''), paste(id, paste(depth_bin, 'm-depthBin', sep=''),  "LISST_size-spectrum_thresholded.png", sep="_"))
    dir.create(dirname(out_path_fig), recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = out_path_fig, plot = size_spectra_graph, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
    message(sprintf("Saved figure: %s", out_path_fig))

    out_path_data <- file.path(base_dir, "data", "LISST_data", "LISST_size_spectra_thresholded", paste(depth_bin, 'm_depthBin', sep=''), paste(id, paste(depth_bin, 'm-depthBin', sep=''), "LISST_size_spectra.csv", sep="_"))
    dir.create(dirname(out_path_data), recursive = TRUE, showWarnings = FALSE) 
    write.csv(df_threshold_grouped, out_path_data, row.names = FALSE)
    message(sprintf("Saved data: %s", out_path_data))
}







