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
library(optparse)
library(RANN)
library(clue)

base_dir<- getwd()
source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))

# Define options
option_list <- list(
  make_option(c("-c", "--Cruise"), type = "character", help = "Cruise: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t", "--Transect"), type = "character", help = "Cruise: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "time_Pacific,Latitude,Longitude,Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-l", "--transectLats"), type = "character", default = "ancillary/transect_lats.csv", help = "Path to CSV with reference transect latitudes", metavar = "character"),
  make_option(c("-b", "--shelfBreak"), type = "character", default = "ancillary/shelf_break_200.csv", help = "Path to CSV with shelf break longitudes per transect", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$Cruise)) {
  stop("Error: --Cruise ID is required")
}

# Variables for data management

#variables to load data:
Cruise <- opt$Cruise # Cruise<- 'w22'
Transect <- opt$Transect # Transect<- 'NH'

size_var<- opt$sizeVar # size_var<-'contour'

ISIIS_NBSS_dir <- file.path(base_dir, 'data', 'size_spectra_allParticles_5-10m-bins_thresholded', size_var, Cruise, Transect) # #use this for testing ISIIS_NBSS_dir<-"/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/size_spectra_allParticles_5-10m-bins_thresholded/contour/w23/NH"
df_list_ISIIS <- list.files(file.path(ISIIS_NBSS_dir), pattern = '\\.csv$', full.names = TRUE, recursive = TRUE)
df_list_ISIIS <- df_list_ISIIS[grepl("by_community", df_list_ISIIS, fixed = TRUE) ]

df_list_LISST<- file.path(base_dir, 'data/LISST_data/LISST_size_spectra_thresholded/5m_depthBin' ,paste('SPECTRA', Cruise, Transect, '5m-depthBin_LISST_size_spectra.csv', sep='_')) #LISST_NBSS_data <- '/nfs4/FW_HMSC/Cowen_Lab/corralma/ISIIS_NBSS_workflow_SPECTRA/data/LISST_data/LISST_size_spectra/5m_depthBin/SPECTRA_w23_NH_5m-depthBin_LISST_size_spectra.csv'

#processing variables
#Variables for post-processing
transect_lat_path <- file.path(base_dir, opt$transectLats) # transect_lat_path <- file.path(base_dir, "ancillary/transect_lats.csv")
shelf_break_path <- file.path(base_dir, opt$shelfBreak) # shelf_break_path <- file.path(base_dir, "ancillary/shelf_break_200.csv")
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "/ancillary/transect_shelf_land.csv")

env_cols <- c("Temperature","Salinity","Oxygen_ml_l","pH","PAR","chl_a_ug_l")

output_dir_figs <- file.path(base_dir, "figures", 'size_spectra', 'LISST-ISIIS_merged')
if (!dir.exists(output_dir_figs)) {
  dir.create(output_dir_figs, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs)
}



dat_ISIIS<- map_dfr( df_list_ISIIS, ~ read_csv( .x, show_col_types = FALSE))
dat_ISIIS <- dat_ISIIS %>%
  # keep the original for reference
  dplyr::rename(group_id_old = group_id) %>%
  # new globally consecutive id: increments whenever group_id_old changes
  dplyr::mutate(
    group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE)),
    ind_L = ind_m3/1000 # cubic meters to liters
  )

#dat_ISIIS <- dat_ISIIS %>%
  #mutate(time_Pacific = time_Pacific - hours(1))

# rad LISST and create Cruise and Transect_ID
dat_LISST <- map_dfr(
  df_list_LISST,
  ~ {
    file_path <- .x
    fname <- basename(file_path)
    
    parts <- strsplit(fname, "_")[[1]]
    
    cruise      <- parts[2]   # ex: "w23"
    transect_id <- parts[3]   # ex: "OB"
    
    read_csv(file_path, show_col_types = FALSE) %>%
      mutate(
        Cruise      = cruise,
        Transect_ID = transect_id,
        ind_L = particle_count_L,
        ind_m3 = ind_L*1000
      )
  }
)

dat_LISST$ind_m3<-dat_LISST$particle_count_L / 0.001 

dat_LISST$count_uncertainty<-dpois(x=dat_LISST$ind_m3, lambda=dat_LISST$ind_m3)


dat_LISST <- dat_LISST %>%
  dplyr::select(
    group_id,
    size_class_mid_um,
    NB,
    cast_id,
    mean_time_Pacific,
    depth_mid,
    bin_width,
    size_class_mid_um3,
    bin_width_um3,
    concentration_uL_L,
    concentration_um3_m3,
    ind_m3,
    Latitude,
    Longitude,
    count_uncertainty,
    log_size_class_mid_um3,
    log_NB,
    n_size_classes
  ) %>%
  dplyr::rename(
    time_Pacific      = mean_time_Pacific, 
    ECD_um_mid        = size_class_mid_um,
    ECD_um_width      = bin_width,
    biovol_um3_mid    = size_class_mid_um3,
    biovol_um3_width  = bin_width_um3,
    logBiovol         = log_size_class_mid_um3,
    logNB             = log_NB
  )

dat_LISST <- dat_LISST %>%
  dplyr::mutate(
    group_id = as.integer(group_id),
    cast_id  = as.character(cast_id),
    
    time_Pacific = as.POSIXct(time_Pacific, tz = "UTC") + lubridate::hours(2),
    
    depth_mid = as.numeric(depth_mid),
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude),
    
    ECD_um_mid       = as.numeric(ECD_um_mid),
    ECD_um_width     = as.numeric(ECD_um_width),
    biovol_um3_mid   = as.numeric(biovol_um3_mid),
    biovol_um3_width = as.numeric(biovol_um3_width),
    
    NB        = as.numeric(NB),
    ind_m3     = as.numeric(ind_m3),
    logNB     = as.numeric(logNB),
    logBiovol = as.numeric(logBiovol),
    
    concentration_uL_L   = as.numeric(concentration_uL_L),
    concentration_um3_m3 = as.numeric(concentration_um3_m3),
    count_uncertainty    = as.numeric(count_uncertainty),
    
    n_size_classes = as.integer(n_size_classes)
  )


dat_LISST <- dat_LISST %>%
  dplyr::mutate(.row_id = dplyr::row_number()) %>%   # preserve original ordering
  dplyr::group_by(group_id) %>%
  dplyr::mutate(
    time_Pacific = as.POSIXct(
      mean(as.numeric(time_Pacific), na.rm = TRUE),
      origin = "1970-01-01",
      tz = attr(time_Pacific, "tzone")
    ),
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude  = mean(Latitude,  na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(.row_id) %>%
  dplyr::select(-.row_id)

#plot LISST and ISIIS by lonigutde
p_lon<- ggplot() +
  geom_point(
    data = dat_LISST,
    aes(x = Longitude, y = depth_mid),
    color = "blue",
    size = 1.5,
    alpha = 0.3
  ) +
  geom_point(
    data = dat_ISIIS,
    aes(x = Longitude, y = depth_mid),
    color = "green",
    size = 1.5,
    alpha = 0.3
  ) +
  scale_y_reverse() +
  labs(x = "Longitude", y = "Depth (m)")

p_lon

out_path <- file.path(output_dir_figs, paste(paste("LISST-ISIIS_data_coverage_longitude", Cruise, Transect, sep="_"), ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(filename = out_path, plot = p_lon, width = 9, height = 7, units = "in", dpi = 300, bg = "white")


#plot lISST and ISIIS by time:

p_time<- ggplot() +
  geom_point(
    data = dat_LISST,
    aes(x = time_Pacific, y = depth_mid),
    color = "blue",
    size = 1.5,
    alpha = 0.7
  ) +
  geom_point(
    data = dat_ISIIS,
    aes(x = time_Pacific, y = depth_mid),
    color = "green",
    size = 1.5,
    alpha = 0.7
  ) +
  scale_y_reverse() +
  labs(x = "Time (UTC)", y = "Depth (m)")

out_path <- file.path(output_dir_figs, paste(paste("LISST-ISIIS_data_coverage_Time", Cruise, Transect, sep="_"), ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(filename = out_path, plot = p_time, width = 9, height = 7, units = "in", dpi = 300, bg = "white")

#NBSS columns
lisst_spec_cols <- c(
  "ECD_um_mid","ECD_um_width","biovol_um3_mid","biovol_um3_width",
  "NB","ind_m3","logNB","logBiovol",
  "concentration_uL_L","concentration_um3_m3",#"count_uncertainty",
  "n_size_classes"
)

isiis_spec_cols <- c(
  "ECD_um_mid","ECD_um_width","biovol_um3_mid","biovol_um3_width",
  "NB","ind_m3","logNB","logBiovol", #"count_uncertainty",
  "n_size_classes"
)

#data frame of one row per isiis sample
isiis_sample_cols <- c("group_id","depth_mid","time_Pacific","Longitude","Latitude")

isiis_samples <- dat_ISIIS %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(isiis_sample_cols)), .keep_all = TRUE) %>%
  dplyr::mutate(isiis_sample_id = dplyr::row_number())

  isiis_meta <- dat_ISIIS %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(
    depth_mid    = dplyr::first(depth_mid),
    Longitude    = dplyr::first(Longitude),
    Latitude     = dplyr::first(Latitude),
    time_Pacific = dplyr::first(time_Pacific),
    Cruise       = dplyr::first(Cruise),
    Transect_ID  = dplyr::first(Transect_ID),
    cast_id      = dplyr::first(cast_id),

    Temperature  = dplyr::first(Temperature),
    Salinity     = dplyr::first(Salinity),
    Oxygen_ml_l  = dplyr::first(Oxygen_ml_l),
    pH           = dplyr::first(pH),
    PAR          = dplyr::first(PAR),
    chl_a_ug_l   = dplyr::first(chl_a_ug_l),

    .groups = "drop"
  )

# lisst metadata
lisst_meta <- dat_LISST %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(
    lisst_group_id = dplyr::first(group_id),
    depth_mid      = dplyr::first(depth_mid),
    Longitude      = dplyr::first(Longitude),
    Latitude       = dplyr::first(Latitude),
    time_Pacific   = dplyr::first(time_Pacific),
    .groups = "drop"
  )

match_rann_complete_one_to_one <- function(isiis_df, lisst_df,
                                           lon_scale = 1e-3,
                                           depth_scale = 5,
                                           time_scale = 60,
                                           max_dlon = Inf,
                                           max_ddepth = Inf,
                                           max_dt_sec = Inf,
                                           big_cost = 1e12) {
  
  isiis_df <- isiis_df %>%
    dplyr::mutate(
      x_lon   = as.numeric(Longitude) / lon_scale,
      x_depth = as.numeric(depth_mid) / depth_scale,
      x_time  = as.numeric(time_Pacific) / time_scale
    ) %>%
    dplyr::filter(is.finite(x_lon), is.finite(x_depth), is.finite(x_time))
  
  lisst_df <- lisst_df %>%
    dplyr::mutate(
      x_lon   = as.numeric(Longitude) / lon_scale,
      x_depth = as.numeric(depth_mid) / depth_scale,
      x_time  = as.numeric(time_Pacific) / time_scale
    ) %>%
    dplyr::filter(is.finite(x_lon), is.finite(x_depth), is.finite(x_time))
  
  if (nrow(lisst_df) < nrow(isiis_df)) {
    stop("Not enough LISST samples for one-to-one matching: ",
         nrow(lisst_df), " LISST vs ", nrow(isiis_df), " ISIIS.")
  }
  
  isiis_mat <- as.matrix(isiis_df[, c("x_lon", "x_depth", "x_time")])
  lisst_mat <- as.matrix(lisst_df[, c("x_lon", "x_depth", "x_time")])
  
  # Use all LISST samples as candidates.
  # This avoids losing matches because k was too small.
  nn <- RANN::nn2(
    data = lisst_mat,
    query = isiis_mat,
    k = nrow(lisst_mat)
  )
  
  cost_mat <- matrix(
    big_cost,
    nrow = nrow(isiis_df),
    ncol = nrow(lisst_df)
  )
  
  for (i in seq_len(nrow(isiis_df))) {
    
    idx <- nn$nn.idx[i, ]
    dist <- nn$nn.dists[i, ]
    
    dlon <- abs(lisst_df$Longitude[idx] - isiis_df$Longitude[i])
    ddepth <- abs(lisst_df$depth_mid[idx] - isiis_df$depth_mid[i])
    dt_sec <- abs(as.numeric(
      difftime(lisst_df$time_Pacific[idx], isiis_df$time_Pacific[i], units = "secs")
    ))
    
    keep <- is.finite(dist) &
      dlon <= max_dlon &
      ddepth <= max_ddepth &
      dt_sec <= max_dt_sec
    
    cost_mat[i, idx[keep]] <- dist[keep]
  }
  
  n_no_candidate <- sum(apply(cost_mat, 1, min, na.rm = TRUE) >= big_cost)
  
  if (n_no_candidate > 0) {
    stop(
      n_no_candidate,
      " ISIIS samples have no valid LISST candidate under current max_dlon/max_ddepth/max_dt_sec constraints."
    )
  }
  
  assignment <- clue::solve_LSAP(cost_mat)
  
  assigned_lisst_row <- as.integer(assignment)
  assigned_cost <- cost_mat[cbind(seq_len(nrow(isiis_df)), assigned_lisst_row)]
  
  if (any(assigned_cost >= big_cost)) {
    stop("Assignment failed: at least one ISIIS sample was assigned only to a forbidden/high-cost LISST candidate.")
  }
  
  tibble::tibble(
    group_id = isiis_df$group_id,
    lisst_group_id = lisst_df$lisst_group_id[assigned_lisst_row],
    dlon = abs(lisst_df$Longitude[assigned_lisst_row] - isiis_df$Longitude),
    ddepth = abs(lisst_df$depth_mid[assigned_lisst_row] - isiis_df$depth_mid),
    dt_sec = abs(as.numeric(
      difftime(lisst_df$time_Pacific[assigned_lisst_row], isiis_df$time_Pacific, units = "secs")
    )),
    cost = assigned_cost
  )
}


# mapping table with unique lisst_group_id
# lisst_meta: lisst_group_id, depth_mid, Longitude, time_Pacific, ...

mapping_time_space <- match_rann_complete_one_to_one(
  isiis_df = isiis_meta,
  lisst_df = lisst_meta,
  lon_scale = 1e-3,
  depth_scale = 5,
  time_scale = 60,
  max_dlon = Inf,
  max_ddepth = Inf,
  max_dt_sec = Inf
)


message("ISIIS groups: ", dplyr::n_distinct(isiis_meta$group_id))
message("LISST groups: ", dplyr::n_distinct(lisst_meta$lisst_group_id))
message("Matched ISIIS groups: ", dplyr::n_distinct(mapping_time_space$group_id))
message("Matched LISST groups: ", dplyr::n_distinct(mapping_time_space$lisst_group_id))

stopifnot(dplyr::n_distinct(mapping_time_space$group_id) == dplyr::n_distinct(isiis_meta$group_id))
stopifnot(!anyDuplicated(mapping_time_space$lisst_group_id))

isiis_to_lisst <- isiis_meta %>%
  dplyr::left_join(mapping_time_space, by = c("group_id")) %>%
  dplyr::left_join(
    lisst_meta %>%
      dplyr::select(lisst_group_id, lisst_Longitude = Longitude, lisst_depth_mid = depth_mid, lisst_time_Pacific = time_Pacific),
    by = "lisst_group_id"
  ) %>%
  dplyr::mutate(
    dt_hours = dt_sec / 3600
  )

summary(mapping_time_space$ddepth)
summary(mapping_time_space$dt_sec)
summary(mapping_time_space$dlon)

#isiis_to_lisst %>%
#dplyr::filter(group_id %in% 75:79) %>%
#dplyr::select(group_id, depth_mid, Longitude, lisst_group_id, lisst_depth_mid, lisst_Longitude, dlon, dt_hours)

# --------------------------------------------------------------------
# Keep ONLY matched pairs (drop ISIIS groups with no LISST match)
# --------------------------------------------------------------------
matched_pairs <- isiis_to_lisst %>%
  dplyr::filter(!is.na(lisst_group_id)) %>%
  dplyr::distinct(group_id, lisst_group_id)

matched_group_ids <- matched_pairs$group_id

# ----- LISST side (already mostly matched, but make it explicit) -----
lisst_long <- dat_LISST %>%
  dplyr::rename(lisst_group_id = group_id) %>%
  dplyr::inner_join(matched_pairs, by = "lisst_group_id") %>%  # ensures LISST kept only if used
  dplyr::left_join(
    isiis_meta,
    by = "group_id",
    suffix = c("_LISST", "_ISIIS")
  ) %>%
  dplyr::mutate(
    depth_mid    = depth_mid_ISIIS,
    Longitude    = Longitude_ISIIS,
    Latitude     = Latitude_ISIIS,
    time_Pacific = time_Pacific_ISIIS,
    cast_id      = cast_id_ISIIS,
    source       = "LISST"
  ) %>%
  dplyr::select(
    -depth_mid_LISST, -Longitude_LISST, -Latitude_LISST, -time_Pacific_LISST,
    -depth_mid_ISIIS, -Longitude_ISIIS, -Latitude_ISIIS, -time_Pacific_ISIIS,
    -cast_id_LISST, -cast_id_ISIIS
  )

# Ensure order of lisst_long (as you had)
lisst_long <- lisst_long %>%
  dplyr::select(
    group_id,
    depth_mid,
    Longitude,
    Latitude,
    time_Pacific,
    Cruise,
    Transect_ID,
    cast_id,
    dplyr::all_of(lisst_spec_cols),
    source
  )

# ----- ISIIS side: FILTER to matched group_ids ONLY -----
isiis_long <- dat_ISIIS %>%
  dplyr::filter(group_id %in% matched_group_ids) %>%   # <<<<<< this is the key change
  dplyr::mutate(
    concentration_uL_L   = NA_real_,
    concentration_um3_m3 = NA_real_
  ) %>%
  dplyr::select(
    group_id,
    depth_mid,
    Longitude,
    Latitude,
    time_Pacific,
    Cruise,
    Transect_ID,
    cast_id,
    dplyr::all_of(lisst_spec_cols)
  ) %>%
  dplyr::mutate(source = "ISIIS")

# Now merged_long truly contains only matched LISST–ISIIS pairs
merged_long <- dplyr::bind_rows(lisst_long, isiis_long) %>%
  dplyr::arrange(group_id)



# the catch now is what to do with the overlapping biovolume and size bins, approach will be:
# 1) find matches between LISST size bins and ISIIS size bins,
# 2) if two or more LISST NB values fall within an ISIIS size bin, sum their BIOVOLUME and then recalculate NB with the  corresponding ISIIS bin width
# 3) then, since the NBs now have the exact same width and midpoint, average them 

# this impies that I need total biovolume per size class. I will not care about ind_L for now, due to the weird way in which it is calculated for LISST data
merged_long2 <- merged_long %>%
  dplyr::mutate(
    biovolume_um3 = NB * biovol_um3_width
  ) %>%
  dplyr::select(-dplyr::any_of("ind_L"))


# Split
# Split
isiis_only <- merged_long2 %>% filter(source == "ISIIS")
lisst_only <- merged_long2 %>% filter(source == "LISST")


# ------------------------------------------------------------------------------
# 1) Determine overlap bounds [b_l, b_h] *per group_id*
# ------------------------------------------------------------------------------
# These are the sizes over which BOTH instruments can, in principle, contribute.
# If a group has no overlap, we will fall back to simply binding rows.
overlap_bounds <- merged_long2 %>%
  group_by(group_id, source) %>%
  summarise(
    min_b = min(ECD_um_mid, na.rm = TRUE),
    max_b = max(ECD_um_mid, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = source,
    values_from = c(min_b, max_b)
  ) %>%
  transmute(
    group_id,
    b_l = pmax(min_b_ISIIS, min_b_LISST, na.rm = TRUE),
    b_h = pmin(max_b_ISIIS, max_b_LISST, na.rm = TRUE),
    has_overlap = is.finite(b_l) & is.finite(b_h) & (b_h > b_l)
  )

# ------------------------------------------------------------------------------
# 2) Rebin LISST biovolume onto ISIIS bins (ONLY for LISST bins inside ISIIS bins)
# ------------------------------------------------------------------------------
# We map each LISST bin midpoint to the ISIIS bin whose extent contains it.
# This converts LISST from its native binning to the ISIIS binning, *conserving*
# biovolume (bin totals) because we sum LISST bin totals falling into each ISIIS bin.
#
# Important: We use the ISIIS bin *extent*:
#    [ECD_um_mid - ECD_um_width/2, ECD_um_mid + ECD_um_width/2]
#
# This step is about consistent binning, not weighting yet.
lisst_to_isiis_pairs <- lisst_only %>%
  dplyr::select(group_id,
         ECD_um_mid_LISST = ECD_um_mid,
         biovolume_um3_LISST = biovolume_um3,
         ind_m3_LISST = ind_m3) %>%
  inner_join(
    isiis_only %>% 
      dplyr::select(group_id, ECD_um_mid, ECD_um_width),   # <-- keep width only if you need in_bin
    by = "group_id",
    relationship = "many-to-many"
  ) %>%
  dplyr::mutate(
    in_bin = (ECD_um_mid_LISST >= (ECD_um_mid - ECD_um_width/2)) &
      (ECD_um_mid_LISST <= (ECD_um_mid + ECD_um_width/2))
  ) %>%
  dplyr::filter(in_bin)


# Sum LISST biovolume that falls into each ISIIS bin
lisst_rebinned_to_isiis <- lisst_to_isiis_pairs %>%
  dplyr::group_by(group_id, ECD_um_mid) %>%
  dplyr::summarise(
    biovolume_um3_LISST_sum = sum(biovolume_um3_LISST, na.rm = TRUE),
    ind_m3_LISST_sum = sum(ind_m3_LISST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    isiis_only %>% distinct(group_id, ECD_um_mid, biovol_um3_width),
    by = c("group_id", "ECD_um_mid")
  ) %>%
  dplyr::mutate(
    NB_LISST_rebinned = biovolume_um3_LISST_sum / biovol_um3_width,
    ind_m3_LISST_rebinned = ind_m3_LISST_sum 
  ) %>%
  dplyr::select(group_id, ECD_um_mid, NB_LISST_rebinned, biovolume_um3_LISST_sum, ind_m3_LISST_rebinned, ind_m3_LISST_sum)


# ------------------------------------------------------------------------------
# 3) Compute Haëntjens size-based weights at the ISIIS bin midpoints
# ------------------------------------------------------------------------------
# We evaluate weights at b = ECD_um_mid (the target bin mid).
#
# In overlap [b_l, b_h]:
#   w_ISIIS(b) = (b - b_l) / (b_h - b_l)
#   w_LISST(b) = (b_h - b) / (b_h - b_l)
#
# Outside overlap:
#   - if b < b_l: LISST only (w_LISST=1, w_ISIIS=0) BUT note: those bins are
#     typically only present in LISST and will be kept separately.
#   - if b > b_h: ISIIS only (w_ISIIS=1, w_LISST=0)
#
# Because we are merging onto ISIIS bins, we only need weights for ISIIS bins.
isiis_with_weights <- isiis_only %>%
  left_join(overlap_bounds, by = "group_id") %>%
  mutate(
    b = ECD_um_mid,
    # Default: no overlap -> NA weights (we will not blend)
    w_ISIIS = if_else(has_overlap, (b - b_l) / (b_h - b_l), NA_real_),
    w_LISST = if_else(has_overlap, (b_h - b) / (b_h - b_l), NA_real_),
    
    # Clamp weights to [0,1] to handle numerical edge cases at boundaries
    w_ISIIS = pmin(pmax(w_ISIIS, 0), 1),
    w_LISST = pmin(pmax(w_LISST, 0), 1)
  )

# ------------------------------------------------------------------------------
# 4) Blend rebinned LISST NB with ISIIS NB using size-based weights
# ------------------------------------------------------------------------------
# This is the direct analogue of Haëntjens Eq. (7), applied to NB on the same bins:
#
#   NB_merged = w_LISST(b) * NB_LISST_rebinned + w_ISIIS(b) * NB_ISIIS
#
# We do this ONLY where both exist (i.e., where NB_LISST_rebinned is non-NA).
# If no LISST contribution for that ISIIS bin, keep ISIIS as-is.
isiis_blended <- isiis_with_weights %>%
  left_join(lisst_rebinned_to_isiis, by = c("group_id", "ECD_um_mid")) %>%
  mutate(
    # Blend NB where LISST rebinned NB exists
    NB = if_else(
      !is.na(NB_LISST_rebinned) & !is.na(w_LISST) & !is.na(w_ISIIS),
      (w_LISST * NB_LISST_rebinned) + (w_ISIIS * NB),
      NB
    ),
    #Blend ind_m3 where LISST rebinned ind_m3 exists
    ind_m3 = if_else(
      !is.na(ind_m3_LISST_rebinned) & !is.na(w_LISST) & !is.na(w_ISIIS),
      (w_LISST * ind_m3_LISST_rebinned) + (w_ISIIS * ind_m3),
      ind_m3
    ),
    #Set source once: mark LISST-ISIIS if EITHER variable was blended
    source = if_else(
      (!is.na(NB_LISST_rebinned) | !is.na(ind_m3_LISST_rebinned)) &
        !is.na(w_LISST) & !is.na(w_ISIIS),
      "LISST-ISIIS",
      source
    ),
    #Recompute totals after blending (still derived from NB)
    biovolume_um3 = NB * biovol_um3_width
  ) %>%
  dplyr::select(
    -NB_LISST_rebinned,
    -biovolume_um3_LISST_sum,
    -ind_m3_LISST_rebinned,
    -ind_m3_LISST_sum,
    -b, -b_l, -b_h, -has_overlap, -w_LISST, -w_ISIIS
  )


# ------------------------------------------------------------------------------
# 5) Keep LISST-only bins that were NOT consumed into ISIIS bins
# ------------------------------------------------------------------------------
# Any LISST bin midpoint that fell into an ISIIS bin was "rebinned" and is now
# represented in the blended ISIIS row for that bin. To avoid double counting,
# remove those LISST rows and keep only the remaining LISST bins.
lisst_consumed <- lisst_to_isiis_pairs %>%
  distinct(group_id, ECD_um_mid = ECD_um_mid_LISST)

lisst_remaining <- lisst_only %>%
  anti_join(lisst_consumed, by = c("group_id", "ECD_um_mid"))

# ------------------------------------------------------------------------------
# 6) Final merged table: LISST-only (outside overlap / unmatched) + ISIIS blended
# ------------------------------------------------------------------------------
merged_long_no_overlap <- bind_rows(
  lisst_remaining,
  isiis_blended
) %>%
  arrange(group_id, depth_mid, ECD_um_mid)

# merged_long_no_overlap is the Haëntjens-style merge on ISIIS binning with
# size-based tapering weights in the overlap region.

merged_long_no_overlap$count_uncertainty<-dpois(x=merged_long_no_overlap$ind_m3, lambda=merged_long_no_overlap$ind_m3)

threshold_count<- 0.3
## --THRESHOLDING! --
#based on counts
merged_long_no_overlap <- merged_long_no_overlap[which(merged_long_no_overlap$count_uncertainty <= threshold_count), ]

# based on Reynolds et al. (2010) particle sizes <3.2 um  are not reliable
#Reynolds, R. A., D. Stramski, V. M. Wright, and S. B.Wozniak. 2010. Measurements and characterization ofparticle size distributions in coastal waters. J. Geophys.Res. 115: C08024. doi:10.1029/2009JC005930
# for our size classes, this applies where ECD_um_mid < 3.66 um, since that size class contains particles from 3.3 to 3.9 um.
merged_long_no_overlap <- merged_long_no_overlap[which(merged_long_no_overlap$ECD_um_mid >3.2), ]

# based on three empty consecutive size classes
merged_long_no_overlap <- merged_long_no_overlap %>%
   mutate(ECD_um_mid = as.numeric(ECD_um_mid)) %>%
   arrange(group_id, ECD_um_mid) %>%
   group_by(group_id) %>%
   group_modify(function(.x, .g) {
     na_vec <- is.na(.x$NB)
     r      <- rle(na_vec)
     pos    <- which(r$values & r$lengths >= 3)[1]        # first 3-NA run
     keep_n <- if (is.na(pos)) nrow(.x) else sum(r$lengths[seq_len(pos - 1)])
     if (keep_n == 0) .x[0, , drop = FALSE] else .x[seq_len(keep_n), , drop = FALSE]
   }) %>%
   ungroup()


#diagnostic plot
src_cols <- c(
  "LISST"        = "blue",
  "ISIIS"        = "darkgreen",
  "LISST-ISIIS"  = "purple"
)

ggplot(
  merged_long %>% dplyr::filter(group_id == 1),
  aes(x = ECD_um_mid, y = NB, color = source)
) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = src_cols, drop = FALSE) +
  labs(
    x = "ECD (µm)",
    y = "NB",
    color = "Source"
  ) +
  theme_bw()


ggplot(
  merged_long_no_overlap %>% dplyr::filter(group_id == 1),
  aes(x = ECD_um_mid, y = NB, color = source)
) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = src_cols, drop = FALSE) +
  labs(
    x = "ECD (µm)",
    y = "NB",
    color = "Source"
  ) +
  theme_bw()


# Fit
merged_long_no_overlap <- merged_long_no_overlap %>%
  lin_fit_func(., filter_min_size_class = TRUE, n_min_size_class = 4)

# assign metadata
merged_long_no_overlap <- merged_long_no_overlap  %>% 
  assign_transect_name(transect_lat_path, lat_col_name="Latitude") %>% 
  assign_shelf_location(shelf_break_path, lon_col_name='Longitude') %>% 
  assign_season(time_col = "time_Pacific") %>% 
  assign_day_night(crepuscular = FALSE, time_col= "time_Pacific", lat_col= "Latitude", long_col= "Longitude", tz = "America/Los_Angeles") %>% 
  assign_relative_distance(land_ref_path=shore_path, lat_col= 'Latitude', long_col='Longitude')

#add env variables back
merged_long_no_overlap <- merged_long_no_overlap %>%
  dplyr::left_join(
    isiis_meta %>% dplyr::select(group_id, dplyr::all_of(env_cols)),
    by = "group_id"
  )


out_path_merged <- file.path(base_dir, 'data', 'ISIIS-LISST_mergedNBSS_5m_allParticles', paste0(paste( Cruise, Transect, "LISST-ISIIS_NBSS", sep="_"), ".csv"))
dir.create(dirname(out_path_merged), recursive = TRUE, showWarnings = FALSE)
write.csv(merged_long_no_overlap, out_path_merged, row.names=FALSE)
