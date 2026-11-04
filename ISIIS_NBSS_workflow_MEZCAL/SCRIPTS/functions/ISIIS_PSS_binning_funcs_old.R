library(data.table)
library(dplyr)
library(magrittr)
library(tidyr)
library(NLP)
library(here)
library(lubridate)
library(runner)
library(suncalc)
library(purrr)

#to get Laf columns
get_LaF_colnames <- function(csv_path) {
  # Read only the header
  header_line <- readLines(csv_path, n = 1)
  col_names <- strsplit(header_line, ",")[[1]]
  
  # Set all column types to "character"
  col_types <- rep("character", length(col_names))
  
  return(list(
    col_names = col_names,
    col_types = col_types
  ))
}

convert_latlong_dms <- function(df, lat_col = "Lat", long_col = "Long") {
  dms_to_dd <- function(dms) {
    # Extract degrees, minutes, seconds, direction
    parts <- stringr::str_match(dms, "\\s*(\\d+)°\\s*(\\d+)'\\s*([0-9.]+)([NSEW])")
    deg <- as.numeric(parts[, 2])
    min <- as.numeric(parts[, 3])
    sec <- as.numeric(parts[, 4])
    dir <- parts[, 5]
    
    dec <- deg + min / 60 + sec / 3600
    dec[dir %in% c("S", "W")] <- -dec[dir %in% c("S", "W")]
    return(dec)
  }
  
  df %>%
    dplyr::mutate(
      !!lat_col := dms_to_dd(.data[[lat_col]]),
      !!long_col := dms_to_dd(.data[[long_col]])
    )
}

assign_transect_name <- function(df, transect_lat_path) {
  # Read reference data
  transect_ref <- readr::read_csv(transect_lat_path, show_col_types = FALSE)
  
  # Ensure proper types
  transect_ref <- transect_ref %>%
    dplyr::mutate(
      Transect_ID = as.character(Transect_ID),
      Lat_ref = as.numeric(Lat)
    ) %>%
    dplyr::select(Transect_ID, Lat_ref)

  # Assign nearest Transect_ID to each row
  df %>%
    rowwise() %>%
    mutate(
      Transect_ID = if (!is.na(Lat)) {
        nearest_idx <- which.min(abs(Lat - transect_ref$Lat_ref))
        transect_ref$Transect_ID[nearest_idx]
      } else {
        NA_character_
      }
    ) %>%
    ungroup()
}

assign_shelf_location <- function(df, shelf_break_path) {
  # Load the shelf break reference
  shelf_ref <- readr::read_csv(shelf_break_path, show_col_types = FALSE)

  # Ensure correct types and keep only what's needed
  shelf_ref <- shelf_ref %>%
    dplyr::mutate(
      Transect_ID = as.character(Transect_ID),
      Long_shelf_break_200 = as.numeric(Long_shelf_break_200)
    ) %>%
    dplyr::select(Transect_ID, Long_shelf_break_200)

  # Join and assign shelf label
  df %>%
    dplyr::left_join(shelf_ref, by = "Transect_ID") %>%
    dplyr::mutate(
      shelf_loc_200 = dplyr::case_when(
        is.na(Long) ~ NA_character_,
        is.na(Long_shelf_break_200) ~ NA_character_,
        Long < Long_shelf_break_200 ~ "off_shelf",
        Long >= Long_shelf_break_200 ~ "on_shelf"
      )
    ) %>%
    dplyr::select(-Long_shelf_break_200) # optionally remove after use
}

filter_low_class_scores<- function (df, project='MEZCAL'){
	## use results of confusion matrix to remove segments with low likelyhood of correct classification
    if (project=='MEZCAL'){
	probs<-read.csv(Sys.glob(file.path(path.expand("~/SPECTRA_proj/ancillary"), 'MEZCAL_FT_AND_GROUPINGv2.csv')), sep=',')
    	names(probs)[names(probs) == 'class'] <- 'predicted_taxon'
    	df<- left_join(x = df, y = probs, by = "predicted_taxon")
    	df$predicted_taxon <- ifelse(df$probability <= df$filter_prob, "unknown", df$group)
    	df$predicted_taxon<-as.factor(df$predicted_taxon)}
    else if(project=='SPECTRA'){
    	df<-df}
    return(df)}

# assign season label
assign_season <- function(df_clean) {
df_clean <- df_clean %>% mutate(time = ymd_hms(time))
df_clean$MONTH<-month(df_clean$time)
df_clean$DAY<-mday(df_clean$time)
df_clean$DAY_PLUS_MONTH<-df_clean$DAY+df_clean$MONTH*100

df_clean <- df_clean %>%
  mutate(season = case_when(
    DAY_PLUS_MONTH %in% 320:620 ~ 'Spring',
    DAY_PLUS_MONTH %in% 621:922 ~ 'Summer',
    DAY_PLUS_MONTH %in% 923:1221 ~ 'Fall',
    TRUE ~ 'Winter'))

df_clean = subset(df_clean, select = -c(MONTH, DAY, DAY_PLUS_MONTH))

return(df_clean)
}

assign_day_night <- function(df) {
  df <- df %>%
    mutate(
      time_parsed = ymd_hms(time, tz = "UTC"),
      date_local = as.Date(with_tz(time_parsed, tz = "UTC"))
    )
  
  # Create lookup for unique lat/lon/date combinations
  unique_geo_dates <- df %>%
    distinct(date_local, Lat, Long)

  # Get sunrise/sunset for each lat/lon/date
  sun_times <- pmap_dfr(
    list(date = unique_geo_dates$date_local,
         lat = unique_geo_dates$Lat,
         lon = unique_geo_dates$Long),
    function(date, lat, lon) {
      times <- getSunlightTimes(date = date, lat = lat, lon = lon, keep = c("sunrise", "sunset"))
      times %>% mutate(date_local = date, Lat = lat, Long = lon)
    }
  )

  # Join back and assign day/night
  df <- df %>%
    left_join(sun_times, by = c("date_local", "Lat", "Long")) %>%
    mutate(
      day_night = case_when(
        time_parsed >= sunrise & time_parsed < sunset ~ "day",
        TRUE ~ "night"
      )
    ) %>%
    select(-sunrise, -sunset, -date_local, -time_parsed)

  return(df)
}

#gets midpoints of bins obtained by the cut() function
midpoints <- function(x, dp=NULL){
  lower <- as.numeric(gsub(',.*','', gsub('\\(|\\[|\\)|\\]','', x)))# two regex searches: first gsub searches for ( or [ or ) or ) and removes it. Operators: |=or, \=silence (escape)
  upper <- as.numeric(gsub('.*,','', gsub('\\(|\\[|\\)|\\]','', x))) #the other regex search is to find any characters after or before the comma, and remove them 
  
  mid <- lower + ((upper - lower) / 2)
  
  if (is.null(dp)) {
    return(mid)  # full precision
  } else {
    return(round(mid, dp))
  }
}

#gets binwidth obtained by the cut() function
bin_width <- function(x){
  lower <- as.numeric(gsub(',.*','', gsub('\\(|\\[|\\)|\\]','', x)))# two regex searches: first gsub searches for ( or [ or ) or ) and removes it. Operators: |=or, \=silence (escape)
  upper <- as.numeric(gsub('.*,','', gsub('\\(|\\[|\\)|\\]','', x))) #the other regex search is to find any characters after or before the comma, and remove them 
  return(upper-lower)
}

depth_bin_calc <-function(df, depth_var='Depth', depth_bin_width=5, max_depth=100){
  # assing depth bins to data.  initial depth binning will be each 5 m, also remove upper 3m because of bubbles
  df<-df  %>% mutate(Depth= as.numeric(.[[depth_var]])) %>% mutate(Depth = replace(Depth, Depth>100, 100))%>% filter(Depth> 3) %>% mutate(depth_binned = cut(Depth, breaks=seq(0, max_depth, by = depth_bin_width)))
  df$depth_mid<-midpoints(df$depth_binned)
  df$depth_binned<-NULL
  return(df)
}

coord_bin_calc_old <-function(df, lat_var='Lat',lon_var='Long', coord_bin_width_lat=0.25, coord_bin_width_lon=0.0001){
  # assing coordinate bins to data.  initial depth binning will be each 0.25 lat/lon
  df<-df %>% mutate(Lat = as.numeric(.[[lat_var]])) %>%mutate(lat_binned = cut(Lat, breaks=seq(-90, 90, by = coord_bin_width_lat)))
  df<-df %>% mutate(Long = as.numeric(.[[lon_var]])) %>%mutate(lon_binned = cut(Long, breaks=seq(-180, 180, by = coord_bin_width_lon)))
  df$lat_mid<-midpoints(df$lat_binned)
  df$lat_binned<-NULL
  df$lon_mid<-midpoints(df$lon_binned)
  return(df)
}

coord_bin_calc <-function(df, lat_var='Lat',lon_var='Long', lat_breaks, lon_breaks){
  # assing coordinate bins to data.  initial depth binning will be each 0.25 lat/lon
  df<-df %>% mutate(Lat = as.numeric(.[[lat_var]]), Long = as.numeric(.[[lon_var]]),
                    lat_binned =cut(Lat, breaks = lat_breaks), lon_binned =cut(Long, breaks = lon_breaks),
                    lat_mid = midpoints(lat_binned, dp=2), lon_mid = midpoints(lon_binned, dp=3)) %>% 
    select(-lat_binned, -lon_binned)
  
  if (any(df$Long < min(lon_breaks) | df$Long > max(lon_breaks))) {
    warning("Some Longitudes are outside the defined bin range.")
  }
  return(df)
}


# size binning, and after that, we might be able to get particle count PER size PER tax NOTE: modified 12/13/20204, to create adjustable size classes depending on size range
size_binning_adjusted <-function(df, size_var='area',keep_ECD=FALSE,  depth_bin_width=5, ISIIS_ver=2){
  # assign size bins to data. we will use Dugennes file for size classes(same as ecopart size bins). Some relevant conversions:
  #ISIIS-2 pixel size= 66um, unclear what the units of 'area' are, pixels squared?
  #ISIIS-2 volume per frame: 13*13*50=8450 cm3, or 8.45 L
  #ISIIS-3 pixel size=55um, definitely in pixels squared
  # detected some NANs in MEZCAL data, they were all detritus_not_pyrosome. Question: should we remove detritus? At least we should remove NAN
   if (ISIIS_ver==2){
     pixel_size=66
   } else if (ISIIS_ver==3) {
     pixel_size=55
   }
  df[[size_var]]<-as.numeric(df[[size_var]])
  df <-df[!is.na(df[[size_var]]),]#check about this
  df <-df[df[[size_var]]!=0,]
  rownames(df) <- NULL # to reset row indices
  df$biovol<-(4/3)* pi * ((((df[[size_var]]*(pixel_size**2))/pi)**0.5) **3)
  #next two lines added 12/13/2024
  lseq <- function(from=min(df$biovol), to=max(df$biovol), length.out=10) {
        exp(seq(log10(from), log10(to), length.out = length.out))}
  breaks <- lseq()
  df$ECD<-2*sqrt((df[[size_var]]*(pixel_size**2))/pi)
  df$biovol<-(1./6.)*pi*(df$ECD**3)
  df<-df %>% mutate(biovol_binned = cut(biovol, breaks=breaks))# added 12/13/2024, notice that breaks come from the function four lines above
  df$biovol_um3_mid<-midpoints(df$biovol_binned)
  df$biovol_um3_width<-bin_width(df$biovol_binned)
  if (keep_ECD==FALSE){
  df$ECD<-NULL}
  return(df)
}


# size binning, and after that, we might be able to get particle count PER size PER tax
size_binning <-function(df, size_var='area',keep_ECD=FALSE,  depth_bin_width=5, ISIIS_ver=2){
  # assign size bins to data. we will use Dugennes file for size classes(same as ecopart size bins). Some relevant conversions:
  #ISIIS-2 pixel size= 66um, unclear what the units of 'area' are, pixels squared?
  #ISIIS-2 volume per frame: 13*13*50=8450 cm3, or 8.45 L
  #ISIIS-3 pixel size=55um, definitely in pixels squared
  # detected some NANs in MEZCAL data, they were all detritus_not_pyrosome. Question: should we remove detritus? At least we should remove NAN
   if (ISIIS_ver==2){
     pixel_size=66
   } else if (ISIIS_ver==3) {
     pixel_size=55
   }
  df[[size_var]]<-as.numeric(df[[size_var]])
  df <-df[!is.na(df[[size_var]]),]#check about this
  rownames(df) <- NULL # to reset row indices
  size_classes<-read.csv(Sys.glob(file.path(path.expand("~/SPECTRA_proj/ancillary"), 'ecopart_size_bins.csv')), sep=',')
  df$ECD<-2*sqrt((df[[size_var]]*(pixel_size**2))/pi)
  df$biovol<-(1./6.)*pi*(df$ECD**3)
  df<-df %>% mutate(biovol_binned = cut(biovol, breaks=size_classes$biovol_um3))
  df$biovol_um3_mid<-midpoints(df$biovol_binned)
  df$biovol_um3_width<-bin_width(df$biovol_binned)
  if (keep_ECD==FALSE){
  df$ECD<-NULL}
  return(df)
}

group_abundance_func<- function(df, grouping_vars=c('season', 'depth_mid', 'lat_mid', 'lon_mid'), env_vars=c('Temp', 'Salinity', 'Fluoro', 'Oxygen', 'Irradiance'), group_by_taxa =TRUE, ISIIS_ver=2, thresholding = TRUE){
	#this function will use any approach to define a plankton community, and assign a unique identifier to each community (Group_ID) and then will estimate the abundance (biovol/L^1) of each size class
   if (ISIIS_ver==2){
     pixel_size<-66
     vol_frame<-8.45
   } else if (ISIIS_ver==3) {
     pixel_size<-55
     vol_frame<- 8.45 # modify this
   }
    # Step 1: Frame separation and initial volume estimate
    df<-df %>% separate(img, into=c('frame_id', 'img'), sep='_crop_', remove=TRUE) %>% group_by(!!!syms(grouping_vars)) %>% mutate(n_frames=n_distinct(frame_id)) %>%as.data.frame() # volume sampled needs to be based on environmental variables only, not on taxa, else the volume sampled would be too low for rare taxa
    df$vol_sampled<- df$n_frames*vol_frame

    # Step 2: Optional taxonomic grouping
    if (group_by_taxa==TRUE){
       grouping_vars<- c(grouping_vars,'predicted_taxon')# generate size spectra per taxa
       df$predicted_taxon<-vapply(strsplit(df$predicted_taxon, "_"), '[', 1, FUN.VALUE=character(1))} #edit this to provide coarse or fine taxonomic resolution
  
  
    # Step 3: Assign group ID and average environmental vars across group
    df <- df %>%
      group_by(!!!syms(grouping_vars)) %>%
      mutate(group_id = cur_group_id()) %>%
      mutate(across(all_of(env_vars), ~mean(.x, na.rm = TRUE)))%>%
      ungroup() %>%
      as.data.frame()
 
    df <- df[order(df$group_id),]
    
    # Step 4.1 Reload the ecotaxa size bins to assign empty size bins
    size_classes <- read.csv(Sys.glob(file.path(path.expand("~/SPECTRA_proj/ancillary"), 'ecopart_size_bins.csv')), sep = ',')
    all_mids <- as.numeric(as.character(midpoints(cut(as.numeric(size_classes$biovol_um3), breaks = as.numeric(size_classes$biovol_um3)))))

    # Step 4.2: Summarize into size classes, keeping env vars (assumed constant per group)
    df <- df %>%
      group_by(!!!syms(c('group_id', grouping_vars, 'biovol_um3_width', 'biovol_um3_mid'))) %>%
      summarise(
        n_individuals = n(),
        total_biovolume = sum(biovol),
        vol_sampled = first(vol_sampled),
        across(all_of(env_vars), ~first(.x, na.rm = TRUE))
      ) %>%
      as.data.frame()
    
    ## Step 4.3: organize size bins across all ecotaxa size classes, to determine empty size bins 
    df$biovol_um3_mid <- as.numeric(df$biovol_um3_mid)
    df <- df %>% 
      complete( nesting(group_id, !!!syms(grouping_vars), vol_sampled, !!!syms(env_vars)),
                        biovol_um3_mid = all_mids,
                        fill = list(n_individuals = 0, total_biovolume = 0))

    # Step 5: Calculate normalized biovolume spectra and uncertainty metrics
    df <- df %>%
        mutate(
            NB = ifelse(is.na(n_individuals) | n_individuals == 0 | is.na(biovol_um3_width), NaN, (total_biovolume / vol_sampled) / biovol_um3_width),
            ind_L = ifelse(is.na(n_individuals) | n_individuals == 0 | is.na(biovol_um3_width), NaN, n_individuals / vol_sampled),
            logNB = ifelse(is.na(n_individuals) | n_individuals == 0 | is.na(biovol_um3_width), NaN, log10((total_biovolume / vol_sampled) / biovol_um3_width)),
            logBiovol = log10(biovol_um3_mid),
            ECD_um_mid = ((biovol_um3_mid * 6) / pi)^(1./3.),
            ECD_um_width = ((biovol_um3_width * 6) / pi)^(1/3), 
            pixel_size_mid = (ECD_um_mid / 1000) * (1 / (pixel_size / 1000)),
            size_uncertainty = ifelse(is.na(n_individuals) | n_individuals == 0 | is.na(biovol_um3_width), NaN, dnorm((1/3)*pi*pixel_size_mid^3, mean = (1/3)*pi*60^3, sd = (1/3)*1*pi)),
            count_uncertainty = ifelse(is.na(n_individuals), NaN, dpois(x = n_individuals, lambda = n_individuals)) 
    )

    #Flags to keep track of empty size bins
    invalid_idx <- which(!is.finite(df$NB) & df$n_individuals > 0)
    if (length(invalid_idx) > 0) {
      message("⚠️  Found ", length(invalid_idx), " non-finite NB values")
      for (i in invalid_idx) {
        message("🚨 group_id: ", df$group_id[i])
        message("  total_biovolume = ", df$total_biovolume[i])
        message("  vol_sampled     = ", df$vol_sampled[i])
        message("  biovol_width    = ", df$biovol_um3_width[i])
        message("  NB              = ", df$NB[i])
      }
    }
    
    # old code included in the mutate command in line 281
    if (thresholding ==TRUE){
        df<-threshold_func(df)
    }
    return(df)
}
## integrate these two lines below to the function above:
# for poisson: https://makemeanalyst.com/statistics-with-r/poisson-functions-in-r-programming/#:~:text=The%20probability%20mass%20function%20(PMF,is%20the%20factorial%20of%20k.
# for gaussian: https://r02pro.github.io/normal-distribution.html, Mathilde did it in pixels per mm, but the size_class_pixel is each size class in pixels

threshold_func<-function(df, threshold_size=0.2, threshold_count=0.2){
 n_before <- nrow(df)

  df <- df[which(df$size_uncertainty <= threshold_size), ]
  df <- df[which(df$count_uncertainty <= threshold_count), ]
  df <- df[which(df$pixel_size_mid >= 20), ]

  n_after <- nrow(df)
  message("📉 Thresholding removed ", n_before - n_after, " rows (", n_after, " remaining)")

  return(df)
# also include this on the function, which is cutting the size spectra based on empty size bins, check https://masterr.org/r/how-to-find-consecutive-repeats-in-r/


}


group_id_func<-function(df, grouping_vars=c('season', 'depth_mid', 'lat_mid', 'lon_mid'), group_by_taxa =TRUE){
   if (group_by_taxa==TRUE){
      grouping_vars<- c(grouping_vars,'predicted_taxon')}
   df<-df %>% group_by(!!!syms(grouping_vars)) %>% mutate(group_id = cur_group_id())%>%as.data.frame() #assign group_id
   df <- df[order(df$group_id),]
    }
    
lin_fit_func<- function(df, filter_min_size_class=FALSE, n_min_size_class=3){  
    fitted_models <- df %>% filter(is.finite(logNB), is.finite(logBiovol)) %>% group_by(group_id) %>% do(model = lm(logNB ~ logBiovol, data = .))
    ## this is the structure of how the data is obstained: fitted_models[[2]][['group_id']]$coefficients[1] for intercept, coefficients[2] for slope
    bad <- df %>% filter(!is.finite(logNB) | !is.finite(logBiovol))
    if (nrow(bad) > 0) {
        message("⚠️ Dropping ", nrow(bad), " rows with non-finite logNB or logBiovol before fitting")
  
        bad_summary <- bad %>% select(group_id, logNB, logBiovol, NB, total_biovolume, vol_sampled, biovol_um3_width, biovol_um3_mid) %>% distinct() %>% arrange(group_id)
  
        message("🔎 Affected group summary:\n", 
          capture.output(print(bad_summary)) %>% paste(collapse = "\n"))
    }
    valid_groups <- fitted_models$group_id
    model_list <- fitted_models$model
    
    for (i in seq_along(valid_groups)) {
      g <- valid_groups[i]
      model_entry <- model_list[[i]]
    
      df$intercept[df$group_id == g] <- model_entry$coefficients[1]
      df$slope[df$group_id == g]     <- model_entry$coefficients[2]
      df$R2[df$group_id == g]        <- summary(model_entry)$r.squared
      df$n_size_classes[df$group_id == g] <- nrow(model_entry$model)
    }
    
    if (filter_min_size_class==TRUE){
       df <- df[order(df$group_id),]  %>% filter(n_size_classes> n_min_size_class)}
    return(df)
}

# Collapse per-taxon columns into broad groups between two named columns.
# Example:
collapse_taxa_cols <- function(
    df,
    start_name,
    end_name,
    patterns = list(
      copepod              = "copepod",
      crustacean           = "crustacean",
      predatory_gelatinous = "(hydromedusae|chaetognath|siphonophore|ctenophore)",
      grazer_gelatinous    = "(appendicularian|tunicate|pyrosome|pteropod)",
      fish                 = "fish",
      detritus             = "detritus"
    )
) {
  stopifnot(is.data.frame(df))
  clean <- function(x) tolower(gsub("\\s+", "_", x))
  nms        <- names(df)
  nms_clean  <- clean(nms)
  
  s_idx <- which(nms_clean == clean(start_name))
  e_idx <- which(nms_clean == clean(end_name))
  if (length(s_idx) != 1 || length(e_idx) != 1) {
    stop(sprintf("Could not uniquely locate start/end columns: '%s' / '%s'", start_name, end_name))
  }
  
  rng <- if (s_idx <= e_idx) s_idx:e_idx else e_idx:s_idx
  taxa_cols       <- nms[rng]
  taxa_cols_clean <- nms_clean[rng]
  
  # --- FIX: start with a 0-column data.frame that has the right number of rows
  collapsed <- df[, 0, drop = FALSE]
  
  rs <- function(cols) {
    if (length(cols) == 0) return(rep(0, nrow(df)))
    rowSums(df[, cols, drop = FALSE], na.rm = TRUE)
  }
  
  idx <- lapply(patterns, function(p) grepl(p, taxa_cols_clean, ignore.case = TRUE))
  detritus_idx  <- idx$detritus %||% rep(FALSE, length(taxa_cols))
  detritus_cols <- taxa_cols[detritus_idx]
  
  group_names <- setdiff(names(patterns), "detritus")
  grp_cols <- lapply(group_names, function(g) taxa_cols[idx[[g]]])
  names(grp_cols) <- group_names
  
  # Add grouped columns
  for (g in group_names) collapsed[[g]] <- rs(grp_cols[[g]])
  
  matched_non_det <- unique(unlist(grp_cols, use.names = FALSE))
  other_cols <- setdiff(setdiff(taxa_cols, matched_non_det), detritus_cols)
  collapsed$other    <- rs(other_cols)
  collapsed$detritus <- rs(detritus_cols)
  
  keep_cols <- setdiff(nms, taxa_cols)
  out <- cbind(df[, keep_cols, drop = FALSE], collapsed)
  rownames(out) <- NULL
  
  message(sprintf(
    "Collapsed %d taxa columns (%s … %s) into %d groups (+other, detritus).",
    length(taxa_cols), taxa_cols[1], taxa_cols[length(taxa_cols)], length(group_names)
  ))
  out
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
