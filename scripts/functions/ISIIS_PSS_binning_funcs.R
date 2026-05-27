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
library(slider)
library(rlang)
library(broom)

biovolume_integral <- function(m, b, x1, x2, tol = 1e-12) {
  # Area under y = b * x^m over [x1, x2]
  if (isTRUE(all.equal(m, -1, tolerance = tol))) {
    # limit as m -> -1
    return(b * log(x2 / x1))
  } else {
    return((b / (m + 1)) * (x2^(m + 1) - x1^(m + 1)))
  }
}



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

#convert degrees minutes seconds to degrees decimal
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

assign_transect_name <- function(df, transect_lat_path, lat_col_name) {
  # Read reference data
  transect_ref <- read_csv(transect_lat_path, show_col_types = FALSE) %>%
    mutate(
      Transect_ID = as.character(Transect_ID),
      Lat_ref = as.numeric(Lat)
    ) %>%
    select(Transect_ID, Lat_ref)

  # Use tidy eval with .data[[ ]] to dynamically reference the latitude column
  df %>%
    rowwise() %>%
    mutate(
      Transect_ID = {
        lat_value <-as.numeric(.data[[lat_col_name]])
        if (!is.na(lat_value)) {
          nearest_idx <- which.min(abs(lat_value - transect_ref$Lat_ref))
          transect_ref$Transect_ID[nearest_idx]
        } else {
          NA_character_
        }
      }
    ) %>%
    ungroup()
}


assign_shelf_location <- function(df, shelf_break_path, lon_col_name='Longitude') {
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
        is.na({{lon_col_name}}) ~ NA_character_,
        is.na(Long_shelf_break_200) ~ NA_character_,
        {{lon_col_name}} < Long_shelf_break_200 ~ "off_shelf",
        {{lon_col_name}} >= Long_shelf_break_200 ~ "on_shelf"
      )
    ) %>%
    dplyr::select(-Long_shelf_break_200) # optionally remove after use
}


# add distance from shore (km) and relative distance
assign_relative_distance_old <-function(df, land_ref_path, lat_col= 'Latitude', long_col='Longitude'){
    lon <- ensym(long_col)
    lat <- ensym(lat_col)  
    land_ref <- readr::read_csv(land_ref_path, show_col_types = FALSE)
    max_abs_long_per_transect <- df %>%
      group_by(Transect_ID) %>%
      summarise(max_abs_Long = -1 * max(abs(!!lon), na.rm = TRUE))

    land_ref_expanded <- max_abs_long_per_transect %>%
      mutate(Transect_base = str_remove(Transect_ID, "_[A-Z]+$")) %>%
      left_join(land_ref, by = c("Transect_base" = "Transect_ID")) %>%
      select(Transect_ID, mean_lat, shelf_long, land_long, max_abs_Long)

      library(geosphere)
      library(dplyr)

    add_shelf_dists <- function(df) {
      df %>%
        rowwise() %>%
        mutate(
          # Distance between shelf break and land (same lat)
          dist_shelf_break_km = distGeo(
            c(shelf_long, mean_lat),  # lon, lat
            c(land_long, mean_lat)
          ) / 1000,
      
          # Distance between farthest point (abs max) and land (same lat)
          dist_total_km = distGeo(
            c(max_abs_Long, mean_lat),
            c(land_long, mean_lat)
          ) / 1000
        ) %>%
        ungroup()
    }
    land_ref_with_dists <- add_shelf_dists(land_ref_expanded)

    df_joined <- df %>% left_join(land_ref_with_dists, by = "Transect_ID")

    add_relative_distance <- function(df) {
      df %>%
        rowwise() %>%
        mutate(
          distance_sample_km = distGeo(
            c(!!lon, !!lat),
            c(land_long, !!lat)
          ) / 1000,
          rel_dist = distance_sample_km / dist_shelf_break_km
        ) %>%
        ungroup()
    }

    df_joined_dist <- add_relative_distance(df_joined)
    return(df_joined_dist)
}

assign_relative_distance <- function(df, land_ref_path, lat_col = "Latitude", long_col = "Longitude"){
  # columns to use
  lon <- rlang::ensym(long_col)
  lat <- rlang::ensym(lat_col)

  # ref table: Transect_ID, mean_lat, shelf_long, land_long
  land_ref <- readr::read_csv(land_ref_path, show_col_types = FALSE)

  # westernmost long per transect (if you still want it)
  max_abs_long_per_transect <- df %>%
    dplyr::group_by(Transect_ID) %>%
    dplyr::summarise(max_abs_Long = -1 * max(abs(!!lon), na.rm = TRUE), .groups = "drop")

  land_ref_expanded <- max_abs_long_per_transect %>%
    dplyr::left_join(land_ref, by = "Transect_ID") %>%
    dplyr::select(Transect_ID, mean_lat, shelf_long, land_long, max_abs_Long)

  # join refs to each row
  df_joined <- df %>% dplyr::left_join(land_ref_expanded, by = "Transect_ID")

  # distances (vectorized; no rowwise)
  # 1) shelf→land distance at transect mean latitude
  dist_shelf_break_km <- geosphere::distGeo(
    cbind(df_joined$shelf_long, df_joined$mean_lat),
    cbind(df_joined$land_long,  df_joined$mean_lat)
  ) / 1000

  # 2) signed distance from SHELF to the sample at the sample’s latitude
  #    sign = + if sample is east of shelf (toward land), - if west (offshore)
  dist_from_shelf_km <- geosphere::distGeo(
    cbind(df_joined$shelf_long, df_joined[[rlang::as_string(lat)]]),
    cbind(df_joined[[rlang::as_string(lon)]], df_joined[[rlang::as_string(lat)]])
  ) / 1000

  sign_dir <- sign(df_joined[[rlang::as_string(lon)]] - df_joined$shelf_long)

  df_joined %>%
    dplyr::mutate(
      dist_shelf_break_km = dist_shelf_break_km,
      distance_sample_km  = dist_from_shelf_km,
      rel_dist            = sign_dir * (distance_sample_km / dist_shelf_break_km)
      # rel_dist: 0 at shelf_long; +1 near land_long; negative offshore
    )
}

filter_low_class_scores<- function (df,base_dir,project='MEZCAL'){
	## use results of confusion matrix to remove segments with low likelyhood of correct classification
    base_dir<- getwd()
    if (project=='MEZCAL'){
	probs<-read.csv(Sys.glob(file.path(path.expand(base_dir), 'ancillary', 'MEZCAL_FT_AND_GROUPINGv2.csv')), sep=',')
    }else if(project=='SPECTRA'){
	probs<-read.csv(Sys.glob(file.path(path.expand(base_dir), 'ancillary', 'SPECTRA_90_prob_filtering_manyDetritus_python.csv')), sep=',')
    }
    	names(probs)[names(probs) == 'class'] <- 'predicted_taxon'
    	df<- left_join(x = df, y = probs, by = "predicted_taxon")
    	df$predicted_taxon <- ifelse(df$probability <= df$filter_prob, "unknown", df$group)
    	df$predicted_taxon<-as.factor(df$predicted_taxon)
    return(df)}

apply_correction_factors <- function(df, base_dir, project='SPECTRA', binning='abundance'){
    base_dir<- getwd()

     if (project=='MEZCAL'){
        probs<-read.csv(file.path(path.expand(base_dir), 'ancillary', 'MEZCAL_FT_AND_GROUPINGv2.csv'), sep=',')
    }else if(project=='SPECTRA'){
        probs<-read.csv(file.path(path.expand(base_dir), 'ancillary', 'classification_metrics_90.0_manyDetritus_python.csv'), sep=',')
    }   
    names(probs)[names(probs) == 'Class'] <- 'group'
    probs$correction_factor <- probs$Precision/probs$Recall
    df<- df %>% left_join(probs %>% select(group, correction_factor), by = "group")
    
    if (binning == 'abundance'){
	df$ind_m3 <- df$ind_m3*df$correction_factor
    }else if (binning == 'size_spectra'){
	df$total_biovolume <- df$total_biovolume*df$correction_factor
	df$n_individuals <- df$n_individuals*df$correction_factor
    }
    return(df)
}

apply_correction_factors_wide <- function(df, base_dir, project = 'SPECTRA', binning = 'abundance') {
  # === Load correction table ===
  probs <- switch(
    project,
    'MEZCAL' = read.csv(file.path(base_dir, 'ancillary', 'MEZCAL_FT_AND_GROUPINGv2.csv')),
    'SPECTRA' = read.csv(file.path(base_dir, 'ancillary', 'classification_metrics_90.0_manyDetritus_python.csv')),
    stop("Unknown project")
  )
  
  # Prep correction data
  probs <- probs %>%
    dplyr::rename(group = Class) %>%
    dplyr::mutate(correction_factor = Precision / Recall)

  # === Apply corrections based on binning type ===
  if (binning == "abundance") {
    # Find all ind_m3_* columns
    ind_cols <- grep("^ind_m3_", names(df), value = TRUE)
    
    for (col in ind_cols) {
      group_name <- sub("^ind_m3_", "", col)
      cf <- probs$correction_factor[probs$group == group_name]
      
      if (length(cf) == 1 && !is.na(cf)) {
        message("Applying correction to: ", col)
        df[[col]] <- df[[col]] * cf
      } else {
        warning("No correction factor found for: ", group_name)
      }
    }

  } else if (binning == "size_spectra") {
    # Placeholder for future logic if biovolume columns follow similar naming
    stop("Size spectra correction not implemented for wide format")
  } else {
    stop("Unknown binning type")
  }

  return(df)
}

# assign season label
assign_season <- function(df_clean, time_col = "time", tz= "America/Los_Angeles") {
  df_clean <- df_clean %>%
    mutate(time = ymd_hms(.data[[time_col]], tz = tz))  # dynamic col name + tz

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

library(dplyr)
library(lubridate)
library(suncalc)
library(purrr)

assign_day_night <- function(df, crepuscular = FALSE, time_col= "time_Pacific", lat_col= "Latitude", long_col= "Longitude", tz = "America/Los_Angeles") {
  df <- df %>%
    mutate(
      time_parsed = ymd_hms(.data[[time_col]], tz = tz),
      date_local = as.Date(with_tz(time_parsed, tz = tz))
    )

  # Create lookup for unique lat/lon/date combinations
  unique_geo_dates <- df %>%
    distinct(date_local, .data[[lat_col]], .data[[long_col]])

  # Get solar times
  keep_times <- if (crepuscular) c("sunrise", "sunset", "dawn", "dusk") else c("sunrise", "sunset")
  sun_times <- pmap_dfr(
    list(date = unique_geo_dates$date_local,
         lat = unique_geo_dates[[lat_col]],
         lon = unique_geo_dates[[long_col]]),
    function(date, lat, lon) {
      times <- getSunlightTimes(date = date, lat = lat, lon = lon, keep = keep_times)
      times %>% mutate(date_local = date, !!sym(lat_col) := lat, !!sym(long_col) := lon)
    }
  )

  # Join and assign time of day
  if (crepuscular){
     df <- df %>%
       left_join(sun_times, by = c("date_local", lat_col, long_col)) %>%
       mutate(
         time_of_day = case_when(
           time_parsed > dusk ~ "night",
           time_parsed < dawn ~ "night",
           time_parsed > sunrise & time_parsed < sunset ~ "day",
           time_parsed > dawn & time_parsed <= sunrise ~ "crepuscular",
           time_parsed >= sunset & time_parsed < dusk ~ "crepuscular",
           TRUE ~ "night"
         )
       ) 
  } else {
    df <- df %>%
      left_join(sun_times, by = c("date_local", lat_col, long_col)) %>%
      mutate(time_of_day = case_when(time_parsed >= sunrise & time_parsed < sunset ~ "day", TRUE ~ "night"))
  }

 df<- df %>% select(-any_of(c("sunrise", "sunset", "dawn", "dusk")), -date_local, -time_parsed) 

  return(df)
}

# Sensor Calibration Functions

library(dplyr)
apply_calibration_factors <- function(df){
  df$Temp_Kelvin <- df$Temperature + 273.15 
  Voffset <-  -0.4978
  soc <- 0.4540
  A <- -0.0031257
  B <- 0.00012009
  C <- -1.6988e-6
  E <- 0.036
  A0 <- 2.00907
  A1 <- 3.22014
  A2 <- 4.0501
  A3 <- 4.94457
  A4 <- -0.256847
  A5 <- 3.88767
  B0 <- -0.00624523
  B1 <- -0.00737614
  B2 <- -0.010341
  B3 <- -0.00817083
  C0 <- -4.88682e-7
  Ts <-  log(((298.15 - df$Temperature) / (273.15 + df$Temperature) ) )
  OxSol <- exp(A0+ A1*Ts + A2*(Ts^2) + A3*(Ts^3) + A4*(Ts^4) + A5*(Ts^5) + df$Salinity*(B0 +  B1*Ts + B2*(Ts^2) + B3*(Ts^3) + C0*(df$Salinity)^2 ))
  df$Oxygen_ml_l <- soc*(df$`Raw Oxygen [V]`+Voffset)*OxSol*(1+A*df$Temperature)+(B*(df$Temperature^2))+C*(df$Temperature^3)+exp(E*df$Pressure/df$Temp_Kelvin)
  
  df$pH <- 7 + (df$`Raw pH [V]` - 2.5553)/(1.98416e-4*(df$Temp_Kelvin)*4.5297)
  
  df$PAR <- (5.22e-6)*(10^df$`Raw PAR [V]`- 10^0.0119)
  
  df$chl_a_ug_l <- 0.0072*(df$Chlorophyll - 49)
    
  df <- df %>%
    dplyr::select(-any_of(c("Raw Oxygen [V]", "Raw PAR [V]", "Raw pH [V]", "Temp_Kelvin", "Chlorophyll"))) %>%
    relocate(c(Oxygen_ml_l,pH, PAR), .before = "Latitude" )
  message("Calibration factors applied to oxygen, PAR, pH, and chlorophyll.")
  return(df)
}

# interpolate missing data + look for erroneous data (negatives)
interpolate_missing_values <- function(df, cols, time_col='time_Pacific'){
  for (col in cols){
  message("checking ", col)
  # For non-coordinate columns, replace negatives with NA
  #if (!col %in% c("Latitude", "Longitude")) {
  if (!grepl("Latitude|Longitude", col)) {
     df[[col]] <- ifelse(df[[col]] < 0, NA, df[[col]])
  }
  if (all(!is.na(df[[col]]))){
    message("No negative or NA values. No interpolation needed.") 
    next 
  } else {
    message("Number of negative or missing values for ", col, ":", sum(is.na(df[[col]])))
    # preserve the original data - we only want to replace the NAs 
    original <- df[[col]]
    # Check if there are at least two unique x values
    if (length(unique(df[[time_col]][!is.na(original)])) < 2) {
      message("Not enough unique timestamps to interpolate for ", col, ". Skipping.")
      next
    }
    # interpolate the missing or negative data based on the time stamp and other values around the missing data 
    interpolated <- approx(
      x = df[[time_col]],
      y = original,
      xout = df[[time_col]],
      method = "linear",
      rule = 2
    )$y
  
    # replace only the missing values with the interpolated values 
    df[[col]][is.na(original)] <- interpolated[is.na(original)]
    message("Done: NAs and negative values interpolated")
    }  
  }
 return(df)
}

despike_depth_1 <- function(df, depth_bin_col = "depth_mid", time_col = "time_Pacific") {
  # Ensure the specified columns exist
  if (!(depth_bin_col %in% names(df)) || !(time_col %in% names(df))) {
    stop("Specified depth_bin_col or time_col not found in the data frame.")
  }

  # Arrange by time
  df <- df %>% arrange(.data[[time_col]])

  # Pull depth vector
  depth_vec <- df[[depth_bin_col]]

  # Identify spikes: depth[i] is different from depth[i-1] and depth[i+1], but depth[i-1] == depth[i+1]
  is_spike <- c(FALSE,
                (depth_vec[-c(1, length(depth_vec))] != depth_vec[-c(1, 2)] &
                 depth_vec[-c(1, length(depth_vec))] != depth_vec[-c(3, length(depth_vec))] &
                 depth_vec[-c(1, length(depth_vec))] == depth_vec[-c(3, length(depth_vec))]),
                FALSE)

  # Replace spikes with surrounding value
  corrected <- depth_vec
  corrected[is_spike] <- depth_vec[which(is_spike) - 1]

  # Replace the column in the data frame
  df[[depth_bin_col]] <- corrected

  return(df)
}



despike_depth_2 <- function(df, depth_bin_col = "depth_mid", time_col = "time_Pacific") {
   # Column exists?
  if (!all(c(depth_bin_col, time_col) %in% names(df))) return(df)
  
  df <- df[order(df[[time_col]]), ]
  depth <- df[[depth_bin_col]]
  

  # Early exit if only one unique depth value
  if (all(is.na(depth))) return(df)
  if (length(na.omit(depth)) <= 1) return(df)

  df <- df[order(df[[time_col]]), ]
  depth <- df[[depth_bin_col]]
  
  r <- rle(depth)
  values <- r$values
  lengths <- r$lengths
  
 # If fewer than 3 runs, no A-B-A or zigzag can exist
  if (length(values) < 3) return(df)

  for (i in 2:(length(values)-1)) {
    # Case 1: single-value spike A-B-A
    if (lengths[i] == 1 && values[i-1] == values[i+1] && values[i] != values[i-1]) {
      values[i] <- values[i-1]
    }

    # Case 2: short alternating zigzag (e.g. A-B-A-B)
    if (i > 2 && i < length(values) - 2) {
      pattern <- values[(i-2):(i+2)]
      if (length(unique(pattern)) == 2 && all(lengths[(i-2):(i+2)] == 1)) {
        # Replace with the prior stable value
        values[i] <- values[i-1]
      }
    }
  }

  df[[depth_bin_col]] <- inverse.rle(list(values = values, lengths = lengths))
  return(df)
}

#this is a neato loop that iterates over the depsike_depth function until it doesnt change
repeat_despike_depth <- function(df, depth_bin_col = "depth_mid", time_col = "time_Pacific", max_iter = 10) {
  for (i in seq_len(max_iter)) {
    before <- df[[depth_bin_col]]
    df <- despike_depth_2(df, depth_bin_col = depth_bin_col, time_col = time_col)
    after <- df[[depth_bin_col]]
    if (identical(before, after)) break
  }
  return(df)
}

read_all_taxa_cols<- function(path, all_cols) {
  df <- vroom(path, show_col_types = FALSE)
  missing <- setdiff(all_cols, names(df))
  if (length(missing) > 0) {
    df[missing] <- 0  # fill with zeroes
  }
  # Reorder columns to match the union
  df <- df[, all_cols]
  return(df)
}

collapse_group_columns <- function(df, target_string) {
  # Get matching columns
  matched_cols <- grep(target_string, names(df), value = TRUE)
  
  if (length(matched_cols) == 0) {
    warning("No columns matched: ", target_string)
    return(df)
  }

  # Helper to collapse by pattern
  collapse_by_type <- function(pattern, fun) {
    cols <- matched_cols[grepl(pattern, matched_cols)]
    if (length(cols) == 0) return(NULL)
    result <- df %>%
      dplyr::select(all_of(cols)) %>%
      purrr::pmap_dbl(~ fun(c(...), na.rm = TRUE))
    return(result)
  }

  # Collapse each metric
  out <- tibble(
    !!paste0("mean_dark_", target_string)              := collapse_by_type("mean_dark", mean),
    !!paste0("mean_contour_", target_string)           := collapse_by_type("mean_contour", mean),
    !!paste0("mean_minor_ellipse_axis_", target_string):= collapse_by_type("mean_minor_ellipse", mean),
    !!paste0("mean_major_ellipse_axis_", target_string):= collapse_by_type("mean_major_ellipse", mean),
    !!paste0("mean_ellipse_area_", target_string)      := collapse_by_type("mean_ellipse_area", mean),
    !!paste0("n_obs_", target_string)                  := collapse_by_type("n_obs", sum)
  )

  # Drop old columns and bind result
  df_cleaned <- df %>% dplyr::select(-all_of(matched_cols))
  bind_cols(df_cleaned, out)
}

assign_cast_id_simple <- function(depth_vec) {
  direction <- sign(diff(depth_vec))
  change_points <- c(TRUE, diff(direction) != 0)
  cast_ids <- cumsum(change_points)
  cast_ids <- c(cast_ids[1], cast_ids)  # Prepend first group to match length
  cast_ids[seq_along(depth_vec)]        # Ensure it's the right length
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

assign_cast_id <- function(df,
                           depth_col = "depth_mid",
                           time_col  = "mean_time_Pacific",
                           tol = 0) {             # set tol > 0 to ignore tiny jitters
  depth <- rlang::ensym(depth_col)
  time  <- rlang::ensym(time_col)

  out <- df %>%
    arrange(!!time) %>%
    mutate(
      delta = !!depth - lag(!!depth),
      trend = case_when(
        delta >  tol ~  1L,   # downcast (deeper)
        delta < -tol ~ -1L,   # upcast (shallower)
        TRUE           ~ NA_integer_
      )
    ) %>%
    fill(trend, .direction = "down") %>%
    mutate(
      change   = trend != lag(trend),
      cast_id  = 1L + cumsum(coalesce(change, FALSE)),   # numeric sequence
      cast_id  = paste0("cast_", cast_id)                # label it
    ) %>%
    select(-delta, -trend, -change) %>%
    relocate(cast_id, .before = 1)                       # make first column

  out
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
  #size_classes<-read.csv(Sys.glob(file.path(path.expand("~/SPECTRA_proj/ancillary"), 'ecopart_size_bins.csv')), sep=',') commented out 12/13/2024
  df$biovol<-(4/3)* pi * ((((df[[size_var]]*(pixel_size**2))/pi)**0.5) **3)
  #next two lines added 12/13/2024
  lseq <- function(from=min(df$biovol), to=max(df$biovol), length.out=10) {
        exp(seq(log10(from), log10(to), length.out = length.out))}
  breaks <- lseq()
  df$ECD<-2*sqrt((df[[size_var]]*(pixel_size**2))/pi)
  df$biovol<-(1./6.)*pi*(df$ECD**3)
  #df<-df %>% mutate(biovol_binned = cut(biovol, breaks=size_classes$biovol_um3))# commented out 12/13/20204
  df<-df %>% mutate(biovol_binned = cut(biovol, breaks=breaks))# added 12/13/2024, notice that breaks come from the function four lines above
  df$biovol_um3_mid<-midpoints(df$biovol_binned)
  df$biovol_um3_width<-bin_width(df$biovol_binned)
  #df$biovol_binned<-NULL
  if (keep_ECD==FALSE){
  df$ECD<-NULL}
  return(df)
}


# size binning, and after that, we might be able to get particle count PER size PER tax
size_binning <-function(df, size_var='area',keep_ECD=FALSE,  depth_bin_width=5, ISIIS_ver=2, work_dir=base_dir){
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
  size_classes<-read.csv(file.path(work_dir, 'ancillary', 'contour_custom_ISIIS-3_size_classes.csv'), sep=',')
  #df$biovol<-(4/3)* pi * ((((df[[size_var]]*(pixel_size**2))/pi)**0.5) **3) 
  #build breaks from lower/upper columns
  breaks <- c(size_classes$biovol_lower[1], size_classes$biovol_upper)

  df$ECD<-2*sqrt((df[[size_var]]*(pixel_size**2))/pi)
  df$biovol<-(1./6.)*pi*(df$ECD**3)
  df<-df %>% mutate(biovol_binned = cut(df$biovol, breaks = breaks, right = FALSE, include.lowest = TRUE))
  
  #index into size_classes by bin id
  bin_id <- as.integer(df$biovol_binned)   # 1..n or NA for out-of-range
  ok <- !is.na(bin_id)
  df$biovol_um3_mid   <- NA_real_
  df$biovol_um3_width <- NA_real_
  df$biovol_um3_mid[ok]   <- size_classes$biovol_um3_mid[bin_id[ok]]
  df$biovol_um3_width[ok] <- size_classes$biovol_um3_width[bin_id[ok]]
  
  #df$biovol_binned<-NULL
  if (keep_ECD==FALSE){
  df$ECD<-NULL}
  return(df)
}

group_abundance_func<- function(df, grouping_vars=c('season', 'depth_mid', 'lat_mid', 'lon_mid'), env_vars=c('Temp', 'Salinity', 'Fluoro', 'Oxygen', 'Irradiance'), group_by_taxa =TRUE, ISIIS_ver=2, thresholding = TRUE, count_thres_prob=3){
	#this function will use any approach to define a plankton community, and assign a unique identifier to each community (Group_ID) and then will estimate the abundance (biovol/L^1) of each size class (volume is in cubic meters)
   if (ISIIS_ver==2){
     pixel_size<-66
     vol_frame<-0.00845
   } else if (ISIIS_ver==3) {
     pixel_size<-55
     vol_frame<- 0.003 # 10 cm * 10 cm * 30 cm  
   }

    # Step 1: Frame separation and initial volume estimate
    df<-df %>% separate(img, into=c('frame_id', 'img'), sep='_crop_', remove=TRUE) %>% group_by(!!!syms(grouping_vars)) %>% mutate(n_frames=n_distinct(frame_id)) %>%as.data.frame() # volume sampled needs to be based on environmental variables only, not on taxa, else the volume sampled would be too low for rare taxa
    df$vol_sampled<- df$n_frames*vol_frame

    grouping_vars <- c(grouping_vars, "group")
  
    # Step 3: Assign group ID and average environmental vars across group
    #df<-df %>% group_by(!!!syms(grouping_vars)) %>% mutate(group_id = cur_group_id())%>%as.data.frame() #assign group_id, now a separate function, since binning already separates files into groups
    df <- df %>% 
          group_by(!!!syms(grouping_vars)) %>% 
          mutate(group_id = cur_group_id()) %>% 
          mutate(
            across(
              all_of(env_vars),
              ~ if (cur_column() %in% c("Latitude", "Longitude")) {
                  median(.x, na.rm = TRUE)
                } else if (cur_column() == "time_Pacific") {
                  sort(.x)[ceiling(dplyr::n() / 2)]
                } else {
                  mean(.x, na.rm = TRUE)
                }
            )
          ) %>% 
          ungroup() %>% 
          as.data.frame() 
    df <- df[order(df$group_id),]
    
        # Step 4.1 Summarize into size classes, keeping env vars (assumed constant per group)
    #df<-df %>% group_by(!!!syms(c('group_id', grouping_vars,'vol_sampled', 'biovol_um3_width', 'biovol_um3_mid')))  %>% summarise(n_individuals=n(), total_biovolume=sum(biovol)) %>%as.data.frame()# this is where the date and time info is lost, consider if there is an advantage of retaining this information somehow
    df <- df %>%
      group_by(!!!syms(c('group_id', grouping_vars, 'biovol_um3_width', 'biovol_um3_mid'))) %>%
      summarise(
        n_individuals = n(),
        total_biovolume = sum(biovol),
        vol_sampled = first(vol_sampled),
        across(all_of(env_vars), ~first(.x))
      ) %>%
      as.data.frame()
    df<- df %>% apply_correction_factors(., base_dir,  project=project, binning = 'size_spectra')
    #df <- df %>% filter(total_biovolume != 0) # some sizes are 0 (especially with contour)
    
    #after applying correction factors, remove grouping_vars:
    grouping_vars <- grouping_vars[grouping_vars != "group"]
 
    if (isTRUE(group_by_taxa)) {
        df$group <- as.character(df$group)
        # keep NA as NA (don’t force to "unknown" unless you want that)
        df$grouped_taxon <- ifelse(is.na(df$group), NA_character_,
                                 stringr::word(df$group, 1, sep = "_"))
        grouping_vars <- c(grouping_vars, "grouped_taxon")
    }

    # redo additions after applying correction factors
    df <- df %>%
      group_by(!!!syms(grouping_vars)) %>%
      mutate(group_id = cur_group_id()) %>%
      mutate(across(all_of(env_vars), ~mean(.x, na.rm = TRUE)))%>%
      ungroup() %>%
      as.data.frame()
    
    df <- df[order(df$group_id),]

    df <- df %>%
      group_by(!!!syms(c('group_id', grouping_vars, 'biovol_um3_width', 'biovol_um3_mid'))) %>%
      summarise(
        n_individuals = sum(n_individuals),
        total_biovolume = sum(total_biovolume),
        vol_sampled = first(vol_sampled),
        across(all_of(env_vars), ~first(.x))
      ) %>%
      as.data.frame()


    # Step 4.2 Reload the ecotaxa size bins to assign empty size bins
    size_classes <- read.csv(file.path(base_dir, "ancillary", "contour_custom_ISIIS-3_size_classes.csv"), sep = ',')
    #all_mids <- as.numeric(midpoints(cut(size_classes$biovol_um3, breaks = size_classes$biovol_um3)))
    #all_mids <- as.numeric(as.character(midpoints(cut(size_classes$biovol_um3, breaks = size_classes$biovol_um3))))
    all_mids <- as.numeric(size_classes$biovol_um3_mid)

    
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
            ind_m3 = ifelse(is.na(n_individuals) | n_individuals == 0 | is.na(biovol_um3_width), NaN, n_individuals / vol_sampled),
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
      #stop("💥 Non-finite NB encountered — check inputs.")
    }
    
    # old code included in the mutate command in line 281
    #df$NB<- (df$total_biovolume/df$vol_sampled)/df$biovol_um3_width
    #df$ECD_um_mid<-((df$biovol_um3_mid*6)/pi)**(1./3.)
    #df$pixel_size_mid<-(df$ECD_um_mid/1000) * (1/(pixel_size/1000))# following Mathilde's method, we get how many pixels per mm there are in each ESD size class 
    #df$size_uncertainty<- dnorm((1/3)*pi*df$pixel_size_mid^3, mean=(1/3)*pi*60^3,sd=(1/3)*60*pi)#mean=(1/3)*pi*min(df$pixel_size_mid, na.rm=T)^3 now I set up the mean to 60 pixels, the smallest segments on training lib
    #df$ECD_um_width<-((df$biovol_um3_width*6)/pi)**(1./3.)
    #df$ind_L<- (df$n_individuals/df$vol_sampled)#/df$ECD_um_width add this if you want instead PSD
    #df$logNB<-log10(df$NB)
    #df$logBiovol<-log10(df$biovol_um3_mid)
    #df$count_uncertainty<-dpois(x=df$n_individuals, lambda=df$n_individuals)
    if (thresholding ==TRUE){
        df<-threshold_func(df, threshold_count=count_thres_prob)
    }
    return(df)
}
## integrate these two lines below to the function above:
#NBS_biovol_df=NBS_biovol_df.assign(count_uncertainty=poisson.pmf(k=NBS_biovol_df['ROI_number_sum'], mu=NBS_biovol_df['ROI_number_sum']),
#                                     size_uncertainty=NBS_biovol_df.astype({'size_class_pixel':float}).groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth','Max_obs_depth']).size_class_pixel.transform(lambda x: norm.pdf((1/3)*np.pi*x**3,loc=(1/3)*np.pi*x.min()**3,scale=(1/3)*1*np.pi)).values.tolist())
# for poisson: https://makemeanalyst.com/statistics-with-r/poisson-functions-in-r-programming/#:~:text=The%20probability%20mass%20function%20(PMF,is%20the%20factorial%20of%20k.
# for gaussian: https://r02pro.github.io/normal-distribution.html, Mathilde did it in pixels per mm, but the size_class_pixel is each size class in pixels

threshold_func<-function(df, threshold_size=0.2, threshold_count=0.3){
  n_before <- nrow(df)

  df <- df[which(df$size_uncertainty <= threshold_size), ]
  df <- df[which(df$count_uncertainty <= threshold_count), ]
  df <- df[which(df$pixel_size_mid >= 1), ]
  # per group_id, find the biovol_mid at the max NB (if ties, take the largest mid)
  cuts <- df %>%
    dplyr::group_by(group_id) %>%
    dplyr::slice_max(NB, n = 1, with_ties = TRUE, na_rm = TRUE) %>%
    dplyr::summarise(cut_mid = max(biovol_um3_mid, na.rm = TRUE), .groups = "drop")
  
  # keep only rows with mid >= that threshold
  df <- df %>%
    dplyr::left_join(cuts, by = "group_id") %>%
    dplyr::filter(biovol_um3_mid >= cut_mid) %>%
    dplyr::select(-cut_mid)
    n_after <- nrow(df)

  #drop rows after three consecutive size classes with NB=NA
  df <- df %>%
    dplyr::mutate(biovol_um3_mid = as.numeric(biovol_um3_mid)) %>%
    dplyr::arrange(group_id, biovol_um3_mid) %>%
    dplyr::group_by(group_id) %>%
    dplyr::group_modify(function(.x, .g) {
      na_vec <- is.na(.x$NB)
      r      <- rle(na_vec)
      pos    <- which(r$values & r$lengths >= 3)[1]  # first 3-NA run
      keep_n <- if (is.na(pos)) nrow(.x) else sum(r$lengths[seq_len(pos - 1)])
      if (keep_n == 0) .x[0, , drop = FALSE] else .x[seq_len(keep_n), , drop = FALSE]
    }) %>%
    dplyr::ungroup()


  message("📉 Thresholding removed ", n_before - n_after, " rows (", n_after, " remaining)")

  return(df)
# also include this on the function, which is cutting the size spectra based on empty size bins, check https://masterr.org/r/how-to-find-consecutive-repeats-in-r/


}

group_id_func <- function(df, grouping_vars = c('season', 'depth_mid', 'lat_mid', 'lon_mid'), group_by_taxa = TRUE) {
  if (group_by_taxa == TRUE) {
    grouping_vars <- c(grouping_vars, 'predicted_taxon')
  }
  
  df <- df %>% group_by(!!!syms(grouping_vars)) %>%
               mutate(group_id = cur_group_id()) %>%
               ungroup() %>%
               arrange(group_id) %>%
               as.data.frame()
  
  return(df)
}
    

lin_fit_func<- function(df, filter_min_size_class=FALSE, n_min_size_class=3){  
    
    fitted_models <- df %>% filter(is.finite(logNB), is.finite(logBiovol)) %>% group_by(group_id) %>% do(model = lm(logNB ~ logBiovol, data = .)) 

    # using augmented linear fit to extract residuals
    augmented <- df %>% filter(is.finite(logNB), is.finite(logBiovol)) %>% group_by(group_id) %>% group_modify(~ augment(lm(logNB ~ logBiovol, data = .x), data = .x)) %>% ungroup()
    bad <- df %>% filter(!is.finite(logNB) | !is.finite(logBiovol))
    if (nrow(bad) > 0) {
        message("⚠️ Dropping ", nrow(bad), " rows with non-finite logNB or logBiovol before fitting")
  
        bad_summary <- bad %>% select(group_id, logNB, logBiovol, NB, total_biovolume, vol_sampled, biovol_um3_width, biovol_um3_mid) %>% distinct() %>% arrange(group_id)
  
        message("🔎 Affected group summary:\n", 
          capture.output(print(bad_summary)) %>% paste(collapse = "\n"))
    }   
    valid_groups <- fitted_models$group_id
    model_list <- fitted_models$model
    
    ## this is the structure of how the data is obstained: fitted_models[[2]][['group_id']]$coefficients[1] for intercept, coefficients[2] for slope
    for (i in seq_along(valid_groups)) {
      g <- valid_groups[i]
      model_entry <- model_list[[i]]
    
      df$intercept[df$group_id == g] <- model_entry$coefficients[1]
      df$slope[df$group_id == g]     <- model_entry$coefficients[2]
      df$R2[df$group_id == g]        <- summary(model_entry)$r.squared
      df$n_size_classes[df$group_id == g] <- nrow(model_entry$model)
    }   
    
    df <- df %>% left_join(
                 augmented %>% select(group_id, biovol_um3_mid, biovol_um3_width, .fitted, .resid, .std.resid),
                   by = c("group_id", "biovol_um3_mid", "biovol_um3_width")
           ) %>% rename( fitted_logNB = .fitted, residual_logNB = .resid, std_residual_logNB = .std.resid)

        
    if (filter_min_size_class==TRUE){
       df <- df[order(df$group_id),]  %>% filter(n_size_classes> n_min_size_class)}
    return(df)
}
