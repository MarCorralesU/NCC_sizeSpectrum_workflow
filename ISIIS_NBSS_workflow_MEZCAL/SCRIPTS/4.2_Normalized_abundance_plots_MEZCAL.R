setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_MEZCAL')
library(ggplot2)
library(colorRamps)
library(plyr)
library(dplyr)
library(reshape2)
library(akima)
library(gstat)
library(raster)
library(randomForest)   # not used unless you enable the RF section
library(scales)
library(metR)
library(kriging)
library(ggthemes)
library(optparse)
library(purrr)
library(readr)
library(fuzzyjoin)
library(RColorBrewer)
library(colorspace)
base_dir<- getwd()

source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


option_list <- list(
  make_option(c("-c", "--Cruise"), type = "character", help = "Cruise: w22/s22/w23/s23", metavar = "character"),
  make_option(c("-t", "--Transect"), type = "character", help = "Cruise: GH/OB/CR/CM/NH/HH/RR", metavar = "character"),
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)

if (is.null(opt$Cruise)) {
  stop("Error: --Cruise is required.")
}

#variables to load data:
Cruise <- opt$Cruise # Cruise<- 's22'
Transect <- opt$Transect # Transect<- 'NH'

size_var<- opt$sizeVar # size_var<-'contour'
instrument <- 'ISIIS'
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "ancillary", "transect_shelf_land.csv")
sizeSpectraDir <- file.path(base_dir,  'data', 'size_spectra_allParticles_5-10m-bins_thresholded', size_var, Cruise, Transect)



output_dir_figs_NB <- file.path(base_dir, "figures", 'size_spectra', instrument, Cruise, Transect)
if (!dir.exists(output_dir_figs_NB)) {
  dir.create(output_dir_figs_NB, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_NB)
}

output_dir_figs_NB_taxa <- file.path(output_dir_figs_NB, instrument, 'NB_taxa')
if (!dir.exists(output_dir_figs_NB_taxa)) {
  dir.create(output_dir_figs_NB_taxa, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_NB_taxa)
}



if (instrument=='ISIIS'){
  df_list <- list.files(sizeSpectraDir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
  dat_ISIIS<- map_dfr( df_list, ~ read_csv( .x, show_col_types = FALSE))
  dat_ISIIS <- dat_ISIIS %>%
    mutate(time_Pacific_orig = time_Pacific,
           mean_time_Pacific = time_Pacific - hours(1))
  
  
  dat_ISIIS <- dat_ISIIS %>%
    # keep the original for reference
    dplyr::rename(group_id_old = group_id) %>%
    # new globally consecutive id: increments whenever group_id_old changes
    dplyr::mutate(
      group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE))
    )
  
  
  #get integrated biovolume:
  x1 <- 12.78849 # minimum range contemplated by  LISST
  x2 <- 4.547442e13 #maximum range contemplated by ISIIS
  
  dat_ISIIS <- dat_ISIIS %>%
    dplyr::mutate(
      slope = as.numeric(slope),
      intercept = as.numeric(intercept),
      b_lin = 10^intercept,
      integrated_biovol = dplyr::if_else(
        is.finite(slope) & is.finite(b_lin),
        ifelse(
          abs(slope + 1) < 1e-12,
          b_lin * log(x2 / x1),
          (b_lin / (slope + 1)) * (x2^(slope + 1) - x1^(slope + 1))
        ),
        NA_real_
      )
    )
  
  # small plot to illustrate differences in slope metrics (ONLY for s22 NH)
  if (Cruise=='s22' && Transect=='NH'){
    isiis_group <- dat_ISIIS %>%
      dplyr::group_by(group_id) %>%
      dplyr::summarise(
        slope_g    = median(as.numeric(slope), na.rm = TRUE),
        rel_dist_g = median(as.numeric(rel_dist), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(is.finite(slope_g), is.finite(rel_dist_g))
    
    # 1) Group closest to slope = -0.3 with rel_dist > 0
    pick_shallow <- isiis_group %>%
      dplyr::filter(rel_dist_g > 0) %>%
      dplyr::mutate(score = abs(slope_g - (-0.5))) %>%
      dplyr::arrange(score) %>%
      dplyr::slice(1)
    
    # 2) Group closest to slope = -1.8 with rel_dist closest to -1.1
    pick_steep <- isiis_group %>%
      dplyr::mutate(
        score = abs(slope_g - (-1.5)) + abs(rel_dist_g - (-1.1))
      ) %>%
      dplyr::arrange(score) %>%
      dplyr::slice(1)
    
    pick_shallow
    pick_steep
    
    id_shallow <- pick_shallow$group_id
    id_steep   <- pick_steep$group_id
    
    
    dat_ISIIS_subset <- dat_ISIIS %>%
      dplyr::filter(group_id %in% c(id_shallow, id_steep))

    # ) Minimal data (log space) + per-group params
    df_obs <- dat_ISIIS_subset %>%
      dplyr::transmute(
        group_id   = group_id,
        logBiovol  = as.numeric(logBiovol),
        logNB      = as.numeric(logNB),
        slope      = as.numeric(slope),
        intercept  = as.numeric(intercept)   # intercept in log10 space
      ) %>%
      dplyr::filter(
        is.finite(logBiovol),
        is.finite(logNB),
        is.finite(slope),
        is.finite(intercept)
      )
    params <- df_obs %>%
      dplyr::group_by(group_id) %>%
      dplyr::summarise(
        slope     = dplyr::first(slope),
        intercept = dplyr::first(intercept),
        .groups = "drop"
      )
    
    # x-ranges in log space for each group
    x_ranges <- df_obs %>%
      dplyr::group_by(group_id) %>%
      dplyr::summarise(
        x_min = min(logBiovol, na.rm = TRUE),
        x_max = max(logBiovol, na.rm = TRUE),
        .groups = "drop"
      )
    
    # fitted lines (log-log linear): logNB_fit = intercept + slope * logBiovol
    df_fit <- x_ranges %>%
      dplyr::left_join(params, by = "group_id") %>%
      dplyr::rowwise() %>%
      do({
        g <- .
        x <- seq(g$x_min, g$x_max, length.out = 200)
        tibble(
          group_id   = g$group_id,
          logBiovol  = x,
          logNB_fit  = g$intercept + g$slope * x,
          slope      = g$slope         # <- map color to continuous slope
        )
      }) %>%
      dplyr::ungroup()
    
    # Optional: equation labels (kept; remove this whole block if you don't want labels)
    fmt_sig <- function(x, sig = 2) format(signif(x, sig), scientific = TRUE)
    eq_labels <- params %>%
      dplyr::transmute(
        group_id,
        label = paste0(
          "NB == ", format(signif(10^intercept, 2), scientific = TRUE),
          " * x^{", format(signif(slope, 3), scientific = TRUE), "}"
        )
      )
    
    label_pos <- df_fit %>%
      dplyr::group_by(group_id) %>%
      dplyr::summarise(
        logBiovol = max(logBiovol, na.rm = TRUE),
        logNB_fit = max(logNB_fit, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(params %>% dplyr::select(group_id, slope), by = "group_id") %>%
      dplyr::left_join(eq_labels, by = "group_id") %>%
      dplyr::mutate(
        logNB_fit = dplyr::if_else(slope > -1, logNB_fit + 0.15, logNB_fit - 0.15)
      )
    
    # Colors for RdBu with white center at -1
    low_col  <- "#D7191C"  # vivid red
    high_col <- "#2C7BB6"  # vivid blue

    # Span the plotted x-range (in log10 space)
    xr <- range(df_fit$logBiovol, na.rm = TRUE)
    
    ideal_df <- tibble(
      logBiovol = xr,
      logNB_fit = 10 + (-1) * xr   # intercept = log10(1e10) = 10
    )
    
    ideal_label_df <- tibble(
      logBiovol = max(xr, na.rm = TRUE),
      logNB_fit = 11.4 + (-1) * max(xr, na.rm = TRUE),
      label = 'paste("Equal biomass distribution (", NB == x^{-1}, ")")'
    )
    
    
    p <- ggplot() +
      # 1) ONLY the fitted functions (no raw data)
      geom_line(
        data = df_fit,
        aes(x = logBiovol, y = logNB_fit, group = group_id, color = slope),
        linewidth = 1.5, linetype = "dashed", alpha = 0.95
      ) +
      # 2) Idealized line: slope = -1, intercept = -10 (in log10 space)
      geom_line(
        data = ideal_df,
        aes(x = logBiovol, y = logNB_fit),
        linetype = "dotted", linewidth = 1.5, color = "black"
      ) + 
      geom_text(
        data = ideal_label_df,
        aes(x = logBiovol, y = logNB_fit, label = label),
        color = "black",
        hjust = 1.1, vjust = -0.4,
        size = 5,
        parse = TRUE,
        inherit.aes = FALSE
      ) +
      # (Optional) equation labels; remove if not desired
      geom_text(
        data = label_pos,
        aes(x = logBiovol, y = logNB_fit, label = label, color = slope),
        parse = TRUE, hjust = 1, size = 5, fontface = "bold", show.legend = FALSE
      ) +
      # 3) RdBu diverging palette centered at -1
      scale_color_gradient2(
        low = high_col,
        mid = "white",
        high = low_col,
        midpoint = -1,
        name = "slope"
      ) +
      labs(
        x = expression(log[10]*"(size, " * mu * "m"^3 * ")"),
        y = expression(log[10]*"(Abundance, " * "m"^-3 * ")")
      ) +
      theme_classic() +
      theme(
        axis.text  = element_text(size = 0),
        axis.title = element_text(size = 15)
      )
    
    p
    
    # plot ISIIS and LISST together (EXCLUSIVE OF ISIIS PLOTTING)
    out_path <- file.path(output_dir_figs_NB, paste(paste("example_specific_samples_equation", Cruise, Transect, sep="_"), ".png"))
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    
    # if your plot object is named `p`:
    ggsave(filename = out_path, plot = p, width = 6, height = 6, units = "in", dpi = 300, bg = "white")
    
    
  }
  

  
  isiis <- dat_ISIIS %>%
    dplyr::transmute(
      source   = "ISIIS",
      group_id = group_id,
      size_num = as.numeric(ECD_um_mid),
      NB       = as.numeric(NB)
    ) %>%
    dplyr::filter(is.finite(size_num), size_num > 0,
                  is.finite(NB), NB > 0)
  
  # 1) Shared ordered factor for size bins
  levels_x <- sort(unique(c(isiis$size_num)))
  isiis <- isiis %>% dplyr::mutate(size_fac = factor(size_num, levels = levels_x))
  
  # 2) Individual observations (points only; no per-group lines)
  pts_df <- dplyr::bind_rows(isiis) %>%
    dplyr::arrange(source, group_id, size_fac)
  
  # 3) Median (or mean) NB per size class, per source
  # I recommend MEDIAN for robustness; switch to mean(.) if you truly want mean.
  med_df <- pts_df %>%
    dplyr::group_by(source, size_fac) %>%
    dplyr::summarise(
      NB_med = median(NB, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(size_num = as.numeric(as.character(size_fac))) %>%
    dplyr::arrange(source, size_fac)
  
  # 4) Equation annotations from median intercept & slope
  fmt_num <- function(x) formatC(x, format = "fg", digits = 3)
  isiis_b <- fmt_num(stats::median(dat_ISIIS$intercept, na.rm = TRUE))
  isiis_s <- fmt_num(stats::median(dat_ISIIS$slope,     na.rm = TRUE))
  
  lbl_isiis <- paste0("y==", isiis_b, " * x^{", isiis_s, "}")
  
  # 5) Y positions (robust quantiles so labels sit near the clouds)
  isiis_y <- stats::quantile(pts_df$NB[pts_df$source == "ISIIS"], 0.20, na.rm = TRUE)
  
  # 6) Plot
  src_cols <- c(ISIIS = "royalblue3")
  
  p <- ggplot() +
    # individual samples: open circles, no connecting lines
    geom_point(
      data = pts_df,
      aes(x = size_num, y = NB, color = source),
      shape = 1, size = 1.3, stroke = 0.45, alpha = 0.3
    ) +
    # solid median lines (two lines, colored by source)
    geom_line(
      data = med_df,
      aes(x = size_num, y = NB_med, color = source),
      linewidth = 1.3, linetype = "solid", alpha = 0.95
    ) +
    # optional: median points as well (often reads nicer)
    geom_point(
      data = med_df,
      aes(x = size_num, y = NB_med, color = source),
      size = 1, alpha = 0.95
    ) +
    # equation annotations
    ggplot2::annotate(
      geom = "text", x = 70, y = isiis_y, label = lbl_isiis,
      color = src_cols["ISIIS"], parse = TRUE, hjust = 0,
      size = 8, fontface = "bold"
    ) +
    scale_color_manual(values = src_cols, name = "Source") +
    scale_x_log10() +
    scale_y_log10(limits = c(1e-5, 1e11)) +
    labs(
      x = expression(Equivalent~spherical~diameter~(mu*m)),
      y = expression(Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))
    ) +
    theme_classic() +
    theme(
      axis.text   = element_text(size = 16),
      axis.title  = element_text(size = 18),
      legend.title= element_text(size = 16),
      legend.text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.6, alpha = 1)))
  
  p
  
  
  out_path <- file.path(output_dir_figs_NB, paste(paste("ISIIS_comparison", Cruise, Transect, sep="_"), ".png"))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  # if your plot object is named `p`:
  ggsave(filename = out_path, plot = p, width = 9, height = 7, units = "in", dpi = 300, bg = "white")
  
  
  tol_sec   <- 60
  tol_depth <- 0.5
  

  
  # -- Mormalized abundance per taxa
  
  #load data
  df_list_taxa <- list.files(sizeSpectraDir, pattern = "by_taxa\\.csv$", full.names = TRUE, recursive = TRUE)
  dat_full_taxa <- map_dfr( df_list_taxa, ~ read_csv( .x, show_col_types = FALSE))
  
  # redefine group_id 
  setDT(dat_full_taxa)[, group_id := .GRP, .(cast_id, depth_mid)]
  # average Lat and Long by group_id (so that we can sum NBs across taxa):
  dat_full_taxa <- dat_full_taxa %>%
    dplyr::group_by(group_id) %>%
    dplyr::mutate(
      Latitude  = mean(as.numeric(Latitude),  na.rm = TRUE),
      Longitude = mean(as.numeric(Longitude), na.rm = TRUE),
    ) %>%
    dplyr::ungroup()
  
  # stacked taxa plots:
  
  for (gid in unique(dat_full_taxa$group_id)){
    
    taxon_cols <- c(
      "appendicularian" = "#F0E442",
      "chaetognath"     = "#E69F00",
      "copepod"         = "#CC79A7",
      "crustacean"      = "#D55E00",
      "ctenophore"      = "#8C564B",
      "diatoms"         = "#0072B2",
      "phytoplankton"   = "#009E73",
      "protist"         = "#9467BD",
      "pteropod"        = "#56B4E9",
      "fish"            = "#999999",
      "polychaete"      = "#B2DF8A",
      "hydromedusae"    = "#FF9896",
      "siphonophore"    = "#17BECF",
      "echinoderm"      = "#C5B0D5",
      "detritus"        = "#8B4513",  # brown
      "tunicate"        = "#66C2A5",  # seafoam green
      "pyrosome"        = "#E7298A",   # magenta
      "foraminifera"    = "#D2B48C"   # tan
    )
    # 1) take the rows for THIS group_id
    df_sub <- dat_full_taxa %>% dplyr::filter(group_id == gid)
    
    # 2) grab the “metadata” for naming
    first_combo <- df_sub %>% dplyr::slice(1)
    
    plot_id <- paste(first_combo$Transect_ID,
                     Cruise,                              
                     paste(first_combo$depth_mid, "m", sep = "-"),
                     round(first_combo$Longitude, 2),
                     round(first_combo$rel_dist, 2),
                     sep = "_")
    
    message("plotting ", plot_id)
    
    # 2) Sum *NB* (linear space) per size bin & taxon
    plot_df <- df_sub %>%
      group_by(ECD_um_mid, grouped_taxon) %>%
      summarise(NB_sum = sum(NB, na.rm = TRUE), .groups = "drop")
    
    plot_df <- plot_df %>%
      dplyr::group_by(ECD_um_mid) %>%
      dplyr::mutate(
        total_NB  = sum(as.numeric(NB_sum),  na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
    
    plot_df$per_totalNB<- (plot_df$NB_sum/plot_df$total_NB)*100
    
    x_levels <- sort(unique(plot_df$ECD_um_mid))
    
    
    # totals: one per bin
    totals <- plot_df |> distinct(ECD_um_mid, total_NB)
    
    # taxon series (one per bin × taxon), keep positives for log-y
    taxa_series <- plot_df |>
      dplyr::group_by(ECD_um_mid, grouped_taxon) |>
      dplyr::summarise(NB_taxon = sum(NB_sum, na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(NB_taxon > 0)
    
    # ---- compute per-bin bar slots on a log-x axis ----
    # geometric spacing ratio between adjacent mids
    mids <- sort(unique(taxa_series$ECD_um_mid))
    r    <- median(mids[-1] / mids[-length(mids)], na.rm = TRUE)      # ~2 if doubling bins
    
    taxa_series <- taxa_series |>
      group_by(ECD_um_mid) |>
      arrange(grouped_taxon, .by_group = TRUE) |>
      mutate(
        n_bin   = n(),                                   # #taxa in this bin
        span    = ECD_um_mid * (sqrt(r) - 1/sqrt(r)) * 0.6,  # total horizontal span per bin
        w       = span / n_bin,                          # bar width inside the bin
        slot    = row_number() - 0.5,                    # 0.5, 1.5, 2.5, ...
        xmin    = ECD_um_mid - span/2 + (slot - 1) * w,  # left edge
        xmax    = xmin + w,                              # right edge
        ymin    = min(NB_taxon[NB_taxon > 0]) * 0.5      # small positive base for log-y
      ) |>
      ungroup()
    
    # y-axis breaks from totals
    yr <- range(totals$total_NB[totals$total_NB > 0], na.rm = TRUE)
    ybrks <- 10^seq(floor(log10(yr[1])), ceiling(log10(yr[2])))
    
    # x-axis breaks you like
    xbrks <- c(100, 500, 1000, 2500, 5000, 10000, 20000)
    
    # pick a floor one decade below the smallest positive segment
    pos_min <- min(taxa_series$NB_taxon[taxa_series$NB_taxon > 0], na.rm = TRUE)
    y_floor <- 10^(floor(log10(pos_min)) - 1)
    
    # make all bars start at the common floor
    taxa_series <- taxa_series %>% mutate(ymin = y_floor)
    
    # recompute breaks/limits using that floor
    ymax  <- max(totals$total_NB, na.rm = TRUE)
    ybrks <- 10^seq(floor(log10(y_floor)), ceiling(log10(ymax)))
    
    x_limits <- c(100, 30000)
    
    missing_taxa <- setdiff(unique(taxa_series$grouped_taxon), names(taxon_cols))
    
    if (length(missing_taxa) > 0) {
      message("Unmapped taxa in ", plot_id, ": ",
              paste(missing_taxa, collapse = ", "))
    }
    
    p <- ggplot() +
      geom_rect(
        data = taxa_series,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = NB_taxon, fill = grouped_taxon),
        alpha = 0.85
      ) +
      geom_line(data = totals, aes(ECD_um_mid, total_NB), color = "black", linewidth = 0.6) +
      geom_point(data = totals, aes(ECD_um_mid, total_NB), color = "black", shape = 16, size = 2.6) +
      scale_x_log10(
        limits = x_limits,                               # <- add this
        breaks = xbrks,
        labels = scales::label_number(big.mark = ","),
        expand = c(0.02, 0)
      ) +
      scale_y_log10(
        limits = c(y_floor, ymax),
        breaks = ybrks,
        labels = scales::label_math(10^.x),
        expand = c(0, 0)
      ) +
      labs(
        x = expression(Equivalent~spherical~diameter~(mu*m)),
        y = expression(Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))
      ) +
      scale_fill_manual(
        values = c(
          "appendicularian" = "#F0E442",
          "chaetognath"     = "#E69F00",
          "copepod"         = "#CC79A7",
          "crustacean"      = "#D55E00",
          "ctenophore"      = "#8C564B",
          "diatoms"         = "#0072B2",
          "phytoplankton"   = "#009E73",
          "protist"         = "#9467BD",
          "pteropod"        = "#56B4E9",
          "fish"            = "#999999",
          "polychaete"      = "#B2DF8A",
          "hydromedusae"    = "#FF9896",
          "siphonophore"    = "#17BECF",
          "echinoderm"      = "#C5B0D5",
          "detritus"        = "#8B4513",  # brown
          "tunicate"        = "#66C2A5",  # seafoam green
          "pyrosome"        = "#E7298A",  # magenta
          "foraminifera"    = "#D2B48C"   # tan
        ),
        name = "Predicted taxon"
      ) +
      theme_classic(base_size = 12)
    
    ggsave( filename = file.path(output_dir_figs_NB_taxa, paste(plot_id, "NB_taxa.png", sep="_")), plot = p,
            width  = 12, height = 6, units = "in", dpi = 300)
    
    
    ## PERCENTAGE BARS!!
    
    # # per-bin × taxon totals and shares
    # 
    # # totals for the black line/dots
    # 
    # # per-bin width on a log-x axis (one bar per bin)
    # 
    # # map Total (log10) into [0,1] so it overlays % bars; right axis shows real values
    # 
    #   # 1) one stacked bar per ECD bin showing relative contribution
    #   # 2) overlay Total NB (black) scaled into [0,1] space
    #   # log-x ticks
    
    
  }
  
  

  
  
}  





## DATA SOURCE ONLY:
# 1) Minimal ISIIS frame (positives only for log axes)
#filter out LISST data with very high slopes
dat_ISIIS<- dat_ISIIS[dat_ISIIS$slope>1,]
isiis <- dat_ISIIS %>%
  dplyr::transmute(
    group_id  = group_id,
    size_num  = as.numeric(ECD_um_mid),
    NB        = as.numeric(NB),
    intercept = as.numeric(intercept),
    slope     = as.numeric(slope)
  ) %>%
  dplyr::filter(is.finite(size_num), size_num > 0,
                is.finite(NB), NB > 0)

# 2) Ordered factor for exact size bins
levels_x <- sort(unique(isiis$size_num))
isiis <- isiis %>% dplyr::mutate(size_fac = factor(size_num, levels = levels_x))

# 3) Collapse to ONE row per (group, size) and compute per-group slope (for color)
isiis_lines <- isiis %>%
  dplyr::group_by(group_id, size_fac) %>%
  dplyr::summarise(NB = median(NB, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(
    isiis %>% dplyr::group_by(group_id) %>%
      dplyr::summarise(group_slope = median(slope, na.rm = TRUE), .groups = "drop"),
    by = "group_id"
  ) %>%
  dplyr::mutate(size_num = as.numeric(as.character(size_fac))) %>%
  dplyr::arrange(group_id, size_fac)

# 4) Median spectrum across groups (per size bin)
isiis_median <- isiis_lines %>%
  dplyr::group_by(size_fac) %>%
  dplyr::summarise(NB_med = median(NB, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(size_num = as.numeric(as.character(size_fac))) %>%
  dplyr::arrange(size_fac)

# 5) Equation label from median intercept & slope (use y and x)
fmt_num <- function(x) formatC(x, format = "fg", digits = 3)
b_med <- fmt_num(stats::median(dat_ISIIS$intercept, na.rm = TRUE))
s_med <- fmt_num(stats::median(dat_ISIIS$slope,     na.rm = TRUE))
lbl_eq <- paste0("y==", b_med, " * x^{", s_med, "}")

## Dynamic x for the annotation (80% of range) and y near top of range
x_ann <- 10^(mean(range(log10(isiis_lines$size_num), na.rm = TRUE)) + 0.3 *
               diff(range(log10(isiis_lines$size_num), na.rm = TRUE)))
y_ann <- 10^(-0.2)  # ~0.63 * 1e-1 on log scale; adjust if needed

p<- ggplot() +
  # slightly more visible per-group lines (still very transparent)
  geom_line(
    data = isiis_lines,
    aes(x = size_num, y = NB, group = group_id, color = group_slope),
    linewidth = 0.35, alpha = 0.12
  ) +
  # black dashed median spectrum
  geom_line(
    data = isiis_median,
    aes(x = size_num, y = NB_med),
    color = "black", linewidth = 1.2, linetype = "dashed"
  ) +
  # black equation text (using y and x)
  ggplot2::annotate("text", x = 5000, y = 1e-2, label = lbl_eq,
           color = "black", parse = TRUE, hjust = 0,
           size = 8, fontface = "bold") +
  scale_color_viridis_c(name = "Slope",
                        guide = guide_colorbar(barheight = grid::unit(70, "pt"),
                                               barwidth  = grid::unit(12, "pt"))) +
  scale_x_log10(limits = c(min(levels_x, na.rm = TRUE), max(levels_x, na.rm = TRUE))) +
  scale_y_log10() +
  labs(
    x = expression(Equivalent~spherical~diameter~(mu*m)),
    y = expression(Normalized~Biovolume~(mu*m^{3}~m^{-3}~mu*m^{-3}))
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 16),
    axis.title  = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.title= element_text(size = 16)
  )

p



out_path <- file.path(output_dir_figs_NB, paste(paste("NBSS", Cruise, Transect, instrument, sep="_"), ".png", sep=""))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

# if your plot object is named `p`:
ggsave(filename = out_path, plot = p, width = 9, height = 7, units = "in", dpi = 300, bg = "white")

