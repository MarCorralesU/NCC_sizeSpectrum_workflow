# summary size spectra 
setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
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
base_dir<- getwd()

source(file.path(base_dir, 'SCRIPTS', 'functions', 'ISIIS_PSS_binning_funcs.R'))


option_list <- list(
  make_option(c("-s", "--sizeVar"), type = "character", default = "contour", help = "the size variable to be used to bin by size", metavar = "character"),
  make_option(c("-e", "--envVars"), type = "character", default = "Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l", help = "Comma-separated environmental variables", metavar = "character"),
  make_option(c("-d", "--distFromShore"), type = "character", default = "ancillary/transect_shelf_land.csv", help = "Path to CSV with longitude shore locations  per transect", metavar = "character"))

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print and validate
print(opt)


size_var<- opt$sizeVar # size_var<-'contour'
env_vars <- strsplit(opt$envVars, ",")[[1]] # env_vars<- strsplit('Temperature,Salinity,Oxygen_ml_l,pH,PAR,chl_a_ug_l', ",")[[1]]
shore_path <- file.path(base_dir, opt$distFromShore) # shore_path  <- file.path(base_dir, "ancillary", "transect_shelf_land.csv")
sizeSpectraDir <- file.path(base_dir,  'data', 'size_spectra') 
LISSTdir <-  file.path(base_dir,  'data', "LISST_data", "LISST_size_spectra", "5m_depthBin") 


output_dir_figs_NB <- file.path(base_dir, "figures", 'size_spectra')
if (!dir.exists(output_dir_figs_NB)) {
  dir.create(output_dir_figs_NB, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", output_dir_figs_NB)
}


df_list_ISIIS <- list.files(sizeSpectraDir, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)
df_list_LISST <- list.files(LISSTdir, pattern = "5m-depthBin_LISST_size_spectra\\.csv$", full.names = TRUE, recursive = TRUE)

dat_ISIIS<- map_dfr( df_list_ISIIS, ~ read_csv( .x, show_col_types = FALSE))
dat_ISIIS <- dat_ISIIS %>%
  mutate(time_Pacific_orig = time_Pacific,
         mean_time_Pacific = time_Pacific - hours(1))

dat_LISST<- map_dfr( df_list_LISST, ~ read_csv( .x, show_col_types = FALSE))


dat_ISIIS <- dat_ISIIS %>%
  # keep the original for reference
  dplyr::rename(group_id_old = group_id) %>%
  # new globally consecutive id: increments whenever group_id_old changes
  dplyr::mutate(
    group_id = cumsum(dplyr::coalesce(group_id_old != dplyr::lag(group_id_old), TRUE))
  )
dat_ISIIS <- dat_ISIIS %>%
  mutate(
    shelf_loc_200 = case_when(
      rel_dist >= 0 ~ "on_shelf",
      rel_dist <  0 ~ "off_shelf",
      TRUE          ~ NA_character_
    ),
    shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf"))
  )






#1- mean slope per transect, but distinguishing inshore vs offshore, one panel per cruise
# Ensure factors are in the right order
# Violin plot

dat_ISIIS <- dat_ISIIS %>%
  mutate(
    Cruise        = factor(Cruise, levels = c("w22", "s22", "w23", "s23")),
    shelf_loc_200 = factor(shelf_loc_200, levels = c("on_shelf", "off_shelf")),
    Transect_ID   = factor(Transect_ID,
                           levels = c("GH", "OB", "CR", "CM", "NH", "HH", "RR"))
  )


#slope summary
p_slope <- ggplot(dat_ISIIS, aes(x = slope, y = Transect_ID, fill = shelf_loc_200)) +
  geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.3) +
  facet_wrap(~ Cruise, ncol = 4) +
  coord_cartesian(xlim = c(-2, 0)) +
  scale_fill_brewer(palette = "Set2") +
  
  # GH → RR order (top → bottom)
  scale_y_discrete(
    limits = rev(c("GH","OB","CR","CM","NH","HH","RR")),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    x = "Size-spectrum slope",
    y = "",
    fill = "Shelf location (200 m)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",                   # ⬅⬅ legend on top
    legend.direction = "horizontal",           # optional: lay out horizontally
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(size = 18),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1.2, "lines")
  )

p_slope

out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_slope_summary", ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = out_path, plot = p_slope, width = 9, height = 15, units = "in", dpi = 300, bg = "white")


#intercept summary
p_intercept <- ggplot(dat_ISIIS, aes(x = intercept, y = Transect_ID, fill = shelf_loc_200)) +
  geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.3) +
  facet_wrap(~ Cruise, ncol = 4) +
  coord_cartesian(xlim = c(0, 20)) +
  scale_fill_brewer(palette = "Set2") +
  
  # GH → RR order (top → bottom)
  scale_y_discrete(
    limits = rev(c("GH","OB","CR","CM","NH","HH","RR")),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    x = "Size-spectrum intercept",
    y = "",
    fill = "Shelf location (200 m)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",                   # ⬅⬅ legend on top
    legend.direction = "horizontal",           # optional: lay out horizontally
    strip.background = element_rect(fill = "gray90"),
    axis.text.y = element_text(size = 18),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(1.2, "lines")
  )

p_intercept

out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_intercept_summary", ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = out_path, plot = p_intercept, width = 9, height = 15, units = "in", dpi = 300, bg = "white")

#GLOBAL SLOPE RELATIONSHIPS:
df_list_combined <- list.files(sizeSpectraDir, pattern = "ISIIS_LISST_paired\\.csv$", full.names = TRUE, recursive = TRUE)
pairs <- map_dfr( df_list_combined, ~ read_csv( .x, show_col_types = FALSE))

mod <- lm(slope_ISIIS ~ slope_LISST, data = pairs)

# Extract coefficients
a  <- coef(mod)[2]   # slope
b  <- coef(mod)[1]   # intercept
r2 <- summary(mod)$r.squared

# Create formatted equation label
eq_label <- sprintf("y = %.3f x + %.3f\nR² = %.3f", a, b, r2)

# Build the plot
p <- ggplot(pairs, aes(x = slope_LISST, y = slope_ISIIS, color = chl)) +
  geom_point(alpha = 0.85) +
  
  # Linear regression line
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # Equation text added to bottom-right corner
  geom_text(
    aes(
      x = max(slope_LISST, na.rm = TRUE),
      y = min(slope_ISIIS, na.rm = TRUE),
      label = eq_label
    ),
    hjust = 1, vjust = 0,
    size = 5,
    inherit.aes = FALSE
  ) +
  
  scale_color_viridis_c(name = expression("Chl-a ("*mu*"g L"^{-1}*")")) +
  labs(x = "LISST slope", y = "ISIIS slope") +
  theme_classic() +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    aspect.ratio = 1
  ) +
  guides(color = guide_colorbar(
    barwidth  = unit(6, "cm"),
    barheight = unit(0.4, "cm")
  ))

p

out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_LISST_ISIIS_slopes.py", ".png"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

ggsave(filename = out_path, plot = p, width = 9, height = 10, units = "in", dpi = 300, bg = "white")


####. FUNCTION and code to explore relationships between environmental variables and coefficients
plot_two_panel_regression <- function(df, predictor_col) {

  
  #---------------------------------------
  # Prepare long-format dataframe
  #---------------------------------------
  df_plot <- df %>%
    dplyr::select(all_of(predictor_col), slope, intercept, season, shelf_loc_200) %>%
    pivot_longer(
      cols = c(slope, intercept),
      names_to = "metric",
      values_to = "value"
    )
  
  # Nicely formatted panel titles
  df_plot$metric <- factor(df_plot$metric,
                           levels = c("slope", "intercept"),
                           labels = c("Slope", "Intercept"))
  
  # Regression equation helper
  lm_eqn <- function(x, y) {
    model <- lm(y ~ x)
    a <- coef(model)[2]
    b <- coef(model)[1]
    r2 <- summary(model)$r.squared
    sprintf("y = %.3f x + %.3f\nR² = %.3f", a, b, r2)
  }
  
  #---------------------------------------
  # Core plot with GLOBAL regression only
  #---------------------------------------
  p_two <- ggplot(
    df_plot,
    aes_string(
      x = predictor_col,     # dynamic aesthetic
      y = "value",
      color = "season",
      shape = "shelf_loc_200"
    )
  ) +
    
    # points
    geom_point(alpha = 0.3, size = 2) +
    
    # ONE global regression per panel (black dashed)
    geom_smooth(
      data = df_plot,
      aes_string(
        x = predictor_col,
        y = "value"
      ),
      method = "lm",
      se = FALSE,
      linewidth = 1,
      color = "black",
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~ metric, ncol = 2, scales = "free_y") +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = predictor_col,
      y = "Statistic",
      color = "Season",
      shape = "Shelf location (200 m)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "right",
      strip.text = element_text(face = "bold")
    )
  
  #---------------------------------------
  # Regression equation labels
  #---------------------------------------
  eq_labels <- df_plot %>%
    group_by(metric) %>%
    summarise(
      eq = lm_eqn(
        x = .data[[predictor_col]],   # values inside this metric group
        y = value
      ),
      x = max(.data[[predictor_col]], na.rm = TRUE),
      y = min(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_two_final <- p_two +
    geom_text(
      data = eq_labels,
      aes(x = x, y = y, label = eq),
      hjust = 1, vjust = 0, size = 3.5,
      inherit.aes = FALSE
    )
  
  #save fig
  out_path <- file.path(output_dir_figs_NB, paste("SPECTRA_NBSS_coefficients_", predictor_col, ".png"))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = out_path, plot = p_two_final, width = 12, height = 9, units = "in", dpi = 300, bg = "white")
  
  print(p_two_final)
  return(p_two_final)
}

plot_regression_2x4 <- function(df, predictor_col,
                                output_dir = output_dir_figs_NB) {
  
  library(tidyverse)
  
  #-----------------------------------------
  # Sanity checks
  #-----------------------------------------
  if (!predictor_col %in% names(df)) {
    stop(paste("Column", predictor_col, "not found in dataframe"))
  }
  if (!all(c("season","shelf_loc_200","slope","intercept") %in% names(df))) {
    stop("df must contain: season, shelf_loc_200, slope, intercept")
  }
  
  #-----------------------------------------
  # Prep dataframe
  #-----------------------------------------
  df_plot <- df %>%
    dplyr::select(all_of(predictor_col), slope, intercept, season, shelf_loc_200) %>%
    pivot_longer(
      cols = c(slope, intercept),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric,
                      levels = c("slope", "intercept"),
                      labels = c("Slope", "Intercept")),
      season = factor(season),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("on_shelf","off_shelf"))
    )
  
  #-----------------------------------------
  # Regression helper: equation + R2 + p-value
  #-----------------------------------------
  lm_eqn <- function(x, y) {
    mod  <- lm(y ~ x)
    summ <- summary(mod)
    
    a    <- coef(mod)[2]            # slope
    b    <- coef(mod)[1]            # intercept
    r2   <- summ$r.squared          # R2
    pval <- summ$coefficients[2, 4] # p-value for slope term
    
    sprintf(
      "y = %.3f x + %.3f\nR² = %.3f\np = %.3f",
      a, b, r2, pval
    )
  }
  
  # Compute text locations and regression per facet
  stat_labels <- df_plot %>%
    group_by(metric, season, shelf_loc_200) %>%
    summarise(
      eq  = lm_eqn(.data[[predictor_col]], value),
      x   = max(.data[[predictor_col]], na.rm = TRUE),
      y   = min(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  #-----------------------------------------
  # Build plot
  #-----------------------------------------
  p_2x4 <- ggplot(df_plot,
                  aes_string(x = predictor_col, y = "value")) +
    
    geom_point(alpha = 0.35, size = 1.8,
               aes(color = season)) +
    
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.9,
      color = "black",
      linetype = "dashed"
    ) +
    
    facet_grid(
      rows = vars(metric),
      cols = vars(season, shelf_loc_200),
      scales = "free_y"
    ) +
    
    labs(
      x = predictor_col,
      y = "Value",
      title = paste("Slope and Intercept vs", predictor_col)
    ) +
    
    scale_color_brewer(palette = "Dark2") +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines")
    ) +
    
    geom_text(
      data = stat_labels,
      aes(label = eq),
      x     = Inf,
      y     = -Inf,
      hjust = 1.1,
      vjust = -0.4,
      size  = 3.5,
      inherit.aes = FALSE
    )
  
  #-----------------------------------------
  # Save figure
  #-----------------------------------------
  outfile <- file.path(
    output_dir,
    paste0("Regression_SeasonXShelf_", predictor_col, ".png")
  )
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = outfile,
    plot = p_2x4,
    width = 16, height = 8,
    dpi = 300, bg = "white"
  )
  
  message("Saved 2×4 regression plot to: ", outfile)
  print(p_2x4)
  return(p_2x4)
}


#filtering based on CTD comparison
dat_ISIIS <- dat_ISIIS %>%
  dplyr::filter(
    Salinity <= 40 | is.na(Salinity),
    between(Temperature, -60, 60) | is.na(Temperature),
    between(Oxygen_ml_l, 0, 10) | is.na(Oxygen_ml_l),
    pH <= 8.25 | is.na(pH)
    
    
  )


plot_regression_2x4(dat_ISIIS, "Temperature")
plot_regression_2x4(dat_ISIIS, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS, "pH")

#with subsetted data

dat_ISIIS_CR <- dat_ISIIS[dat_ISIIS$Transect_ID=='CR',]
dat_ISIIS_NH <- dat_ISIIS[dat_ISIIS$Transect_ID=='NH',]
dat_ISIIS_RR <- dat_ISIIS[dat_ISIIS$Transect_ID=='RR',]

plot_regression_2x4(dat_ISIIS_CR, "Temperature")
plot_regression_2x4(dat_ISIIS_CR, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS_CR, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS_CR, "pH")

plot_regression_2x4(dat_ISIIS_NH, "Temperature")
plot_regression_2x4(dat_ISIIS_NH, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS_NH, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS_NH, "pH")

plot_regression_2x4(dat_ISIIS_RR, "Temperature")
plot_regression_2x4(dat_ISIIS_RR, "chl_a_ug_l")
plot_regression_2x4(dat_ISIIS_RR, "Oxygen_ml_l")
plot_regression_2x4(dat_ISIIS_RR, "pH")



plot_slope_vs_intercept_2x4 <- function(df, output_dir = output_dir_figs_NB) {
  
  library(tidyverse)
  
  req_cols <- c("season", "shelf_loc_200", "slope", "intercept")
  if (!all(req_cols %in% names(df))) {
    stop(paste("df must contain:", paste(req_cols, collapse=", ")))
  }
  
  df_plot <- df %>%
    filter(
      !is.na(slope),
      !is.na(intercept),
      is.finite(slope),
      is.finite(intercept)
    ) %>%
    dplyr::select(season, shelf_loc_200, slope, intercept) %>%
    rename(
      intercept_val = intercept,
      slope_val     = slope
    ) %>%
    mutate(
      season = factor(season),
      shelf_loc_200 = factor(shelf_loc_200,
                             levels = c("on_shelf","off_shelf"))
    )
  
  lm_eqn <- function(x, y) {
    mod <- lm(y ~ x)
    sprintf(
      "y = %.3f x + %.3f\nR² = %.3f",
      coef(mod)[2], coef(mod)[1], summary(mod)$r.squared
    )
  }
  
  df_labels <- df_plot %>%
    group_by(season, shelf_loc_200) %>%
    summarise(
      eq = lm_eqn(intercept_val, slope_val),
      x = Inf,
      y = -Inf,
      .groups = "drop"
    )
  
  p <- ggplot(df_plot,
              aes(x = intercept_val, y = slope_val)) +
    geom_point(alpha = 0.35, size = 1.8,
               aes(color = season)) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.9,
      color = "black",
      linetype = "dashed"
    ) +
    facet_grid(
      rows = vars(season),
      cols = vars(shelf_loc_200),
      scales = "free"
    ) +
    labs(
      title = "Slope as a function of Intercept",
      x = "Intercept",
      y = "Slope"
    ) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines")
    ) +
    geom_text(
      data = df_labels,
      aes(label = eq, x = x, y = y),
      hjust = 1.1,
      vjust = -0.4,
      size = 3.5,
      inherit.aes = FALSE
    )
  
  out <- file.path(output_dir, "Slope_vs_Intercept_bySeasonShelf.png")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  
  ggsave(out, plot = p, width = 12, height = 8,
         bg = "white", dpi = 300)
  
  message("Saved slope-vs-intercept plot to: ", out)
  print(p)
  return(p)
}




plot_slope_vs_intercept_2x4(dat_ISIIS)