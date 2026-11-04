#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
})
setwd('~/SCRIPTS_backups/ISIIS_NBSS_workflow_SPECTRA')
# ============================================================
# Plot one representative/irregular plankton NBSS
#
# Selection:
#   1. Combine all *by_community.csv plankton files
#   2. Assign globally unique group_id values
#   3. Count unique ECD_um_mid size classes per group_id
#   4. Fit logNB ~ log10(ECD_um_mid) for every group_id
#   5. Calculate residual_logNB
#   6. Select the group_id with:
#        a) the most unique size classes
#        b) among ties, the greatest residual variability
#   7. Plot ONLY that group_id, including observations
#      and its linear fit
# ============================================================

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_dir <- getwd()

sizeSpectraDirPlankton <- file.path(base_dir, "data", "size_spectra_plankton_5-10m-bins_thresholded")

output_dir <- file.path(base_dir, "figures", "single_plankton_NBSS")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Find all plankton by_community files
# ------------------------------------------------------------

df_list_plankton <- list.files(sizeSpectraDirPlankton, pattern = "by_community\\.csv$", full.names = TRUE, recursive = TRUE)

if (length(df_list_plankton) == 0) stop("No *by_community.csv files found under:\n", sizeSpectraDirPlankton)

message("Found ", length(df_list_plankton), " by_community.csv files.")

# ------------------------------------------------------------
# Read and combine all files
#
# group_id is only unique within an individual file.
# Retain source_file and group_id_old first.
# ------------------------------------------------------------

dat_plankton <- map_dfr(df_list_plankton, function(f) {
  read_csv(f, show_col_types = FALSE, col_select = c("group_id", "ECD_um_mid", "logNB", dplyr::any_of(c("slope", "intercept", "Cruise", "Transect_ID", "depth_mid", "time_Pacific", "Latitude", "Longitude", "residual_logNB")))) %>% mutate(source_file = f)
})

# ------------------------------------------------------------
# Check required columns
# ------------------------------------------------------------

required_cols <- c("group_id", "ECD_um_mid", "logNB")

missing_cols <- required_cols[!required_cols %in% names(dat_plankton)]

if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# ------------------------------------------------------------
# Assign globally unique group_id
#
# Unique sample identity:
# source_file + original group_id
# ------------------------------------------------------------

dat_plankton <- dat_plankton %>% rename(group_id_old = group_id) %>% mutate(ECD_um_mid = as.numeric(ECD_um_mid), logNB = as.numeric(logNB)) %>% mutate(sample_key = paste(source_file, group_id_old, sep = "::")) %>% mutate(group_id = match(sample_key, unique(sample_key)))

message("Total number of globally unique plankton samples: ", n_distinct(dat_plankton$group_id))

# ------------------------------------------------------------
# Retain valid NBSS observations
# ------------------------------------------------------------

dat_valid <- dat_plankton %>% filter(is.finite(ECD_um_mid), ECD_um_mid > 0, is.finite(logNB))

# ------------------------------------------------------------
# Fit each individual NBSS and calculate residuals
#
# NBSS relationship:
#
# log10(NB) = slope * log10(size) + intercept
# ------------------------------------------------------------

fit_one_spectrum <- function(df) {
  if (nrow(df) < 2 || n_distinct(df$ECD_um_mid) < 2) {
    df$residual_logNB <- NA_real_
    df$fitted_logNB <- NA_real_
    return(df)
  }
  mod <- lm(logNB ~ log10(ECD_um_mid), data = df)
  df <- df %>% mutate(fitted_logNB = fitted(mod), residual_logNB = residuals(mod))
  return(df)
}

dat_with_residuals <- dat_valid %>% group_by(group_id) %>% group_modify(~ fit_one_spectrum(.x)) %>% ungroup()

# ------------------------------------------------------------
# Summarize every spectrum
#
# residual_sd = variability in residual_logNB
# ------------------------------------------------------------

spectrum_summary <- dat_with_residuals %>% group_by(group_id) %>% summarise(n_size_classes = n_distinct(ECD_um_mid[is.finite(ECD_um_mid)]), residual_sd = if (sum(is.finite(residual_logNB)) >= 2) sd(residual_logNB, na.rm = TRUE) else NA_real_, residual_variance = if (sum(is.finite(residual_logNB)) >= 2) var(residual_logNB, na.rm = TRUE) else NA_real_, group_id_old = first(group_id_old), source_file = first(source_file), Cruise = if ("Cruise" %in% names(cur_data())) first(Cruise) else NA_character_, Transect_ID = if ("Transect_ID" %in% names(cur_data())) first(Transect_ID) else NA_character_, depth_mid = if ("depth_mid" %in% names(cur_data())) first(depth_mid) else NA_real_, time_Pacific = if ("time_Pacific" %in% names(cur_data())) first(time_Pacific) else NA_character_, .groups = "drop")

# ------------------------------------------------------------
# Select ONE spectrum
#
# PRIMARY:
# maximum number of unique size classes
#
# SECONDARY:
# maximum SD of residual_logNB
# ------------------------------------------------------------

selected_summary <- spectrum_summary %>% filter(is.finite(residual_sd)) %>% arrange(desc(n_size_classes), desc(residual_sd)) %>% slice(1)

if (nrow(selected_summary) == 0) stop("Could not identify a spectrum with sufficient data to calculate residual variability.")

selected_group_id <- selected_summary$group_id[[1]]

# ------------------------------------------------------------
# Extract ONLY selected group_id
# ------------------------------------------------------------

selected_NBSS <- dat_with_residuals %>% filter(group_id == selected_group_id) %>% arrange(ECD_um_mid)

# ------------------------------------------------------------
# Refit selected spectrum explicitly
# ------------------------------------------------------------

selected_model <- lm(logNB ~ log10(ECD_um_mid), data = selected_NBSS)

model_coef <- coef(selected_model)

selected_intercept <- unname(model_coef["(Intercept)"])

selected_slope <- unname(model_coef["log10(ECD_um_mid)"])

selected_r2 <- summary(selected_model)$r.squared

# ------------------------------------------------------------
# Calculate fitted values at each observed size class
# ------------------------------------------------------------

selected_NBSS <- selected_NBSS %>% arrange(ECD_um_mid) %>% mutate(fit_logNB = predict(selected_model, newdata = .), residual = logNB - fit_logNB, logECD = log10(ECD_um_mid))

# ------------------------------------------------------------
# Build polygons between observed NBSS and fitted NBSS
#
# Red  = observed data > fit
# Blue = observed data < fit
#
# Crossings are interpolated in log10(ECD) space.
# ------------------------------------------------------------

polygon_list <- list()

polygon_id <- 0L

for (i in seq_len(nrow(selected_NBSS) - 1)) {
  
  x1 <- selected_NBSS$logECD[i]
  x2 <- selected_NBSS$logECD[i + 1]
  
  y1_data <- selected_NBSS$logNB[i]
  y2_data <- selected_NBSS$logNB[i + 1]
  
  y1_fit <- selected_NBSS$fit_logNB[i]
  y2_fit <- selected_NBSS$fit_logNB[i + 1]
  
  r1 <- y1_data - y1_fit
  r2 <- y2_data - y2_fit
  
  if (r1 == 0 && r2 == 0) next
  
  if (r1 * r2 >= 0) {
    
    polygon_id <- polygon_id + 1L
    
    fill_type <- ifelse(mean(c(r1, r2)) >= 0, "Data > fit", "Data < fit")
    
    polygon_list[[polygon_id]] <- tibble(polygon_id = polygon_id, fill_type = fill_type, logECD = c(x1, x2, x2, x1), logNB = c(y1_data, y2_data, y2_fit, y1_fit))
    
  } else {
    
    frac <- abs(r1) / (abs(r1) + abs(r2))
    
    x_cross <- x1 + frac * (x2 - x1)
    
    y_data_cross <- y1_data + frac * (y2_data - y1_data)
    
    y_fit_cross <- y1_fit + frac * (y2_fit - y1_fit)
    
    y_cross <- mean(c(y_data_cross, y_fit_cross))
    
    polygon_id <- polygon_id + 1L
    
    fill_type_1 <- ifelse(r1 > 0, "Data > fit", "Data < fit")
    
    polygon_list[[polygon_id]] <- tibble(polygon_id = polygon_id, fill_type = fill_type_1, logECD = c(x1, x_cross, x_cross, x1), logNB = c(y1_data, y_cross, y_cross, y1_fit))
    
    polygon_id <- polygon_id + 1L
    
    fill_type_2 <- ifelse(r2 > 0, "Data > fit", "Data < fit")
    
    polygon_list[[polygon_id]] <- tibble(polygon_id = polygon_id, fill_type = fill_type_2, logECD = c(x_cross, x2, x2, x_cross), logNB = c(y_cross, y2_data, y2_fit, y_cross))
  }
}

polygon_dat <- bind_rows(polygon_list) %>% mutate(ECD_um_mid = 10^logECD)

# ------------------------------------------------------------
# Print selected sample information
# ------------------------------------------------------------

message("")

message("====================================================")

message("SELECTED PLANKTON NBSS")

message("====================================================")

message("New global group_id: ", selected_group_id)

message("Original group_id: ", selected_summary$group_id_old)

message("Source file: ", selected_summary$source_file)

message("Unique ECD size classes: ", selected_summary$n_size_classes)

message("SD residual_logNB: ", signif(selected_summary$residual_sd, 5))

message("Slope: ", signif(selected_slope, 5))

message("Intercept: ", signif(selected_intercept, 5))

message("R2: ", signif(selected_r2, 5))

message("====================================================")

message("")

# ------------------------------------------------------------
# Generate fitted line
# ------------------------------------------------------------

fit_line <- tibble(ECD_um_mid = 10^seq(log10(min(selected_NBSS$ECD_um_mid, na.rm = TRUE)), log10(max(selected_NBSS$ECD_um_mid, na.rm = TRUE)), length.out = 500)) %>% mutate(logNB = predict(selected_model, newdata = .))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(selected_NBSS, aes(x = ECD_um_mid, y = logNB)) +
  geom_polygon(data = polygon_dat, aes(x = ECD_um_mid, y = logNB, group = polygon_id, fill = fill_type), inherit.aes = FALSE, alpha = 0.30, color = NA) +
  geom_line(aes(group = 1), linetype = "dashed", linewidth = 0.8, color = "black") +
  geom_point(shape = 21, size = 3.5, stroke = 0.8, fill = "white", color = "black") +
  geom_line(data = fit_line, aes(x = ECD_um_mid, y = logNB), inherit.aes = FALSE, linewidth = 1.1, color = "black") +
  scale_fill_manual(values = c("Data > fit" = "red", "Data < fit" = "blue")) +
  scale_x_log10() +
  labs(x = expression(ECD~(mu*m)), y = expression(log[10]~Normalized~biovolume), fill = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"), axis.text = element_text(size = 15))

# ------------------------------------------------------------
# Save plot
# ------------------------------------------------------------

output_file <- file.path(output_dir, paste0("selected_plankton_NBSS_group_", selected_group_id, ".png"))

ggsave(filename = output_file, plot = p, width = 7, height = 6, dpi = 300, bg = "white")

# ------------------------------------------------------------
# Save selected spectrum and ranking table
# ------------------------------------------------------------

write_csv(selected_NBSS, file.path(output_dir, paste0("selected_plankton_NBSS_group_", selected_group_id, ".csv")))

write_csv(spectrum_summary %>% arrange(desc(n_size_classes), desc(residual_sd)), file.path(output_dir, "plankton_NBSS_selection_ranking.csv"))

message("Saved figure: ", output_file)