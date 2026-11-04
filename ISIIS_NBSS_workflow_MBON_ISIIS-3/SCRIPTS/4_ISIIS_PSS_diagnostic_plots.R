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
library(cowplot)

base_dir<- getwd()


## CODE FOR PLOTS BELOW


##  PART 1, size spectra

#create binned data if it doesn't exist
if (!dir.exists(file.path(base_dir,"figures"))) {
    dir.create(file.path(base_dir, "figures"), recursive = TRUE)
}

df_full<-read.csv(file.path(base_dir, "data", "size_spectra",'ISIIS_NBSS_by_community.csv'))
df_full_noThres<-read.csv(file.path(base_dir, "data", "size_spectra",'ISIIS_NBSS_by_community_no_threshold.csv'))
df_taxa<-  read.csv(file.path(base_dir, "data", "size_spectra",'ISIIS_NBSS_by_taxa.csv'))


# Calculate limits with buffer
x_min <- min(df_full$ECD_um_mid, na.rm = TRUE)
x_max <- max(df_full$ECD_um_mid, na.rm = TRUE) * 1.1
x_breaks <- sort(unique(c(200, 2000, 20000, signif(x_min, digits = 2))))

x_min_noThres <- min(df_full_noThres$ECD_um_mid, na.rm = TRUE)
x_max_noThres <- max(df_full_noThres$ECD_um_mid, na.rm = TRUE) * 1.1
x_breaks_noThres <- sort(unique(c(200, 2000, 20000, signif(x_min_noThres, digits = 2))))

# Thresholded plot
size_spectra_graph_community <- ggplot() +
  geom_line(data = df_full, aes(x = ECD_um_mid, y = NB, group = group_id, color = depth_mid), alpha = 0.5, size = 0.05) +
  geom_point(data = df_full, aes(x = ECD_um_mid, y = NB, color = depth_mid), size = 1, alpha = 0.5, shape = 'o') +
  geom_line(data = df_full %>%
              mutate(ECD_um_mid = as.factor(ECD_um_mid)) %>%
              group_by(shelf_loc_200, ECD_um_mid) %>%
              summarise(NB = median(NB), .groups = "drop") %>%
              mutate(ECD_um_mid = as.numeric(as.character(ECD_um_mid))),
            aes(x = ECD_um_mid, y = NB, group = shelf_loc_200),
            color = "black", alpha = 1, size = 0.5) +
  scale_x_log10(limits = c(x_min, x_max), breaks = x_breaks) +  
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_color_viridis_c(option = "D", direction = -1) +
  facet_wrap(~ shelf_loc_200, scales = "fixed") +
  labs(
    y = expression(Normalized~Biovolume~(mu*l~L^{-1}~u*L^{-1})),
    x = expression(Equivalent~spherical~diameter~(mu*m)),
    color = "Depth (m)"
  ) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.text.x = element_text(family = "Times", size = 15),
    axis.title.x = element_text(family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 15),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey")
  )
ggsave(size_spectra_graph_community, file = file.path(base_dir, 'figures', 'ISIIS_size_spectra_shelf_depth.png'), width = 8, height = 6)

# No-threshold plot
size_spectra_graph_community_noThres <- ggplot() +
  geom_line(data = df_full_noThres, aes(x = ECD_um_mid, y = NB, group = group_id, color = depth_mid), alpha = 0.5, size = 0.05) +
  geom_point(data = df_full_noThres, aes(x = ECD_um_mid, y = NB, color = depth_mid), size = 1, alpha = 0.5, shape = 'o') +
  geom_line(data = df_full_noThres %>%
              mutate(ECD_um_mid = as.factor(ECD_um_mid)) %>%
              group_by(shelf_loc_200, ECD_um_mid) %>%
              summarise(NB = median(NB), .groups = "drop") %>%
              mutate(ECD_um_mid = as.numeric(as.character(ECD_um_mid))),
            aes(x = ECD_um_mid, y = NB, group = shelf_loc_200),
            color = "black", alpha = 1, size = 0.5) +
  scale_x_log10(limits = c(x_min_noThres, x_max_noThres), breaks = x_breaks_noThres) +  
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_color_viridis_c(option = "D", direction = -1) +
  facet_wrap(~ shelf_loc_200, scales = "fixed") +
  labs(
    y = expression(Normalized~Biovolume~(mu*l~L^{-1}~u*L^{-1})),
    x = expression(Equivalent~spherical~diameter~(mu*m)),
    color = "Depth (m)"
  ) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.text.x = element_text(family = "Times", size = 15),
    axis.title.x = element_text(family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 15),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey")
  )
ggsave(size_spectra_graph_community_noThres, file = file.path(base_dir, 'figures', 'ISIIS_size_spectra_shelf_depth_noThres.png'), width = 8, height = 6)


  ## PART 2 coefficients by shelf location

slope_graph <- ggplot(data = df_full, aes(x = shelf_loc_200, y = slope)) +
  geom_violin(aes(fill = shelf_loc_200, color = shelf_loc_200)) +
  geom_crossbar(stat = "summary", fun = mean, fatten = 2, width = .5) +
  ylim(-2, 0) +
  labs(x = '', y = expression(NBSS ~ slope ~ (L^{-1}~mu*L^{-1}))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(family = "Times", size = 15),
    axis.title.x = element_text(family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 15),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey")
  )

# INTERCEPT plot by shelf location
intercept_graph <- ggplot(data = df_full, aes(x = shelf_loc_200, y = intercept)) +
  geom_violin(aes(fill = shelf_loc_200, color = shelf_loc_200)) +
  geom_crossbar(stat = "summary", fun = mean, fatten = 2, width = .5) +
  labs(x = '', y = expression(NBSS ~ intercept ~ (mu*L~L^{-1}~mu*L^{-1}))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(family = "Times", size = 15),
    axis.title.x = element_text(family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 15),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey")
  )

# R2 plot by shelf location
R2_graph <- ggplot(data = df_full, aes(x = shelf_loc_200, y = R2)) +
  geom_violin(aes(fill = shelf_loc_200, color = shelf_loc_200)) +
  labs(x = '', y = expression(Determination ~ coefficient ~ (R^2))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(family = "Times", size = 15),
    axis.title.x = element_text(family = "Times", size = 15),
    axis.text.y = element_text(family = "Times", size = 15),
    axis.title.y = element_text(family = "Times", size = 15),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey")
  )

coeffs_plots <- plot_grid(slope_graph, intercept_graph, R2_graph, ncol = 3, labels = "AUTO")
print(coeffs_plots)

ggsave(
  coeffs_plots,
  filename = file.path(base_dir, "figures", "ISIIS_coefficients_size_spectra_shelfLoc.png"),
  width = 16, height = 4, dpi = 300
)

## PART 4: slope as a function of temperature and depth

# Panel A: Temperature vs. Slope
p1 <- ggplot(df_full, aes(x = Temp, y = slope)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Temperature (°C)", y = expression(NBSS ~ slope)) +
  theme_minimal(base_family = "Times") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "grey")
  )

# Panel B: Depth vs. Slope
p2 <- ggplot(df_full, aes(x = depth_mid, y = slope)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Depth (m)", y = expression(NBSS ~ slope)) +
  theme_minimal(base_family = "Times") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "grey")
  )

# Combine panels
slope_env_plot <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)

# Save
ggsave(
  filename = file.path(base_dir, "figures", "Slope_vs_Temp_and_Depth.png"),
  plot = slope_env_plot,
  width = 12, height = 5, dpi = 300
)


 ## PART 5: size spectra per taxa, coded by inshore/offshore

##size spectra graphs for major taxa
# Build list of plots
plot_list <- list()
n <- 0

for (i in unique(df_taxa$predicted_taxon)) {
  n <- n + 1
  df_sub <- df_taxa %>% filter(predicted_taxon == i)

  # Summary curve per shelf location
  summary_df <- df_sub %>%
    mutate(ECD_um_mid = as.factor(ECD_um_mid)) %>%
    group_by(shelf_loc_200, ECD_um_mid) %>%
    summarise(NB = median(NB, na.rm = TRUE), .groups = "drop") %>%
    mutate(ECD_um_mid = as.numeric(as.character(ECD_um_mid)))

  # Plot
  p <- ggplot() +
    geom_line(data = df_sub,
              aes(x = ECD_um_mid, y = NB, group = group_id, color = shelf_loc_200),
              alpha = 0.7, size = 0.1) +
    geom_point(data = df_sub,
               aes(x = ECD_um_mid, y = NB, color = shelf_loc_200),
               size = 0.5, alpha = 0.7, shape = 'o') +
    geom_line(data = summary_df,
              aes(x = ECD_um_mid, y = NB, group = shelf_loc_200, color = shelf_loc_200),
              alpha = 1, size = 0.5) +
    scale_x_log10(limits = c(50, 35000), breaks = c(200, 2000, 20000)) +
    scale_y_log10(labels = label_log(digits = 2)) +
    ggtitle(i) +
    labs(y = expression(Normalized~Biovolume~(mu*l~L^{-1}~u*L^{-1})),
         x = expression(Equivalent~spherical~diameter~(mu*m))) +
    theme(axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          axis.text.x = element_text(family = "Times", size = 15),
          axis.title.x = element_text(family = "Times", size = 15),
          axis.text.y = element_text(family = "Times", size = 15),
          axis.title.y = element_text(family = "Times", size = 15),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"))

  plot_list[[paste0("plot", n)]] <- p
}

# Combine all into one figure
all_taxa_plots <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 6)

# Save figure using base_dir
output_path <- file.path(base_dir, "figures", "taxa_ISIIS_size_spectra_shelfLoc.png")

ggsave(filename = output_path, plot = all_taxa_plots, width = 19, height = 25, dpi = 300)

message("✅ Saved taxa-wise size spectra figure to: ", output_path) 



