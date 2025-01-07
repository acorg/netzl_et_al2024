# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(ggplot2)

working_dir <- getwd()

target_map <- "omicron_map_preprocessed_wSAVE_lessWD"

map_dir <- file.path(working_dir, "data", "maps")
metadata_dir <- file.path(working_dir,'data', 'metadata')
save_dir <- file.path(working_dir, "data","maps", target_map, "validation")
figure_dir <- file.path(working_dir, "figures","maps", target_map, "validation")
suppressWarnings(dir.create(figure_dir, recursive = T))


results <- readRDS(file.path(save_dir, "cross_validation.rds" ))

# Read the map
map <- read.acmap(file.path(map_dir, paste0(target_map, '_conv_ag_sub.ace')))

agFillScale <- function(map) {
  fills <- agFill(map)
  names(fills) <- agNames(map)
  fills
}


# Set detectable results subset
detectable_results <- results %>%
  filter(!is.infinite(predicted_logtiter))
detectable_rmse <- sqrt(mean(detectable_results$residual^2))

# do histogram of pred - measured
mean <- round(mean(detectable_results$residual, na.rm = T),2)
sd <- round(sd(scale(detectable_results$residual, scale = F), na.rm = T), 2)

hist_diff <- ggplot(detectable_results) +
  geom_histogram(aes(x = residual), fill = "grey50", alpha = 0.8, bins = 100) +
  xlim(c(-15, 15)) +
  geom_vline(xintercept = mean, linetype = "dashed") +
  labs(x= "Measured - predicted log2 titers", y = "Count", title = paste0("Mean = ", mean, "; SD = ", sd)) +
  theme_bw()

ggsave(plot = hist_diff, filename = file.path(figure_dir, "histogram_residuals.png"), width = 5, height = 4, dpi = 300)

do_hist_single_sr_group <- function(data, sr_group_spec, color) {
  
  data %>%
    filter(sr_group == sr_group_spec) -> data_sub
  
  mean <- round(mean(data_sub$residual, na.rm = T),2)
  sd <- round(sd(scale(data_sub$residual, scale = F), na.rm = T), 2)
  
  data_sub %>%
    ggplot() +
    geom_histogram(aes(x = residual), fill = color, alpha = 0.8, bins = 100) +
    geom_vline(xintercept = mean, linetype = "dashed") +
    xlim(c(0, 300)) +
    ylim(c(0, 3000)) +
    labs(x= "Measured - predicted log2 titers", y = "Count", title = paste0("Mean = ", mean, "; SD = ", sd)) +
    theme_bw() -> plot
  
  return(plot)
}

stop()
detectable_results %>%
  arrange(abs(residual)) %>%
  filter(abs(residual) >= 5) -> high_res

View(titerTable(map)[,high_res$sr_num])

ag_pretty <- c("D614G" = "D614G",
               "WT" = "WT",
               "B.1.1.7" = "Alpha",
               "B.1.351" = "Beta",
               "P.1" = "Gamma",
               "B.1.617.2" = "Delta",
               "BA.1" = "BA.1",
               "BA.2" = "BA.2")
# have to change it to include new antigens and sera
detectable_results$ag_pretty <- factor(ag_pretty[as.character(detectable_results$ag_name)], levels = ag_pretty)
detectable_results$sr_pretty <- detectable_results$sr_group #factor(sr_pretty[as.character(detectable_results$sr_group),], levels = sr_pretty$val)

    # Antigen and serum group tab
        detectable_results %>%
          ggplot(
            aes(
              x = predicted_logtiter,
              y = measured_logtiter,
              color = ag_name
            )
          ) +
          labs(x = "Predicted log2 titer",
               y = "Measured log2 titer") +
          # geom_smooth() +
          geom_point(
            alpha = 0.1
          ) +
          geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed"
          ) +
          scale_color_manual(
            values = agFillScale(map)
          ) +
          xlim(c(-10,10))+
          ylim(c(-10, 10))+
          facet_grid(
            cols = vars(sr_pretty),
            rows = vars(ag_pretty)
          ) +
          theme_bw() +
          coord_fixed() +
          theme(legend.position = "none",
               strip.text.x = element_text(size = 6),
               strip.text.y = element_text(size = 6))-> gp


ggsave(plot = gp, filename = file.path(figure_dir, "scatter_pred_vs_measured.png"), width = 8, height = 8, dpi = 300)
        
