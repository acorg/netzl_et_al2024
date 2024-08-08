#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(ablandscapes)
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

# ----------------------- set directories and load functions
working_dir <- getwd()

code_dir <- file.path(working_dir, "code")
data_dir <- file.path(working_dir, "data")
map_dir <- file.path(data_dir, "maps")
figure_dir_b <- file.path(working_dir, "figures", "landscapes", "gmt_landscapes")


source(file.path(code_dir, "utility", 'map_longinfo.R'))
source(file.path(code_dir, "utility", 'landscape_functions_util.R'))

# Read in data from Wilks et al. for slope comparison
wilks_full <- readRDS("data/maps/slope_comparison_wt.rds")
wilks_df <- wilks_full$slope_factors

sr_group_names <- c("2x mRNA-1273" = "2x Vax ",
                    "3x mRNA-1273 BD01" = "2x Vax ",
                    "3x mRNA-1273 BD29" = "3x Vax ",
                    "3x mRNA-1273 (6 month)" = "3x Vax ")

sr_group_times <- c("2x mRNA-1273" = 1,
                    "3x mRNA-1273 BD01" = 6,
                    "3x mRNA-1273 BD29" = 1,
                    "3x mRNA-1273 (6 month)" = 6)

wilks_df %>%
  mutate(month = sr_group_times[sr_group],
         sr_group = sr_group_names[sr_group]) -> wilks_slopes

# ---------------------- do fitting per map
ags_to_exclude <- c("BA.2.75")

map_names <- list('omicron_map_preprocessed_wSAVE_lessWD_conv_ag_sub.ace')


sr_colors <- read.csv(file = file.path(data_dir,"metadata", "map-colors.csv"), header = TRUE, stringsAsFactors = FALSE, sep = ",", row.names = "Variable")


for(map_name in map_names) {
  
  figure_dir <- file.path(figure_dir_b, gsub(".ace", "", map_name))
  suppressWarnings(dir.create(figure_dir, recursive = T))
  
  # set name for fit
  fit_name <- gsub("_conv_ag_sub", "_conv_ag_sub_since_exposure", map_name)
  fit_name <- gsub(".ace", ".rds", map_name)
  
  # read base map
  map <- read.acmap(file.path(map_dir, map_name))
  
  map <- removeAntigens(map, ags_to_exclude)
  lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
  
  ags_to_fit_lndscp <- agNames(map)
  # read in since exposure data
  
  data_since_exp <- read.csv(file.path(data_dir, "titer_tables", "omicron_folddrops_preprocessed_wSAVE_lessWD_titers_since_exposure.csv")) 
 
  # group exposures by month: 0.5, 1, 3, 6, 9, 12
  data_since_exp$mean_month <- sapply(data_since_exp$mean_days, function(x){
    
    m <- round(x/30, 2)
    
    if(m > (9+ 3/2)){
      12
    } else if(m > (6 + 3/2)){
      9
    } else if(m > (3 + 3/2)){
      6
    } else if(m > (1 + 2/2)){
      3
    } else if(m > (0.5 + 0.5/2)){
      1
    } else {
      0.5
    }
    
  })
  
  data_since_exp$sr_group <- paste(data_since_exp$standardise_encounters, " - ", data_since_exp$mean_month)
  
  # set rownames for antigens
  ag_names <- c("Delta" = "B.1.617.2", "Alpha" = "B.1.1.7", "Beta" = "B.1.351", "Gamma" = "P.1")
  
  data_since_exp$ag_name <- sapply(data_since_exp$OmicronVariant, function(x){
    if(x %in% names(ag_names)){
      ag_names[x]
    } else {
      x
    }
  })
  
  data_since_exp <- data_since_exp %>%
    select(!Comparator.antigen) %>%
    unique()
  
  # subset the map target sera and ags
  data_since_exp %>%
     mutate(titer = 2^Log2Omi*10,
            sr_name = paste(Study , Sera.details.long, standardised_pseudo, mean_days, sr_group, sep = "_")) %>%
     filter(ag_name %in% ags_to_fit_lndscp) %>%
     select(titer, ag_name, sr_name, sr_group) %>%
     unique() -> titerdata
 
  # subset to sera with at least three titrations
  titerdata %>%
    group_by(sr_name) %>%
    mutate(n_detectable = length(titer[!(titer %in% c("*", NA))])) %>%
    ungroup() %>%
    filter(n_detectable > 1) %>%
    group_by(
      sr_group
    ) -> titerdata
  
  titerdata %>%
    group_map(
      get_titertable
    ) -> titertables
  
  lndscp_fits <- lapply(
    titertables,
    
    function(titertable) {
      
       if(!is.null(dim(titertable))){
         ags_t_table <- colnames(titertable)
    
         ablandscape.fit(
           titers = titertable[,ags_t_table, drop = TRUE], 
           bandwidth = 1,
           degree = 1,
           method = "cone",
           error.sd = 1.1,
           acmap = map,
           control = list(
             optimise.cone.slope = TRUE,
             max.titer.possible = 20
           )
         )
       }
      
      
      
    }
  )
  
  titertables_groups <- group_data(titerdata)
 
  slope_df <- list()
  for(i in 1:length(titertables_groups$sr_group)){
    if(!is.null(lndscp_fits[[i]])){
      log_table <- lndscp_fits[[i]]$logtiters
      
      n_measurable <- nrow(log_table)*length(unique(titerdata$ag_name))
      n_measured <- length(log_table[!is.na(log_table)])
      
      slope_df[[paste0(i)]] <- data.frame("sr_group_b" = titertables_groups$sr_group[i],
                                          "sr_group" = strsplit(titertables_groups$sr_group[i], " - ")[[1]][1],
                                          "month" = strsplit(titertables_groups$sr_group[i], " - ")[[1]][2],
                                          "slope" = lndscp_fits[[i]]$cone$cone_slope,
                                          "n" = nrow(log_table),
                                          "fraction_measured" = n_measured/n_measurable)
      
    }
    
  }
  
  slope_df <- do.call(rbind, slope_df)
  
  slope_df %>%
    group_by(sr_group) %>%
    count() %>% 
    filter(n > 1) %>%
    pull(sr_group) -> time_sr_groups
  
  
  # do subplot with only 2x, 3x Vax
  slope_df$Data <- "Collected data"
  wilks_slopes$Data <- "Wilks et al. (2023)"
  
  sr_levels <- c("2x Vax ", "3x Vax ")
  slope_df %>%
    filter(sr_group %in% sr_levels) %>%
    mutate(month = as.numeric(month),
           confidence = n*fraction_measured,
           "Scaled confidence" = confidence/max(confidence),
           sr_group = factor(sr_group, levels =sr_levels)) %>%
    ggplot(aes(x = month, y = slope, group = sr_group, shape = Data)) + 
    geom_pointrange(data = wilks_slopes, aes(y = estimate, ymin = lower, ymax = upper), color = "grey60") +
    geom_line(data = wilks_slopes, aes(y = estimate), color = "grey60") +
    geom_line(color = "skyblue") +
    geom_point(aes(color = fraction_measured, size = `Scaled confidence`)) +
    geom_text(aes(label = n, y = ifelse(month == 0.5, 1.15, 1.2)), size = 3) +
    geom_text(label = "n = ", y = 1.2, x = -0.75, size = 3) +
    scale_x_continuous(name = "Time since exposure (months)",
                       breaks = seq(0, 10, 2),
                       limits = c(-1.2, 10)) +
    scale_y_continuous(name = "Fitted landscape slope",
                       limits = c(0, 1.25),
                       breaks = seq(0, 1.2, 0.2)) +
    scale_shape_manual(name = "Data",
                         values = c("Wilks et al. (2023)" = 1, "Collected data" = 19)) +
    scale_size(name = "Scaled data amount") +
    scale_color_binned(name = "Fraction of titrated variants",
                       low = "skyblue",
                       high = "darkblue") +
    facet_wrap(~sr_group, ncol = length(time_sr_groups)) + 
    theme_bw() + 
    theme(strip.background.x = element_blank(),
          strip.text.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) -> p
 
  ggsave(file.path(figure_dir, "main_slopes_since_exposure.png"), plot = p, units = "in", dpi = 300, width = 6, height = 4)
  
  
  sr_levels <- c("WT conv ", "2x Vax ", "3x Vax ", "Inf + Vax ", "Vax + BA.1 ")
  
  
  # now compare calculated and fitted gmt
  # Add impulses
  titerdata %>%
    group_by(
      sr_group,
      ag_name
    ) %>%
    summarize(gmt = mean(log2(as.numeric(titer)/10), na.rm = TRUE)) -> gmt_data
  
  null_landscpaes <- sapply(lndscp_fits, is.null)
  
  non_null_lndscps <- lndscp_fits[!null_landscpaes]
  
  sr_levels_time <- sapply(sr_levels, function(x){
    paste0(x, " -  ", c(0.5, 1, 3, 6, 9, 12))
  })
  
  comb <- combine_landscape_and_calculated_gmt(lndscp_fits = non_null_lndscps, gmt_data = gmt_data, sr_group_fields = 5)
  
  comb <- comb %>%
    mutate(sr_group = factor(sr_group, levels = sr_levels_time)) %>%
    filter(sr_group %in% unique(slope_df$sr_group_b))
  
  plot_lndscp_calculated_gmt_lineplot(comb) + 
    theme(legend.position = "top") -> p_gmts
  
  # do diff to GMT, not residuals
  comb %>%
    group_by(ag_name, sr_group) %>%
    mutate(gmt_diff = logtiter[Data == "Calculated GMT"] - logtiter,
           sr_group = factor(sr_group, levels = sr_levels_time)) %>%
    filter(!is.na(sr_group)) -> gmt_diff
  
  plot_gmt_diff(gmt_diff) -> p_diff
  
  p_gmts/p_diff + plot_annotation(tag_levels = 'A') -> p_comb
  p_comb
  ggsave(file.path(figure_dir, "landscapes_gmt_comb.png"), plot = p_comb, units = "in", dpi = 300, width = 8, height = 12)
  

  
  
}

