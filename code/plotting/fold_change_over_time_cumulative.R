# create forest plots for the different tables and variants
# setup page and load metadata
rm(list = ls())

library(meantiter)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)
library(gtable)
library(patchwork)
library(ggpubr)


# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "over_time")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
#----------------------------------------------- set path to save -----------------------------------------------
fileext <- "png"

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path(utility_dir, "plot_over_time_functions.R"))

#-------------------------------------  SET TABLE NAME
table_names = list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')

#for (table_name in table_names){
for (table_name in table_names){
  
  
  # create path to save for each table
  tab_name <- strsplit(table_name, "\\.")[[1]][1]
 
  path_to_save <- file.path(figures_dir, tab_name)
  
  # create subfolder for boxplots
  suppressWarnings(dir.create(path_to_save, recursive = TRUE))
  
  
  table_path <- file.path(google_sheets_dir,table_name)
  
  
#----------------- load and prepare data
  forest_data <- read.csv(table_path)
  
  
  forest_data <- format_table_forest_plots(forest_data) %>%
    filter(Webplotdigitizer != "y")

  forest_data <- standardise_time_format(forest_data)
  
  forest_data$standardise_time[forest_data$Study == "HKU"] <- as.Date("2021/12/12")
  forest_data$standardise_time[forest_data$Study == "Krumbholz"] <- as.Date("2022/03/29")
  forest_data %>%
    mutate(Study = gsub("Suthar \\(SAVE as publihed\\)", "Suthar", Study),
           Time = standardise_time,
           'Assay type' = standardised_pseudo,
           'Cell type' = standardised_cell) %>%
    select(Study, Sourcelink, Time, `Assay type`, `Cell type`) %>%
    arrange(Time) %>%
    unique() -> main_table
  
  write.csv(main_table, "data/summary_tables/main_table.csv")  
  
  forest_data %>%
    filter(!is.na(standardise_time)) %>%
    select(all_of(c("standardise_time", "first_record", "time_since_first", "Comparator antigen", "OmicronVariant", "Sera_details_no_time", "Log2HAg", "Log2Omi",
           "Standardised_sera_names", "standardise_encounters", "vacc_type_het", "vacc_type","vaccine_manufacturer", "standardised_assay", "standardised_pseudo",
           "log_fold_change"))) %>%
    mutate(shape = "sample") %>%
    mutate(`Comparator antigen` = gsub("WT", "D614G", `Comparator antigen`)) %>%
    mutate(log_fold_change = -log_fold_change) -> forest_data_sub
    
  #Add mean in 2 week time slots to plot
  dates <- seq(from = min(forest_data_sub$standardise_time, na.rm = T), to = max(forest_data_sub$standardise_time, na.rm = T), by = 14)
  
  # cumulative means per day
  dates <- sort(as.numeric(unique(forest_data_sub$time_since_first)))
  
  forest_data_sub <- add_time_group(forest_data_sub, dates, since_first = TRUE)
  
  all_sr_groups <- levels(forest_data_sub$standardise_encounters)
  
# forest_data_plot <- calculate_mean_over_time(forest_data_sub, dates, since_first = TRUE)
 
 cum_means <- calculate_cumulative_mean(forest_data_sub, dates)
 

  for(comp_ag in c("D614G")) {
    
    
    target_sr_groups <- c("WT conv","2x Vax", "3x Vax", "Vax + BA.1")
    
    lowess <- plot_over_time_lowess(forest_data_sub %>% filter(time_since_first < 44) , comp_ag, target_sr_groups, ymax = 0.5, ymin = -7,
                                    span = 0.75, cum_means = cum_means %>% filter(time_since_first < 44),
                                    omicron_variant_levels = c("BA.1"), show_loess = FALSE)
    
    lowess_long <- plot_over_time_lowess(forest_data_sub, comp_ag, target_sr_groups, ymax = 0.5, ymin = -7,
                                         span = 0.75, cum_means = cum_means,
                                         omicron_variant_levels = c("BA.1"), show_loess = FALSE)
    
    ggsave(file.path(path_to_save, paste0(comp_ag, "_sub_cumMean_44d_fc.", fileext)), lowess, dpi = 300, width = 7, height = 3)
    ggsave(file.path(path_to_save, paste0(comp_ag, "_sub_cumMean_long_fc.", fileext)), lowess_long, dpi = 300, width = 7, height = 3)
    
    lowess <- plot_over_time_lowess(forest_data_sub %>% filter(time_since_first < 44) , comp_ag, target_sr_groups, ymax = 0.5, ymin = -7,
                                    span = 0.75, cum_means = cum_means %>% filter(time_since_first < 44),
                                    omicron_variant_levels = c("BA.1"), show_loess = FALSE, cum_means_long = cum_means)
    
    lowess_long <- plot_over_time_lowess(forest_data_sub, comp_ag, target_sr_groups, ymax = 0.5, ymin = -7,
                                         span = 0.75, cum_means = cum_means,
                                         omicron_variant_levels = c("BA.1"), show_loess = FALSE, cum_means_long = cum_means)
    
    ggsave(file.path(path_to_save, paste0(comp_ag, "_sub_cumMean_44d_toLong_fc.", fileext)), lowess, dpi = 300, width = 7, height = 3)
    ggsave(file.path(path_to_save, paste0(comp_ag, "_sub_cumMean_long_toLong_fc.", fileext)), lowess_long, dpi = 300, width = 7, height = 3)
    

   
  }
    

}
