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
  forest_data_b <- read.csv(table_path)
  
  
  forest_data_b <- format_table_forest_plots(forest_data_b) %>%
    filter(Webplotdigitizer != "y")

  forest_data <- standardise_time_format(forest_data_b)
  
  forest_data$standardise_time[forest_data$Study == "HKU"] <- as.Date("2021/12/12")
  forest_data$standardise_time[forest_data$Study == "Krumbholz"] <- as.Date("2022/03/29")
  
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
  target_sr_groups <- c("WT conv","2x Vax", "3x Vax", "Vax + BA.1")
  
 cum_means <- calculate_cumulative_mean(forest_data_sub %>%
                                          filter(`Comparator antigen` == "D614G") %>%
                                          filter(standardise_encounters %in% target_sr_groups) %>%
                                          filter(OmicronVariant == "BA.1"), dates)
 
 final_mean <- cum_means %>%
   filter(date_group == max(cum_means$date_group)) %>%
   mutate(Study = "Average")
 
 # filter to select studies
 forest_data_b %>%
   mutate(`Comparator antigen` = gsub("WT", "D614G", `Comparator antigen`)) %>%
   mutate(log_fold_change = -log_fold_change) %>%
   filter(Study %in% c("Kimpel/Roessler", "Montefiori/Doria-Rose")) %>%
   filter(standardise_encounters %in% target_sr_groups) %>%
   filter(`Comparator antigen` == "D614G") %>%
   filter(OmicronVariant == "BA.1") %>%
   select(Study, `Sera details long`, standardise_encounters, `Comparator antigen`, OmicronVariant, log_fold_change) %>%
   mutate(shape = "sample") -> target_studies
   
target_studies %>%
  ggplot(aes(x = standardise_encounters, y = log_fold_change, color = Study, group = Study)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey20") +
  geom_pointrange(data = final_mean %>%
                    filter(standardise_encounters != "Vax + BA.1"), aes(ymin = lower, ymax = upper), alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) + 
  scale_color_discrete(name = "Data") +
  xlab("Serum group") + 
  ylab("Log2 titer fold change from 614D/G to BA.1") +
  ylim(c(-7,0)) +
  theme_bw() -> p

  ggsave(file.path(path_to_save, "mean_final_fold_drop_to_single_study.png"), p, dpi = 300, width = 4, height = 5)


}
