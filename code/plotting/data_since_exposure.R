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
figures_dir <- file.path(working_dir, "figures", "since_exposure")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
#----------------------------------------------- set path to save -----------------------------------------------
fileext <- "png"
plot_titers <- TRUE

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path(utility_dir, "plot_since_exposure_functions.R"))

#-------------------------------------  SET TABLE NAME

# here for multiple tables
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
    filter(Webplotdigitizer != "y") %>%
    standardise_time_to_mean_days(.)
 
  forest_data %>%
    select(all_of(c("Comparator antigen", "OmicronVariant", "Sera_details_no_time", "Log2HAg", "Log2Omi",
                    "Standardised_sera_names", "standardise_encounters", "vacc_type_het", "vacc_type","vaccine_manufacturer", "standardised_assay", "standardised_pseudo",
                    "log_fold_change", "mean_days"))) %>%
    mutate(shape = "sample") %>%
    mutate(`Comparator antigen` = gsub("WT", "D614G", `Comparator antigen`)) %>%
    mutate(log_fold_change = -log_fold_change) %>%
    filter(!is.na(mean_days)) -> forest_data_sub
  
  all_sr_groups <- levels(forest_data_sub$standardise_encounters)
  
  target_sr_groups <- c("WT conv", "2x Vax", "3x Vax", "Vax + BA.1", "Inf + Vax")
  
  fc_exposure <- fold_change_since_exposure(forest_data_sub %>%
                                              mutate(OmicronVariant = paste0("FC to ", OmicronVariant)), variants = c("FC to BA.1"), target_sr_groups = target_sr_groups, comp_antigen = "D614G",
                             ymax = 0.5, ymin = -7)
  
  ggsave(file.path(path_to_save, paste0("fc_since_exposure.", fileext)), fc_exposure, dpi = 300, width = 8, height = 3.5)
  
  # now titers
  forest_data %>%
    filter(Webplotdigitizer == "n") %>%
    select(all_of(c("Comparator antigen","Study", "OmicronVariant", "Sera_details_no_time", "Log2HAg", "Log2Omi",
                    "Standardised_sera_names", "standardise_encounters", "vacc_type_het", "vacc_type","vaccine_manufacturer", "standardised_assay", "standardised_pseudo",
                    "log_fold_change", "mean_days", "Sera details long"))) %>%
    mutate(shape = "sample") %>%
    mutate(`Comparator antigen` = gsub("WT", "D614G", `Comparator antigen`)) %>%
    mutate(log_fold_change = -log_fold_change) %>%
    filter(!is.na(mean_days)) -> forest_data_sub
  
  
  forest_data_titer <- forest_data_sub %>%
    mutate(Log2Omi = Log2HAg,
           #    `Comparator antigen` = "D614G", 
           OmicronVariant = `Comparator antigen`) %>%
    filter(OmicronVariant %in% c("D614G", "Beta", "Delta", "Alpha")) %>%
    rbind(., forest_data_sub) %>%
    mutate(log_fold_change = Log2Omi) %>%
    unique()
  
  write.csv(forest_data_titer, file.path(data_dir, "titer_tables", paste0(tab_name, "_titers_since_exposure.csv")), row.names = FALSE)
  
  forest_data_titer <- forest_data_titer %>%
    select(!`Sera details long`)
  
  titers <- titer_since_exposure(forest_data_titer %>%
                                   mutate(OmicronVariant = gsub("D614G", "614D/G", OmicronVariant)) %>%
                                   mutate(OmicronVariant = paste(OmicronVariant, "GMT")), ymin = -1, ymax = 13, target_sr_groups = c(target_sr_groups), variants = c("614D/G GMT", "BA.1 GMT"))
  
  ggsave(file.path(path_to_save, paste0("titers_since_exposure.", fileext)), titers, dpi = 300, width = 8, height = 6)

}
