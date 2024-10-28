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


# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "forest_plots")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
#----------------------------------------------- set path to save -----------------------------------------------
fileext <- "png"

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path(working_dir, "code", "plotting", "sublineage_forest_plots.R"))
source(file.path(utility_dir, "plot_over_time_functions.R"))

#-------------------------------------  SET TABLE NAME
table_names <- list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')


sublineages_of_interest <- c("BA.1")

do_month_one <- TRUE

for (table_name in table_names){
  
  
  # create path to save for each table
  tab_name <- strsplit(table_name, "\\.")[[1]][1]
 
  path_to_save <- file.path(figures_dir, tab_name)
  
  # create subfolder for fold changes and titer drops
  suppressWarnings(dir.create(file.path(path_to_save, "fold_drops"), recursive = TRUE))
  suppressWarnings(dir.create(file.path(path_to_save, "titers"), recursive = TRUE))
  
  table_path <- file.path(google_sheets_dir,table_name)
  
  
#----------------- load and prepare data
  forest_data <- read.csv(table_path)
  
  if(grepl("lessWD", table_name)){
    forest_data <- forest_data %>%
      filter(Webplotdigitizer != "y")
  }
  
  forest_data$SAVE_lab <- "n"

  forest_data <- format_table_forest_plots(forest_data)
  
  for(sl in sublineages_of_interest) {
    
    if(do_month_one & sl == "BA.1"){
      
      forest_data <- standardise_time_format(forest_data)
      
      forest_data$standardise_time[forest_data$Study == "HKU"] <- as.Date("2021/12/12")
      forest_data$standardise_time[forest_data$Study == "Krumbholz"] <- as.Date("2022/03/29")
      
      
      forest_data %>%
        mutate(year = as.Date(gsub("\\.", "\\/", year), format = "%Y/%m/%d")) %>%
        filter(year < "2021/12/15") %>%
        filter(`Comparator antigen` %in% c("WT", "D614G")) %>%
        filter(standardise_encounters %in% c("2x Vax", "3x Vax", "WT conv")) %>%
        mutate(standardise_encounters = factor(standardise_encounters, levels = rev(c("WT conv", "2x Vax", "3x Vax")))) -> forest_data_sub
      
      path_to_save_sub <- file.path(path_to_save, "2_weeks")
      
      suppressWarnings(dir.create(file.path(path_to_save_sub, "fold_drops"), recursive = TRUE))
      suppressWarnings(dir.create(file.path(path_to_save_sub, "titers"), recursive = TRUE))
      
      forest_data_sub <- forest_data_sub %>%
        arrange(standardise_encounters) %>%
        mutate(Row_long = 1:nrow(forest_data_sub))
      
      plot_forest_plot_by_sublineage(forest_data_sub, sublineage = c(sl), group_by = c("standardise_encounters"), 
                                     path_save = path_to_save_sub, axis_text_size = 8, max_plot_height = 6)
      
    }
    
    if(sl == "BA.1"){
      plot_forest_plot_by_sublineage(forest_data, sublineage = c(sl), group_by = c("standardise_encounters"), 
                                     path_save = path_to_save, axis_text_size = 4, max_plot_height = 16)
      
    } else {
      plot_forest_plot_by_sublineage(forest_data, sublineage = c(sl), group_by = c("standardise_encounters"), 
                                     path_save = path_to_save)
    }
   
  }
  
  
}
