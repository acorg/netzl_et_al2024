# plot study by country of corresponding author to show data distribution
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

set.seed(100)

# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "over_time")
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
    filter(Webplotdigitizer != "y") 
  
  # select only 2x Vax, 3x Vax and BA.2 Omicron
  
  forest_data %>% 
    filter(standardise_encounters %in% c("2x Vax", "3x Vax")) %>%
    filter(`Comparator antigen` %in% c("D614G", "WT")) %>%
    filter(OmicronVariant == "BA.1") -> forest_data
  
  # select randomly subset of data of size n with replacement
  # do this 100 times per n to calculate confidence interval
  data_sub <- forest_data %>%
    select(Study, `Sera details long`, OmicronVariant, standardise_encounters, log_fold_change)
  
  n_studies <- data_sub %>%
    select(standardise_encounters) %>%
    table()
  
  calc_mean_ci_for_random_subset <- function(data, n){
    
    sub_rows <- sample(length(data), n, replace = TRUE)
    
    summary_stats <- Rmisc::CI(data[sub_rows])
    
    df <- data.frame(n = n,
                     mean = summary_stats["mean"],
                     lower = summary_stats["lower"],
                     upper = summary_stats["upper"])
    
    return(df)
  }
  
  random_subset_n_times <- function(data, n_times, n_resamples){
    
    res <- lapply(1:n_times, function(x){
      calc_mean_ci_for_random_subset(data, n_resamples) %>%
        mutate(n_time = x)
    })
    
    res <- do.call(rbind, res)
    
    return(res)
  }
  
  bootstrap_list <- list()
  for(srg in c("2x Vax", "3x Vax")){
    
    vax_data <- data_sub %>%
      filter(standardise_encounters == srg) %>%
      pull(log_fold_change)
    
    temp_res <- lapply(1:length(vax_data), function(x){
      random_subset_n_times(vax_data, 100, x) %>%
        mutate(serum_group = srg,
               n_times = 100)
    })
    
    temp_res_sub <- lapply(1:length(vax_data), function(x){
      random_subset_n_times(vax_data, round(0.1*length(vax_data)), x) %>%
        mutate(serum_group = srg,
               n_times = round(0.1*length(vax_data)))
    })
    
    bootstrap_list[[srg]] <- rbind(do.call(rbind, temp_res),
                                   do.call(rbind, temp_res_sub))
    
  }
  
  
  bootstrap_data <- do.call(rbind, bootstrap_list)
  
  bootstrap_summary <- bootstrap_data %>%
    group_by(n, serum_group, n_times) %>%
    summarize(lower = Rmisc::CI(mean)["lower"],
              upper = Rmisc::CI(mean)["upper"],
              mean = Rmisc::CI(mean)["mean"],
              CI_width = upper - lower)
  
  # plot it 
  bootstrap_summary %>%
    filter(n_times < 100) %>%
    ggplot(aes(x = n, y = -mean)) + 
    geom_ribbon(aes(ymin = -lower, ymax = -upper), alpha = 0.2, color = NA) + 
    geom_line() + 
    facet_grid(~serum_group) + 
    scale_y_continuous(labels = function(x) paste0("-",round(2^abs(x), 2), "x"),
                       name = "Fold change from WT") +
    scale_x_continuous(name = "Number of subsamples", 
                       breaks = seq(0, 110, 10)) +
    theme_bw() + 
    theme(strip.background = element_blank())-> mean_bootstrap
  
  bootstrap_summary %>%
    filter(n_times < 100) %>%
    ggplot(aes(x = n, y = CI_width)) + 
    geom_line() + 
    facet_grid(~serum_group) + 
    scale_y_continuous(labels = function(x) paste0(round(2^x, 2), "x"),
                       name = "95%CI of the mean width") +
    scale_x_continuous(name = "Number of subsamples", 
                       breaks = seq(0, 110, 10)) +
    theme_bw() +
    theme(strip.background = element_blank()) -> ci_width
  
  mean_bootstrap / ci_width + plot_annotation(tag_levels = "A") -> p
  
  ggsave(file.path(figures_dir, gsub(".csv", "", table_name), "bootrstrap_n_studies.png"), p, dpi = 300, width = 8, height = 6)
  
}
