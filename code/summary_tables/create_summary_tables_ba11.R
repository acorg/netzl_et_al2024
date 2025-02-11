# create forest plots for the different tables and variants
# setup page and load metadata
library(meantiter)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(huxtable)
library(ggplot2)

rm(list = ls())

# set target variants, target comparator antigen, grouping variable, and serum groups
target_variants <-  c("BA.1", "BA.1.1")
target_sr_groups <- c("2x Vax", "3x Vax", "Inf + Vax", "Vax + Inf",
                      "Vax + BA.1", "Vax + BA.2", "Delta conv")
comp_antigen <- "WT"
grouping_var <- "standardise_encounters"

# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
table_dir <- file.path(data_dir, "summary_tables")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
#----------------------------------------------- set path to save -----------------------------------------------

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path(working_dir,"code", "summary_tables", "table_utility_functions.R"))


#-------------------------------------  SET TABLE NAME
# here for multiple tables
table_names = list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')


for (table_name in table_names){
  
  # create path to save for each table
  tab_name <- strsplit(table_name, "\\.")[[1]][1]
  
  path_to_save <- file.path(table_dir, tab_name)
  
  figures_dir <- file.path(working_dir, "figures", "forest_plots", tab_name, "titers")
  
  # create subfolder for boxplots
  suppressWarnings(dir.create(path_to_save, recursive = TRUE))
  
  
  table_path <- file.path(google_sheets_dir,table_name)
  
  summary_tab_name <- paste0(comp_antigen,"-", paste0(target_variants, collapse = "_"),"-", paste0(grouping_var, collapse = "_"), ".docx") %>%
    gsub("\\/", "\\.",.)
  
  #----------------- load and prepare data
  forest_data <- read.csv(table_path)
  
  if(grepl("lessWD", table_name)){
    forest_data <- forest_data %>%
      filter(Webplotdigitizer != "y")
  }
  
  forest_data <- format_table_forest_plots(forest_data)
  
  
  forest_data <- melt_omicron_data(forest_data)
  
  # match sera name long
  forest_data$`Sera details long` <- gsub("woR346K|wR346K|wo R346K|w R346K", "", forest_data$`Sera details long`)
  
  forest_ba1 <- forest_data %>%
    mutate(`Sera details long` = paste(Study, `Sera details long`))
  
  ba1_ba11_samples <- Reduce(intersect, list(forest_ba1 %>% filter(OmicronVariant == "BA.1") %>% pull(`Sera details long`),
                                             forest_ba1 %>% filter(OmicronVariant == "BA.1.1") %>% pull(`Sera details long`)))
  
  forest_ba1 %>%
    filter(`Sera details long` %in% ba1_ba11_samples) -> forest_ba1
  stop()
  # ----------------------------- do the same for titers
  
  if(grepl("lessWD", table_name)){
    forest_data <- forest_data %>%
      filter(Webplotdigitizer == "n")
    
    forest_ba1 <- forest_ba1 %>%
      filter(Webplotdigitizer == "n")
  }
  
  full_titer <- wider_titer_variant_tables(forest_ba1, variants = target_variants, group_var = grouping_var,
                                           comp_antigen = comp_antigen, match_sera = FALSE)
  
  full_titer <- full_titer %>%
    filter(Lineage == "BA.1")
  formatted_titer <- format_huxtable_variant_table(full_titer, target_variants, column_names = c("Comparator antigen",
                                                                                                 "Serum group"))
  
  save_summary_table(formatted_titer, paste0("mean_titer-", summary_tab_name), table_dir, tab_name)
  
  forest_ba1 %>%
    filter(`Comparator antigen` == comp_antigen) %>%
    select(standardise_encounters, OmicronVariant) %>%
    table
  # make plot
  forest_ba1 %>%
    filter(`Comparator antigen` == comp_antigen) %>%
    filter(OmicronVariant %in% target_variants) %>%
    ggplot(aes(x = OmicronVariant, y = Log2Omi, group = `Sera details long`)) + 
    geom_line() +
    geom_point() + 
    xlab("Omicron lineage") +
    scale_y_continuous(labels = function(x) round(2^x*10),
                       breaks = -1:8,
                       name = "Neutralization titer") +
    facet_wrap(~ standardise_encounters, ncol = length(target_sr_groups)) +
    theme_bw() + 
    theme(strip.background.x = element_blank()) -> p
  

    ggsave(file.path(figures_dir, paste0("paired_mean_titer_sub-", gsub("docx", "png", summary_tab_name))), p, width = 12, height = 4, units = "in", dpi = 300, ) 
}
