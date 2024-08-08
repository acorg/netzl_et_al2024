#setup page and load metadata
rm(list = ls())

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
google_sheets_dir <- file.path(data_dir, 'google_sheet_tables')
titer_tables_dir <- file.path(data_dir, 'titer_tables')

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path("code","summary_tables", 'table_utility_functions.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path(utility_dir, "plot_over_time_functions.R"))

#-------------------------------------  SET TABLE NAME

# here for multiple tables
table_names <- list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')

for(early_only in c("only30d","only60d", "only90d", "")){
  
  for (table_name in table_names){
    
    table_path <- file.path(google_sheets_dir,table_name)
   
    forest_data <- read.csv(table_path)
    
    forest_data <- format_table_forest_plots(forest_data, remove_4xVax_ba45bt_samples = FALSE) %>%
      filter(Webplotdigitizer == "n")
    
    
    if(early_only == "only30d"){
      
      forest_data <- standardise_time_format(forest_data, group_by_sr_group = FALSE, group_not_by_comparator = TRUE ,group_by_omicron = FALSE)
      
      forest_data <- forest_data %>%
        filter(time_since_first <= 30)
    } else if(early_only == "only60d"){
      
      forest_data <- standardise_time_format(forest_data, group_by_sr_group = FALSE, group_not_by_comparator = TRUE ,group_by_omicron = FALSE)
      
      forest_data <- forest_data %>%
        filter(time_since_first <= 60)
    } else if(early_only == "only90d"){
      
      forest_data <- standardise_time_format(forest_data, group_by_sr_group = FALSE, group_not_by_comparator = TRUE ,group_by_omicron = FALSE)
      
      forest_data <- forest_data %>%
        filter(time_since_first <= 90)
    }
    
    
    # get omicron titre
    forest_data <- melt_omicron_data(forest_data)
    
    forest_data <- forest_data %>% filter(!is.na(TitersHAg))
    forest_data$Study <- gsub("Khan2", "Khan", forest_data$Study) # Make them the same study so that it is one serum that contains BA.1 and BA.4/5
    forest_data$serum_name <- paste(forest_data$standardise_encounters, forest_data$standardised_assay, forest_data$standardised_pseudo, forest_data$vaccine_manufacturer,forest_data$Standardised_sera_names,forest_data$time, forest_data$Sera_details_no_time, forest_data$Study, sep = "_")
    

    # ------------------------------------------------------- Make titer table 
    # take mean of titers if more than one omicron titer present, eg against both WT and D614G
    mean_titer <- function(titers) {
      if(length(titers) == 0) {
        "*"
      } else if(length(titers) >1) {
        mean_t <- mean(log2(titers/10))
        2^mean_t*10
      } else {
        titers
      }
    }
    
    # create table
    forest_data %>% 
      select(`Comparator antigen`, serum_name, TitersHAg) %>%
      pivot_wider(names_from = serum_name, values_from = TitersHAg, values_fn = mean_titer) %>%
      column_to_rownames("Comparator antigen")-> omicron_table
    
    omicron_table <- as.data.frame(omicron_table)
    omicron_table[is.na(omicron_table)] <- "*"
    
    
    # set rownames for antigens
    ag_names <- c("Delta" = "B.1.617.2", "Alpha" = "B.1.1.7", "Beta" = "B.1.351", "Gamma" = "P.1")
    
    ag_names_map <- function(x){
      if (x %in% names(ag_names)){ag_names[x]}
      else {x}
    }
    
    rownames(omicron_table) <- unlist(lapply(rownames(omicron_table), ag_names_map), 
                                      use.names=FALSE)
    
    omicron_table[,grepl("Sigal/Khan", colnames(omicron_table))]
    
  
    # save titer table
    table_name <- str_replace(table_name, 'folddrops', paste0('titer_',early_only))
    
    write.csv(x = omicron_table, file = file.path(titer_tables_dir, table_name))
  }
  
}
