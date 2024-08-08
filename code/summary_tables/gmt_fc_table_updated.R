# ------------------------------------------------------------------------------
#-------------------------------- IMPORTANT ------------------------------------
# ------------------------------------------------------------------------------
# Restart R before running this script as other packages and Rmisc might interfere
# and result in inaccurate grouping 
# ------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# Setup workspace
rm(list = ls())
library(Rmisc)
library(tidyverse)
library(ggplot2)
library(meantiter)
library(huxtable)
library(patchwork)
library(rstatix)



# --------------------------------------------- load required functions -----------------------------------------------
working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
# figures_dir <- file.path(working_dir, "figures", "since_exposure")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
#----------------------------------------------- set path to save -----------------------------------------------
fileext <- "png"
plot_titers <- TRUE

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source(file.path("code", "summary_tables", "table_utility_functions.R"))

#-------------------------------------  SET TABLE NAME
# table_save
# table_no_save
# table_no_WD
# table_WD_save, etc...

# here for multiple tables
table_names = list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')

#for (table_name in table_names){
for (table_name in table_names){
  
  
  # create path to save for each table
  tab_name <- strsplit(table_name, "\\.")[[1]][1]
  
  #----------------------------------------------- set path to save -----------------------------------------------
  
  path_to_save <- file.path("data", "summary_tables", tab_name)
  dir.create(path_to_save, recursive = TRUE)
  
  
  table_path <- file.path(google_sheets_dir,table_name)
  
  
  #----------------- load and prepare data
  forest_data <- read.csv(table_path)
  
  
  forest_data <- format_table_forest_plots(forest_data) %>%
    filter(Webplotdigitizer != "y") 
  
  
  # ----------------------------------------------- prepare data -----------------------------------------------
  # rename all WT-like antigens 
  forest_data$`Comparator antigen`[forest_data$`Comparator antigen` %in% c("Wu-1", "D614G","B.1", "Wu-1?", "WT", "WA1")] <- "WT"
  
  #-------------------- Prepare data frame for mean fold drops
  
  mean_drop_df <- forest_data
  
  # Add fold drop to WT for other antigens than Omicron
  mean_drop_df %>%
    filter(OmicronVariant == "BA.1") %>%
    select(Study, Sourcelink, `Sera details long`,`Comparator antigen`, time, TitersHAg, `Assay Type`, standardised_assay, standardise_encounters) -> fold_drop_wt
  
  na_studies <- fold_drop_wt[is.na(fold_drop_wt$TitersHAg), "Study"]
  
  fold_drop_wt %>% filter(!(Study %in% na_studies)) %>%
    unique() %>%
    pivot_wider(names_from = `Comparator antigen`, values_from = TitersHAg) -> wide_fold_drop
  
  fold_drop_wt %>% filter(!(Study %in% na_studies)) %>%
    unique() %>%
    dplyr::group_by(Study, Sourcelink, `Sera details long`,`Assay Type`, standardised_assay, time, `Comparator antigen`, standardise_encounters) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n > 1L)
  
  wide_fold_drop %>%
    pivot_longer(cols = c("Alpha", "Beta","Gamma", "Delta")) %>%
    mutate(drop_to_wt = as.numeric(WT)/as.numeric(value)) %>%
    filter(!is.na(drop_to_wt))-> wide_fold_drop
  
  # calculate significant difference between PV and LV fold drop to wt
  wide_fold_drop %>%
    mutate(log2_fc = log2(drop_to_wt)) %>%
    select(name, log2_fc, standardised_assay, standardise_encounters) -> stat_wt_fc
  
  
  # match fold drop to WT for other variant antigens
  mean_drop_df$fold_drop_to_WT <- log2(wide_fold_drop$drop_to_wt[match(interaction(mean_drop_df[,c("Study", "Sera details long","Comparator antigen")]), 
                                                                       interaction(wide_fold_drop[,c("Study", "Sera details long","name")]))])
  
  # match the rows to add Omicron, WT comparison in wide format
  mean_data_wt <- mean_drop_df %>% filter(`Comparator antigen` == "WT")
  mean_data_wt[,c("fold_drop_to_WT", "TitersHAg", "Log2HAg")] <- mean_data_wt[,c("log_fold_change", "TitersOmicron", "Log2Omi")]
  mean_data_wt$`Comparator antigen` <- mean_data_wt$OmicronVariant
  
  # bind the wide variant to WT fold drop 
  mean_drop_df <- rbind(mean_drop_df, mean_data_wt)
  

  #------------------------------ Calculate mean tables ------------------------------
  # -------------- Means by serum group --------------------
  # not the same for WT when doing mean of fold drops and GMT WT/ GMT Omicron because some studies
  # have only fold drops and no titers (e.g. Suzuki)
  saveRDS(mean_drop_df, "code/summary_tables/mean_drop_df.rds")
  
  mean_drop_df %>% 
    # calculate now the ones for the OmicronVariants
    group_by(`Comparator antigen`,standardise_encounters, OmicronVariant) %>%
    mutate(mean_fold_drop_m = round(2^CI(log_fold_change)[["mean"]],1), 
           mean_fold_drop_lower = round(2^CI(log_fold_change)[["lower"]],1), 
           mean_fold_drop_upper = round(2^CI(log_fold_change)[["upper"]],1), 
           gmt_o = round(2^CI(na.omit(Log2Omi[Webplotdigitizer == "n"]))[["mean"]]*10),
           gmt_lower_o = round(2^CI(na.omit(Log2Omi[Webplotdigitizer == "n"]))[["lower"]]*10),
           gmt_upper_o = round(2^CI(na.omit(Log2Omi[Webplotdigitizer == "n"]))[["upper"]]*10),
           mean_fold_drop = paste0(mean_fold_drop_m, " (", mean_fold_drop_lower, "; ",mean_fold_drop_upper,")\n n=", length(log_fold_change[!is.na(log_fold_change)])),
           gmt_Omic = paste0(gmt_o," (", gmt_lower_o, "; ", gmt_upper_o, ")\n n=", length(Log2Omi[!is.na(Log2Omi) & Webplotdigitizer == "n"]))) %>%
    ungroup() %>%
    # calculate here the ones for the Comparator antigens
    group_by(`Comparator antigen`,standardise_encounters) %>%
    mutate(mean_fold_drop_to_WT_m = round(2^CI(na.omit(fold_drop_to_WT[OmicronVariant == "BA.1"]))[["mean"]],1), 
           mean_fold_drop_to_WT_lower = round(2^CI(na.omit(fold_drop_to_WT[OmicronVariant == "BA.1"]))[["lower"]],1),
           mean_fold_drop_to_WT_upper = round(2^CI(na.omit(fold_drop_to_WT[OmicronVariant == "BA.1"]))[["upper"]],1),
           gmt = round(2^CI(na.omit(Log2HAg[OmicronVariant == "BA.1" & Webplotdigitizer == "n"]))[["mean"]]*10),
           gmt_lower = round(2^CI(na.omit(Log2HAg[OmicronVariant == "BA.1" & Webplotdigitizer == "n"]))[["lower"]]*10),
           gmt_upper = round(2^CI(na.omit(Log2HAg[OmicronVariant == "BA.1" & Webplotdigitizer == "n"]))[["upper"]]*10),
           mean_fold_drop_to_WT = paste0(mean_fold_drop_to_WT_m, " (", mean_fold_drop_to_WT_lower, "; ",mean_fold_drop_to_WT_upper,")\n n=", length(fold_drop_to_WT[!is.na(fold_drop_to_WT) & OmicronVariant == "BA.1"])),
           gmt_hAG = paste0(gmt, " (", gmt_lower, "; ", gmt_upper, ")\n n=", length(Log2HAg[!is.na(Log2HAg) & OmicronVariant == "BA.1" & Webplotdigitizer == "n"]))) %>%
    select(standardise_encounters, OmicronVariant, `Comparator antigen`, mean_fold_drop, mean_fold_drop_to_WT,
           gmt_hAG, gmt_Omic) %>% 
    unique()-> mean_data_encounter
  
  
  mean_data_encounter <- mean_data_encounter[order(mean_data_encounter$`Comparator antigen`),]
  
  target_vars <- "BA.1"
  my_doc <- officer::read_docx()
  for(i in 1:length(target_vars)){
    
    my_doc <- mean_ecounter_table_per_ag(mean_data_encounter, target_vars[i], my_doc, i, only_omicron = i>1)  
  }
  
  print(my_doc, target =
          paste0(path_to_save,paste0("/mean_encounters.docx")))
  
  
}
