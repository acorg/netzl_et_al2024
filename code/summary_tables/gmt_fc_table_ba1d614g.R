# create forest plots for the different tables and variants
# setup page and load metadata
rm(list = ls())
library(Rmisc)
library(tidyverse)
library(ggplot2)
library(meantiter)
library(huxtable)
library(patchwork)
library(rstatix)

# pkgs <- c("huxtable", "effectsize", "flextable", "broom", "report")
# install_if_not_installed(pkgs)


# to suppress NA removal warnings from ggplot
options(warn=-1)

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
source("code/summary_tables/table_utility_functions.R")

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
 
  path_to_save <- file.path("data", "summary_tables", tab_name)
  dir.create(path_to_save, recursive = TRUE)
  
  # create subfolder for boxplots
#  suppressWarnings(dir.create(path_to_save, recursive = TRUE))
  
  
  table_path <- file.path(google_sheets_dir,table_name)
  
  
#----------------- load and prepare data
  forest_data <- read.csv(table_path)
  
  
  forest_data <- format_table_forest_plots(forest_data) %>%
    filter(Webplotdigitizer != "y") 
  
  # ------------------ assay and cell type impact on fold changes
  forest_data %>%
    select(all_of(c("Comparator antigen", "OmicronVariant", "Sera_details_no_time", "Log2HAg", "Log2Omi",
                    "Standardised_sera_names", "standardise_encounters", "vacc_type_het", "vacc_type","vaccine_manufacturer", "standardised_assay", "standardised_pseudo",
                    "standardised_cell", "log_fold_change", "Webplotdigitizer"))) %>%
    mutate(`Comparator antigen` = gsub("WT", "D614G", `Comparator antigen`)) -> forest_data_sub_b
  
  
  # Enough data only for WT vs BA.1 in these sr groups
  target_sr_groups <- c("Inf + Vax", "3x Vax", "2x Vax", "WT conv", "Vax + BA.1")
  
  # Most data for BA.1 and D614G -> do these as comparison for fold change
  forest_data_sub_b %>%
    filter(`Comparator antigen` == "D614G") %>%
    filter(OmicronVariant == "BA.1") %>%
    filter(standardised_assay %in% c("LV", "PV")) -> forest_data_wt_omi
  
  forest_data_wt_omi %>%
    group_by(standardised_assay, standardise_encounters) %>%
    mutate(count_n = n()) %>%
    ungroup() -> forest_data_sub
  
  # Fold change from D614G for other variants in gmt_fc_table_updated.R
  
  # check for normality
  forest_data_sub %>%
    filter(count_n > 2) %>%
    group_by(standardise_encounters, standardised_assay) %>%
    mutate(shapiro_p = shapiro_test(log_fold_change)$p.value) %>%
    group_by(standardise_encounters) %>%
    mutate(shapiro_comb = min(shapiro_p),
           `Fold change from WT` = log_fold_change,
           both_assays = ("PV" %in% standardised_assay) & ("LV" %in% standardised_assay)) -> forest_data_sub
  
  forest_data_sub %>%
    filter(both_assays) %>%
    filter(shapiro_comb > 0.05) %>%
    filter(count_n > 1) %>%
    group_by(standardise_encounters) %>%
    t_test(log_fold_change~standardised_assay) %>%
    add_significance() %>%
    select("Serum group" = standardise_encounters,
           "Variable" = .y.,
           "Group 1" = group1,
           "Group 2" = group2,
           n1, n2, statistic, df, p,
           "Significance" = p.signif) %>%
    mutate(Variable = gsub("log_fold_change", "Fold change from WT", Variable)) %>%
    arrange(p)-> fc_from_wt_ttest
  
  # Show below
  forest_data_sub %>%
    filter(both_assays) %>%
    filter(shapiro_comb <= 0.05) %>%
    filter(count_n > 1) %>%
    group_by(standardise_encounters) %>%
    wilcox_test(log_fold_change~standardised_assay) %>%
    add_significance() %>%
    select("Serum group" = standardise_encounters,
           "Variable" = .y.,
           "Group 1" = group1,
           "Group 2" = group2,
           n1, n2, statistic, p,
           "Significance" = p.signif) %>%
    mutate(Variable = gsub("log_fold_change", "Fold change from WT", Variable)) %>%
    arrange(p) -> fc_from_wt_wilcox
  
  fc_ttest <- format_stat_table_flex(fc_from_wt_ttest)
  fc_wilcox <- format_stat_table_flex(fc_from_wt_wilcox)
  
  
  stat_doc <- officer::read_docx() #"data/summary_tables/omicron_folddrops_preprocessed_wSAVE_lessWD/stat_test_assay_FCfromWT.docx")
  
  stat_doc %>%
    officer::body_add_par("Supplementary Table XX: Assay based statistical comparison of fold changes from WT/D614G to BA.1. A t-test was performed to compare pseudovirus (PV) and live virus (LV) assessed fold changes in different serum groups. Normality was checked with a Shapiro-Wilk test.") %>%
    flextable::body_add_flextable(fc_ttest) %>%
    officer::body_add_par("Supplementary Table XX: Assay based statistical comparison of fold changes from WT/D614G to BA.1. A Wilcoxon-test was performed to compare pseudovirus (PV) and live virus (LV) assessed fold changes in different serum groups. Normality was checked with a Shapiro-Wilk test.") %>%
    flextable::body_add_flextable(fc_wilcox) -> stat_doc
  

  print(stat_doc, target =
          paste0(path_to_save,paste0("/stat_test_assay_BA1FCfromWT.docx")))
 
 }

