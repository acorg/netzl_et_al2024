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
library(rstatix)
library(rempsyc)
# library(huxtable)
# library(flextable)
library(knitr)
library(kableExtra)

pkgs <- c("effectsize", "flextable", "broom", "report", "coin")
install_if_not_installed(pkgs)


# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")
#----------------------------------------------- set path to save -----------------------------------------------

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))
source("code/summary_tables/table_utility_functions.R")

#-------------------------------------  SET TABLE NAME
# here for multiple tables
table_names = list('omicron_folddrops_preprocessed_wSAVE_lessWD.csv')

#for (table_name in table_names){
for (table_name in table_names){
  
  
  # create path to save for each table
  tab_name <- strsplit(table_name, "\\.")[[1]][1]
 
  path_to_save <- file.path("data", "summary_tables", tab_name)
  dir.create(path_to_save, recursive = TRUE)

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
  
  
  do_stat_test_and_effsize <- function(data_sub, test_formula="log_fold_change~standardised_assay", stat_test = "t"){
    
    if(stat_test == "t"){
      data_sub %>%
        filter(both_assays) %>%
        filter(shapiro_comb > 0.05) %>%
        filter(count_n > 1) %>%
        group_by(standardise_encounters) %>%
        t_test(as.formula(test_formula), detailed = TRUE) %>%
        add_significance() %>%
        mutate("95%CI (lower;upper)" = paste0(round(conf.low, 2), ";", round(conf.high, 2)),
               statistic = round(statistic, 2),
               p = round(p, 2),
               df = round(df, 1)) %>%
        select("Serum group" = standardise_encounters,
               "n(PV)" = n1,
               "n(LV)" = n2,
               statistic, `95%CI (lower;upper)`, df, p,
               "Significance" = p.signif) %>%
        arrange(p)-> fc_from_wt_ttest
      
      # we are not assuming same variance across the groups, as sample size and collection times vary. We are also dealing with small <50 sample sizes, 
      # so we do hedges correction
      data_sub %>%
        filter(both_assays) %>%
        filter(shapiro_comb > 0.05) %>%
        filter(count_n > 1) %>%
        group_by(standardise_encounters) %>%
        cohens_d(as.formula(test_formula), var.equal = FALSE, hedges.correction = TRUE) %>%
        mutate(effsize = round(effsize, 2)) %>%
        select("Serum group" = standardise_encounters,
               "Effect size" = effsize,
               "Magnitude" = magnitude) %>%
        left_join(fc_from_wt_ttest, ., by = "Serum group") -> tab_return
      
      return(tab_return)
    } else if(stat_test == "wilcoxon"){
      
      data_sub %>%
        filter(both_assays) %>%
        filter(shapiro_comb <= 0.05) %>%
        filter(count_n > 1) %>%
        group_by(standardise_encounters) %>%
        wilcox_test(as.formula(test_formula), detailed = TRUE) %>%
        add_significance() %>%
        mutate("95%CI (lower;upper)" = paste0(round(conf.low, 2), ";", round(conf.high, 2)),
               statistic = round(statistic, 2),
               p = round(p, 2)) %>%
        select("Serum group" = standardise_encounters,
               "n(PV)" = n1,
               "n(LV)" = n2,
               statistic,`95%CI (lower;upper)`, p,
               "Significance" = p.signif) %>%
        arrange(p) -> fc_from_wt_wilcox
      
      data_sub %>%
        filter(both_assays) %>%
        filter(shapiro_comb <= 0.05) %>%
        filter(count_n > 1) %>%
        group_by(standardise_encounters) %>%
        wilcox_effsize(as.formula(test_formula)) %>%
        mutate(effsize = round(effsize, 2)) %>%
        select("Serum group" = standardise_encounters,
               "Effect size" = effsize,
               "Magnitude" = magnitude) %>%
        left_join(fc_from_wt_wilcox, ., by = "Serum group") -> tab_return
      
      return(tab_return)
      
    } else {
      warning("Please specify 't' or 'wilcoxon' as stat test")
      return(NULL)
    }
    
  }
  
  fc_from_wt_ttest <- do_stat_test_and_effsize(forest_data_sub, test_formula = "log_fold_change~standardised_assay", stat_test = "t")
  fc_from_wt_wilcox <- do_stat_test_and_effsize(forest_data_sub, test_formula = "log_fold_change~standardised_assay", stat_test = "wilcoxon")
  
  kable(fc_from_wt_ttest, "latex") %>%
    row_spec(0,bold=TRUE)
  
  kable(fc_from_wt_wilcox, "latex") %>%
    row_spec(0,bold=TRUE)
  
  stat_doc <- officer::read_docx()
  # stat_doc <- officer::read_docx(path = paste0(path_to_save,paste0("/stat_test_assay_FCfromWT.docx")))
  
  stat_doc %>%
    officer::body_add_par("Supplementary Table XX: Assay based statistical comparison of fold changes from WT/D614G to BA.1. A t-test was performed to compare
pseudovirus (PV) and live virus (LV) assessed Geometric Mean Titers
(GMTs) in different serum groups. Normality was checked with a
Shapiro-Wilk test and these serum groups were found to not differ
significantly from a normal distribution. Results are ordered by
increasing p-value. A statistic above 0 indicates higher values in Group
1 than in Group 2. The effect size was calculated with Cohen's d for unequal variances and Hedge's correction due to small sample sizes. All tests were performed with the rstatix package.") %>%
    officer::body_add_table(.,fc_from_wt_ttest, style = "table_template") %>%
    officer::body_add_par("Supplementary Table XX: Assay based statistical comparison of fold changes from WT/D614G to BA.1. 
                           A
Wilcoxon-test was performed to compare pseudovirus (PV) and live virus
(LV) assessed fold changes in different serum groups. Normality was
checked with a Shapiro-Wilk test and these serum groups were found to
differ significantly from a normal distribution. Results are ordered by
increasing p-value. A statistic above 0 indicates higher values in Group
1 than in Group 2. The effect size was calculated with Wilcoxon effect size. All tests were performed with the rstatix package") %>%
    officer::body_add_table(.,fc_from_wt_wilcox, style = "table_template")-> stat_doc
  
  print(stat_doc, target =
          paste0(path_to_save,paste0("/stat_test_assay_FCfromWT.docx")))
  
  
  ## ------------ For TITERS melt it so that it is all in one column, and group by Comparator antigen 
  comp_df <- forest_data_sub_b %>%
    mutate(Log2HAg = Log2Omi,
           `Comparator antigen` = OmicronVariant)
  
  forest_data_comb <- rbind(comp_df, forest_data_sub_b) %>%
    filter(Webplotdigitizer == "n") %>%
    select(!OmicronVariant) %>%
    select(!Log2Omi) %>%
    select(!log_fold_change) %>%
    filter(standardised_assay %in% c("LV", "PV")) %>%
    filter(!(is.na(Log2HAg))) %>%
    unique()
  
  forest_data_comb %>%
    group_by(`Comparator antigen`, standardise_encounters, standardised_assay) %>%
    mutate(count_n = n()) %>%
    ungroup() %>%
    filter(`Comparator antigen` %in% c("D614G", "BA.1", "Alpha", "Beta", "Delta")) -> forest_data_comb
  
  forest_data_comb %>%
    filter(count_n > 2) %>%
    group_by(`Comparator antigen`, standardise_encounters, standardised_assay) %>%
    mutate(shapiro_p = shapiro_test(Log2HAg)$p.value) %>%
    group_by(`Comparator antigen`, standardise_encounters) %>%
    mutate(shapiro_comb = min(shapiro_p),
           `Titer` = Log2HAg,
           both_assays = ("PV" %in% standardised_assay) & ("LV" %in% standardised_assay)) %>%
    ungroup() -> forest_data_comb

  forest_data_comb %>%
    filter(both_assays) %>%
    filter(count_n > 1) %>%
    filter(shapiro_comb > 0.05) %>%
    group_by(`Comparator antigen`, standardise_encounters) %>%
    t_test(Titer~standardised_assay, detailed = TRUE) %>%
    add_significance() %>%
    mutate("95%CI (lower;upper)" = paste0(round(conf.low, 2), ";", round(conf.high, 2)),
           statistic = round(statistic, 2),
           p = round(p, 4),
           df = round(df, 1)) %>%
    select("Variant" = `Comparator antigen`,
          "Serum group" = standardise_encounters,
          "n(PV)" = n1,
          "n(LV)" = n2,
          statistic, `95%CI (lower;upper)`, df, p,
           "Significance" = p.signif) %>%
    mutate(Variant = gsub("D614G", "WT", Variant)) %>%
    arrange(p)-> titer_ttest
  
  
  # add effect size
  forest_data_comb %>%
    filter(both_assays) %>%
    filter(count_n > 1) %>%
    filter(shapiro_comb > 0.05) %>%
    group_by(`Comparator antigen`, standardise_encounters) %>%
    cohens_d(Titer~standardised_assay, var.equal = FALSE, hedges.correction = TRUE) %>%
    mutate(effsize = round(effsize, 2),
           `Comparator antigen` = gsub("D614G", "WT", `Comparator antigen`)) %>%
    select("Serum group" = standardise_encounters,
           "Effect size" = effsize,
           "Magnitude" = magnitude,
           "Variant" = `Comparator antigen`) %>%
    left_join(titer_ttest, ., by = c("Variant", "Serum group")) -> titer_ttest
  
  # wilcoxon test
  forest_data_comb %>%
    filter(both_assays) %>%
    filter(count_n > 1) %>%
    filter(shapiro_comb <= 0.05) %>%
    group_by(`Comparator antigen`, standardise_encounters) %>%
    wilcox_test(Titer~standardised_assay, detailed = TRUE) %>%
    add_significance() %>%
    mutate("95%CI (lower;upper)" = paste0(round(conf.low, 2), ";", round(conf.high, 2)),
           statistic = round(statistic, 2),
           p = round(p, 4)) %>%
    select("Variant" = `Comparator antigen`,
           "Serum group" = standardise_encounters,
           "n(PV)" = n1,
           "n(LV)" = n2,
           statistic,`95%CI (lower;upper)`, p,
           "Significance" = p.signif) %>%
    mutate(Variant = gsub("D614G", "WT", Variant)) %>%
    arrange(p)-> titer_wilcox
  
  forest_data_comb %>%
    filter(both_assays) %>%
    filter(count_n > 1) %>%
    filter(shapiro_comb <= 0.05) %>%
    group_by(`Comparator antigen`, standardise_encounters) %>%
    wilcox_effsize(Titer~standardised_assay) %>%
    mutate(effsize = round(effsize, 2),
           `Comparator antigen` = gsub("D614G", "WT", `Comparator antigen`)) %>%
    select("Serum group" = standardise_encounters,
           "Effect size" = effsize,
           "Magnitude" = magnitude,
           "Variant" = `Comparator antigen`) %>%
    left_join(titer_wilcox, ., by = c("Serum group", "Variant")) -> titer_wilcox
  
  kable(titer_ttest, "latex") %>%
    row_spec(0,bold=TRUE)
  
  kable(titer_wilcox, "latex") %>%
    row_spec(0,bold=TRUE)
  
  stat_doc <- officer::read_docx(path = paste0(path_to_save,paste0("/stat_test_assay_FCfromWT.docx")))
  
  stat_doc %>%
    officer::body_add_break() %>%
    officer::body_add_par("Supplementary Table XX: 
                          Assay-based statistical comparison of live
virus vs pseudo virus variant GMTs.} A t-test was performed to compare
pseudovirus (PV) and live virus (LV) assessed Geometric Mean Titers
(GMTs) in different serum groups. Normality was checked with a
Shapiro-Wilk test and these serum groups were found to not differ
significantly from a normal distribution. Results are ordered by
increasing p-value. A statistic above 0 indicates higher values in Group
1 than in Group 2. The effect size was calculated with Cohen's d for unequal variances and Hedge's correction due to small sample sizes. All tests were performed with the rstatix package") %>%
    officer::body_add_table(.,titer_ttest, style = "table_template") %>%
    officer::body_add_break() %>%
    officer::body_add_par("Supplementary Table XX: Assay-based statistical comparison of live
virus vs pseudovirus variant GMTs.. 
                          A Wilcoxon-test was performed to
compare pseudovirus (PV) and live virus (LV) assessed Geometric Mean
Titers (GMTs) in different serum groups. Normality was checked with a
Shapiro-Wilk test and these serum groups were found to differ
significantly from a normal distribution. Results are ordered by
increasing p-value. A statistic above 0 indicates higher values in Group
1 than in Group 2. The effect size was calculated with Wilcoxon effect size. All tests were performed with the rstatix package ") %>%
    officer::body_add_table(.,titer_wilcox, style = "table_template") -> stat_doc
  
  print(stat_doc, target =
          paste0(path_to_save,paste0("/stat_test_assay_FCfromWT.docx")))
  
  
 
 }

