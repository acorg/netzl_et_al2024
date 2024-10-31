# plot study by country of corresponding author to show data distribution
rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
theme_set(theme_bw())


# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "by_sex")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "plot_over_time_functions.R"))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))

#-------------------------------------  Read in studies table
studies <- read.csv(file.path(google_sheets_dir, "Netzl et al. - Collected Omicron antigenic data.csv")) %>%
  select(Study, Sourcelink, Country, year, Number.of.sera, N.female., Sera.details.long) %>%
  unique() %>%
  mutate(Country = gsub("The Netherlands", "Netherlands", Country)) %>%
  filter(Study != "Moore/Richardson") # not included in analysis as they had BA.4/5

studies <- standardise_time_format(studies, add_time_since_first = FALSE)

studies$Study[grep("nch", studies$Study)] <- "Münch"

studies$standardise_time[studies$Study == "HKU"] <- as.Date("2021/12/12")
studies$standardise_time[studies$Study == "Krumbholz"] <- as.Date("2022/03/29")

studies_f <- studies %>%
  mutate(frac_female = as.numeric(N.female.)/as.numeric(Number.of.sera)) %>%
  filter(!is.na(frac_female)) %>%
  group_by(Study) %>%
  mutate(mean_frac = mean(frac_female)) %>%
  ungroup() %>%
  arrange(mean_frac)

studies_f$Study <- factor(studies_f$Study, levels = unique(studies_f$Study))

studies_f %>%
  ggplot(aes(x = Study, y = frac_female)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_point(aes(group = Sera.details.long), position = position_dodge(width = 0.6), alpha = 0.6) +
  geom_point(aes(y = mean_frac), shape = 21, fill = "transparent", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylim(c(0,1)) + 
  ylab("Fraction of female subjects") -> p

ggsave(file.path(figures_dir, "fraction_of_female.png"), p, dpi = 300, width = 12, height = 5)


studies %>%
  filter(is.na(N.female.)) %>%
  select(Sera.details.long) %>%
  unique() %>%
  count()

studies %>%
  filter(is.na(N.female.)) %>%
  select(Study) %>%
  unique() %>%
  count()
