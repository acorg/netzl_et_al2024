# plot study by country of corresponding author to show data distribution
rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(rnaturalearthdata)
library(ggnewscale)
theme_set(theme_bw())


# to suppress NA removal warnings from ggplot
options(warn=-1)

working_dir <- getwd()

utility_dir <- file.path(working_dir, "code", "utility")
data_dir <- file.path(working_dir,'data')
figures_dir <- file.path(working_dir, "figures", "by_country")
google_sheets_dir <- file.path(data_dir, "google_sheet_tables")

source(file.path(utility_dir,'plot_functions_auto_label.R'))
source(file.path(utility_dir, "plot_over_time_functions.R"))
source(file.path(utility_dir, "prepare_table_for_forest_plots.R"))

#-------------------------------------  Read in studies table
studies <- read.csv(file.path(google_sheets_dir, "Netzl et al. - Collected Omicron antigenic data.csv")) %>%
  select(Study, Sourcelink, Country, year) %>%
  unique() %>%
  mutate(Country = gsub("The Netherlands", "Netherlands", Country)) %>%
  filter(Study != "Moore/Richardson") # not included in analysis as they had BA.4/5

studies <- standardise_time_format(studies, add_time_since_first = FALSE)

studies$standardise_time[studies$Study == "HKU"] <- as.Date("2021/12/12")
studies$standardise_time[studies$Study == "Krumbholz"] <- as.Date("2022/03/29")

world_coordinates <- map_data("world") 

# split multiple country codes
double_country_studies <- studies[grep(",", studies$Country),]

split_countries <- lapply(1:nrow(double_country_studies), function(x){
  
  temp_countries <- strsplit(double_country_studies[x, "Country"], ", ")[[1]]
  
  rbind(double_country_studies[x,] %>%
          mutate(Country = temp_countries[1]),
        double_country_studies[x,] %>%
          mutate(Country = temp_countries[2]))
  
})

studies <- rbind(studies[!grepl(",", studies$Country),],
                 do.call(rbind, split_countries))

# count how many studies per country
studies %>%
  select(Country) %>%
  group_by(Country) %>%
  count() -> country_count

# add it to world coordinates
world_coordinates$study_count <- country_count$n[match(world_coordinates$region, country_count$Country)]

# get coordinates for points per country
country_points <- world_coordinates %>%
  filter(region %in% country_count$Country) %>%
  select(region, long, lat) %>%
  group_by(region) %>%
  summarize(lat = mean(lat),
            long = mean(long)) %>%
  mutate(lat = ifelse(region == "USA", 40, lat),
         long = ifelse(region == "USA", -100, long))

# add coordinates + jitter to studies
get_jitter_country_coordinates <- function(world, country, count){
  
  world <- world %>%
    filter(region == country) %>%
    mutate(long = long + runif(count, min = -4, 4),
           lat = lat + runif(count, min = -4, 4)) + 
    select(long, lat)
  
  return(world)
  
}

studies$lat <- country_points$lat[match(studies$Country, country_points$region)]
studies$long <- country_points$long[match(studies$Country, country_points$region)]
studies$n <- country_count$n[match(studies$Country, country_count$Country)]
# add jitter 
studies %>%
  mutate(jitter = rnorm(nrow(studies), 0, n)) -> studies

time_legend <- unique(studies$standardise_time)
time_vals <- seq(from=min(time_legend), by=30, to=max(time_legend))
# create world map using ggplot() function 
ggplot() + 
  # geom_map() function takes world coordinates  
  # as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(long, lat, map_id = region, fill = study_count),
    alpha = 0.6
  ) + 
  scale_fill_steps(breaks = unique(country_count$n),
                   name = "# Reports per country") + 
  new_scale_fill() + 
  geom_jitter(data = studies %>%
                filter(n >= 4), aes(x=long, y= lat, group = Study, fill = standardise_time), width = 7, height = 5, size = 1.5, shape = 21, alpha = 0.7) +
  geom_jitter(data = studies %>%
                filter(n < 4), aes(x=long, y= lat, group = Study, fill = standardise_time), width = 2, height = 2, size = 1.5, shape = 21, alpha = 0.7) +
  scale_fill_viridis_c(breaks = time_vals,
                       labels = time_vals,
                       name = "Date of report") + 
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(legend.position = c(0.1, 0.7)) -> w_plot
 

ggsave(file.path(figures_dir, "studies_by_country.png"), w_plot, dpi = 300, width = 10, height = 7)
  
# make stacked bar chart

studies %>%
  mutate(Count = 1) %>%
  arrange(standardise_time) %>%
  ggplot(aes(x = Country, y = Count)) + 
  geom_col(aes(fill = standardise_time), color = "white") + 
  scale_fill_viridis_c(breaks = time_vals,
                       labels = time_vals,
                       name = "Date of report") + 
  geom_text(data = country_count, aes(y = n+0.5, label = n)) + 
  scale_x_discrete(limits = rev(country_count %>%
                                  arrange(n) %>%
                                  pull(Country))) +
  theme(legend.position = c(0.85, 0.85),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("# Reports") -> studies_bp

ggsave(file.path(figures_dir, "studies_by_country_bp.png"), studies_bp, dpi = 300, width = 6, height = 7)

