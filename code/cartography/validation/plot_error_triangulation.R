# setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
set.seed(100)

# simple functions

remove_na_sera <- function(map){
  
  na_sera <- srCoords(map)
  na_sera <- rownames(na_sera)[is.na(na_sera)[,1]]
  
  map <- removeSera(map, na_sera)
  
  return(map)
}

remove_na_antigens <- function(map){
  
  na_sera <- agCoords(map)
  na_sera <- rownames(na_sera)[is.na(na_sera)[,1]]
  
  map <- removeAntigens(map, na_sera)
  
  return(map)
}

remove_na_coords <- function(map){
  
  map <- remove_na_sera(map)
  map <- remove_na_antigens(map)
  
  return(map)
  
}

# ---------------- set up directories
working_dir <- getwd()

utility_dir <- file.path(working_dir, "code","utility")
data_dir <- file.path(working_dir,'data')
map_dir <- file.path(data_dir, "maps")


xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x


base_map_names <- list('omicron_map_preprocessed_wSAVE_lessWD')

for(base_map in base_map_names){
  
  figure_dir <- file.path(working_dir, "figures", "maps", base_map, "validation", "error_triangulation")
  suppressWarnings(dir.create(figure_dir, recursive = T))
  
  # here for multiple tables
  map_names <- list(paste0(base_map, '_conv_ag_sub.ace'))
  
  for (map_name in map_names){
    
    
    map <- remove_na_coords(read.acmap(file.path(map_dir, map_name)))
    
    fig_name <- gsub(".ace", ".png", map_name)
    
    png(file.path(figure_dir, paste0("error_triangulation_", fig_name)), width = 6.5, height = 3, units = 'in', res=300, pointsize = 18)
    layout(matrix(c(1:2), ncol = 2, byrow = T))
    par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))
    plot(map, show_error = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
         grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
    text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)
    plot(triangulationBlobs(relaxMap(map), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
         grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
    text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)
    
    dev.off()
    
  }
  
  
}


