# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
set.seed(Sys.time())

working_dir <- getwd()

metadata_dir <- file.path(working_dir,'data', 'metadata')

xlim <- read.csv(file.path(metadata_dir, "xlim_no_zoom.csv"))$x
ylim <- read.csv(file.path(metadata_dir, "ylim_no_zoom.csv"))$x

base_map <- 'omicron_map_preprocessed_wSAVE_lessWD'
bootstraps_dir <- file.path(working_dir,'data','maps', base_map, 'validation')

map_names <- list.files(path = bootstraps_dir, 
                        pattern = ".ace")

bootstrap_blobs_figure_dir <- file.path('figures', 'maps', base_map, 'validation')


confidence_levels = list(0.6)
smoothing = 6
gridspacing = 0.05



  for (map_name in map_names){  
    
    bootstrap_map_path = file.path(bootstraps_dir, map_name)
    
    for (confidence_level in confidence_levels){
      
      map_with_bootstrap_data = read.acmap(bootstrap_map_path)
      
      map_with_blobs = bootstrapBlobs(
          map_with_bootstrap_data,
          conf.level = confidence_level,
          smoothing = smoothing,
          gridspacing = gridspacing,
          method = "ks"
        )
      
      
      bootstrap_plot_name <- gsub(".ace", ".png", map_name)
  
      bootstrap_blobs_figure_path <- file.path(bootstrap_blobs_figure_dir, 
                                              bootstrap_plot_name)
      
      
      png(file=bootstrap_blobs_figure_path, width = 4, height = 4, units = 'in', res=300)
      par(oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
      plot(map_with_blobs, xlim = xlim, ylim = ylim)
      dev.off()
     
      
    }
        
  }

