# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
set.seed(Sys.time())

working_dir <- getwd()
maps_dir <- file.path(working_dir,'data/maps')
bootstraps_dir <- maps_dir

base_map <- 'omicron_map_preprocessed_wSAVE_lessWD'
map_names <- list(paste0(base_map, '_conv_ag_sub'))

bootstrap_repeats <- 500
optimizations_per_repeat <- 100

for (map_name in map_names){
  
  map_path = file.path(maps_dir, paste0(map_name, '.ace'))
  
  map <- read.acmap(map_path)

  for (method in list('bayesian', 'noisy', 'resample')){
    
    map_with_bootstrap_data = bootstrapMap(
                                map,
                                method=method,
                                bootstrap_repeats = bootstrap_repeats,
                                bootstrap_ags = TRUE,
                                bootstrap_sr = TRUE,
                                reoptimize = TRUE,
                                optimizations_per_repeat = optimizations_per_repeat,
                                ag_noise_sd = 0.7,
                                titer_noise_sd = 0.7,
                                options = list()
                              )
  
    ext <- paste0(method,'_bootstrap.ace')
    bootstrap_map_name <- paste(map_name, ext, sep='_')
    bootstrap_path <- file.path(bootstraps_dir,base_map, 'validation', bootstrap_map_name)
    
    save.acmap(map_with_bootstrap_data, bootstrap_path)
    ## there is a plot_bootstrap_blobs.R script inside code/cartography/validation
    ## that does the actual plotting
    
  }

}