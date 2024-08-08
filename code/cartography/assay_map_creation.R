##' Antigenic cartography, do authentic and pseudovirus maps

# setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
set.seed(100)


# ---------------- set up directories
working_dir <- getwd()

utility_dir <- file.path(working_dir, "code/utility")
data_dir <- file.path(working_dir,'data')
map_dir <- file.path(data_dir, "maps")

#----------------- Create maps for PV and LV antigen assays
# here for multiple tables
map_names <- list('omicron_map_preprocessed_wSAVE_lessWD_conv_ag_sub.ace')

for (map_name in map_names){
  
  
  map_full <- read.acmap(file.path(map_dir, map_name))
  alignment_map <- map_full
  
  for(assay in c("LV", "PV")) {
    map_assay_name <- gsub(".ace", paste0("_", assay,".ace"), map_name)
    
    sr_names <- srNames(map_full)
    if(assay == "LV") {
      target_names <-  sr_names[grepl("_LV_", sr_names)]
    } else {
      target_names <-  sr_names[grepl("_PV_", sr_names)]
    }
    
    map_conv <- subsetMap(map_full, antigens = TRUE, sera = target_names)
    
    map_conv <- optimizeMap(map_conv, number_of_dimensions = 2, number_of_optimizations = 1000,
                            options = list(ignore_disconnected = TRUE))
    
    dilutionStepsize(map_conv) <- 0
    map_conv <- realignMap(map_conv, alignment_map)
    
    save.acmap(map_conv, filename = file.path(map_dir, map_assay_name))
    
  }
  
  
}

