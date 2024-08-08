#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(ablandscapes)
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

# ----------------------- set directories and load functions
working_dir <- getwd()

code_dir <- file.path(working_dir, "code")
data_dir <- file.path(working_dir, "data")
map_dir <- file.path(data_dir, "maps")
figure_dir_b <- file.path(working_dir, "figures", "landscapes", "gmt_landscapes")


source(file.path(code_dir, "utility", 'map_longinfo.R'))
source(file.path(code_dir, "utility", 'landscape_functions_util.R'))


# ---------------------- do fitting per map

# if you want to exclude some antigens from the fit
# not positioned in base map
ags_to_exclude <- c("BA.2.75")

map_names <- list('omicron_map_preprocessed_wSAVE_lessWD_conv_ag_sub.ace',
                  'omicron_map_only30d_preprocessed_wSAVE_lessWD_conv_ag_sub.ace',
                  'omicron_map_only90d_preprocessed_wSAVE_lessWD_conv_ag_sub.ace')


sr_colors <- read.csv(file = file.path(data_dir,"metadata", "map-colors.csv"), header = TRUE, stringsAsFactors = FALSE, sep = ",", row.names = "Variable")


for(map_name in map_names) {
  
  figure_dir <- file.path(figure_dir_b, gsub(".ace", "", map_name))
  suppressWarnings(dir.create(figure_dir, recursive = T))
  
  # set name for fit
  fit_name <- gsub("_conv_ag_sub", "_conv_ag_sub_fit", map_name)
  fit_name <- gsub(".ace", ".rds", map_name)
  
  # read base map
  map <- read.acmap(file.path(map_dir, map_name))
  base_map_sr_groups <- unique(as.character(srGroups(map)))
  
  # read the full map and subset to sr groups that are not in base map
  full_map_name <- gsub("conv_ag_sub", "full", map_name)
  full_map <- read.acmap(file.path(map_dir, full_map_name))

    # Remove ags from fit
  map <- removeAntigens(map, ags_to_exclude)
  lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
  
  ags_to_fit_lndscp <- agNames(map)
  
  ags_to_fit_lndscp
 
  # read the full map
  map_orig <- full_map
  
  # subset the map target sera and ags
   map_long <- long_map_info(map_orig) %>%
    filter(!(ag_name %in% ags_to_exclude))
  
  map_long %>%
    select(titer, ag_name, sr_name, sr_group) -> titerdata
  
  titerdata %>%
    group_by(
      sr_group
    ) -> titerdata
  
  titerdata %>%
    group_map(
      get_titertable
    ) -> titertables
  
 
  lndscp_fits <- lapply(
    titertables,
    
    function(titertable) {
      
       if(!is.null(dim(titertable))){
         
         ablandscape.fit(
           titers = titertable[,ags_to_fit_lndscp, drop = TRUE],
           bandwidth = 1,
           degree = 1,
           method = "cone",
           error.sd = 1,
           acmap = map,
           control = list(
             optimise.cone.slope = TRUE
           )
         )
       }
      
      
      
    }
  )
  
  titertables_groups <- group_data(titerdata)

  
  # Add impulses
  titerdata %>%
    group_by(
      sr_group,
      ag_name
    ) %>%
    summarize(gmt = mean(log2(as.numeric(titer)/10), na.rm = TRUE)) -> gmt_data
  
  # angle for html page
  angle <- list(
    rotation = c(-1.5329, -0.0015, -0.0222), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
    translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
    zoom = 2
    # zoom = 1.1646 # higher is more zoomed out
  )
  
  
  data3js <- base_plot_data3js(map, lndscp_fits[2], agNames(map), lims, agNames(map))
  
  
  lndscp_list <- list()
  # set the target exposure groups
  exposure_sr_groups <- c("2x Vax", "3x Vax")
  
  target_ind <- match(exposure_sr_groups, titertables_groups$sr_group)
  
  lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups[target_ind,], lndscp_fits[target_ind], map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors, 
                                          show_gmts = TRUE,
                                          hide_buttons = TRUE)
  lndscp <-r3js(
    lndscp_3js,
    rotation = angle$rotation,
    zoom = angle$zoom
  )
  
  lndscp_list[["gmt_vax"]] <- lndscp
  
  save_name <- file.path(figure_dir, "main_gmt_vax.png")
  plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)
  
  # set the target exposure groups
  
  if(!grepl("30d", fit_name)){
    exposure_sr_groups <- c("Vax + BA.1", "2x Vax", "3x Vax")
    
    target_ind <- match(exposure_sr_groups, titertables_groups$sr_group)
    
    lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups[target_ind,], lndscp_fits[target_ind], map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors, show_gmts = TRUE,
                                            hide_buttons = TRUE)
    lndscp <-r3js(
      lndscp_3js,
      rotation = angle$rotation,
      zoom = angle$zoom
    )
    
    save_name <- file.path(figure_dir, "main_gmt_inf.png")
    plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)
    
    lndscp_list[["breakthrough"]] <- lndscp
  }
  
  
}

