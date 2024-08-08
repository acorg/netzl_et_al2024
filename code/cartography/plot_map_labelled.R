# labelled map 
rm(list = ls())
library(Racmacs)

working_dir <- getwd()
map_dir <- file.path(working_dir, "data", "maps")
metadata_dir <- file.path(working_dir,'data', 'metadata')

roessler_map <- read.acmap(file.path(map_dir, "roessler_et_al_2023.ace"))
wilks_map <- read.acmap(file.path(map_dir, "Wilks_et_al_map_ndsubset_no_outliers_slope_adjusted.ace"))
muehlemann_map <- read.acmap(file.path(map_dir, "muehlemann_et_al_merged_duplicated_antigens_only.ace"))

agNames(roessler_map)[agNames(roessler_map) == "BA.5.3.2"] <- "BA.4/5"
agNames(roessler_map)[agNames(roessler_map) == "P.1.1"] <- "P.1"
agNames(muehlemann_map)[agNames(muehlemann_map) == "614D"] <- "WT"

# change color
agFill(wilks_map)[agNames(wilks_map) == "BA.2"] <- agFill(roessler_map)[agNames(roessler_map) == "BA.2"]

target_map <- "omicron_map_preprocessed_wSAVE_lessWD"

figure_dir <- file.path(working_dir, "figures","maps", target_map)
suppressWarnings(dir.create(file.path(figure_dir, "validation", "procrustes"), recursive = T))

xlim_zoom <- read.csv(file.path(metadata_dir, "xlim_zoom.csv"))$x
ylim_zoom <- read.csv(file.path(metadata_dir, "ylim_zoom.csv"))$x
ylim_zoom[1] <- ylim_zoom[1]

xlim <- read.csv(file.path(metadata_dir, "xlim_no_zoom.csv"))$x
ylim <- read.csv(file.path(metadata_dir, "ylim_no_zoom.csv"))$x
ylim[1] <- ylim[1] + 1


# Setup plotting function
doplot <- function(map, xlims, ylims, plot_labels = TRUE, plot_target_labels = NULL,
                   alpha_right = FALSE) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  # Plot the regular map
  srOutlineWidth(map) <- 1
  srSize(map) <- 6
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.6)
  
  agFill(map)<- adjustcolor(agFill(map), alpha.f = 0.9)
  
  plot(map, xlim = xlims, ylim =ylims, plot_stress = TRUE)
  
  if(plot_labels){
    x_adj <- 1
    y_adj <- 1
    # Plot labels
    label_adjustments <- matrix(0, numAntigens(map), 2)
    rownames(label_adjustments) <- agNames(map)
    label_adjustments["B.1.351",] <- c(x_adj+0.1, 0)
    label_adjustments["P.1",] <- c(0, y_adj-0.25)
    label_adjustments["B.1.1.7",] <- c(ifelse(alpha_right, x_adj, -x_adj), 0)
    label_adjustments["D614G",] <- c(-x_adj, 0)
    
    label_adjustments["B.1.617.2",] <- c(0,-y_adj+0.2)
    label_adjustments["BA.1",] <- c(0, -y_adj + 0.2)
    if("BA.2" %in% agNames(map)){
      label_adjustments["BA.2",] <- c(0, y_adj - 0.4)  
    }
    
    if("WT" %in% agNames(map)){
      label_adjustments["WT",] <- c(x_adj-0.3, 0)
    }
    
    labels <- rep("", length(agNames(map)))
    names(labels) <- agNames(map)
    if(is.null(plot_target_labels)){
      labels <- agNames(map)
      names(labels) <- agNames(map)
    } else {
      labels[plot_target_labels] <- plot_target_labels
      label_adjustments["P.1",] <- c(-x_adj, 0)
    }
    
    labels["B.1.351"] <- "Beta\n(B.1.351)"
    labels["P.1"] <- "Gamma\n(P.1)"
    labels["B.1.617.2"] <- "Delta\n(B.1.617.2)"
    labels["B.1.1.7"] <- "Alpha\n(B.1.1.7)"
    
    if("BA.1.1" %in% agNames(map) & "BA.4/5" %in% agNames(map)){
       label_adjustments["BA.1.1",] <- c(0, 0.5)
       label_adjustments["BA.4/5",] <- c(-0.7, 0)
       label_adjustments["BA.3",] <- c(0, 0.5)
       label_adjustments["BA.2.12.1",] <- c(0.8, 0)
       label_adjustments["BA.2.75",] <- c(0, 0.75)
    }
  
    if("BA.2.75" %in% agNames(map)){
      label_adjustments["BA.2.75",] <- c(0, 0.75)
    }
    label_size <- rep(1, numAntigens(map))
    names(label_size) <- agNames(map)
    
    
    
    text(
      agCoords(map) + label_adjustments,
      cex = label_size,
      label = labels,
      font = 1,
      col = "grey30"
    )
  }
 
  
}


for(map_specifics in c("_conv_ag_sub")){
  
  map_conv <- read.acmap(file.path(map_dir,paste0(target_map, map_specifics, '.ace')))
  map_early <- read.acmap(file.path(map_dir,paste0(gsub("map", "map_only30d", target_map), map_specifics, '.ace')))
  map_names <- list()
  
  map_names[[map_specifics]] <- map_conv
  
  if(grepl("sub", map_specifics)){
    
    map_pv <- read.acmap(file.path(map_dir, paste0(target_map, map_specifics, '_PV.ace')))
    map_lv <- read.acmap(file.path(map_dir, paste0(target_map, map_specifics, '_LV.ace')))
    
    map_names[[paste0(map_specifics, "_PV")]] <- map_pv
    map_names[[paste0(map_specifics, "_LV")]] <- map_lv
                      
    
  }
  
  
  png(file.path(figure_dir, paste0(map_specifics, "_early_comparison.png")), 10, 5, units = 'in', res=300, pointsize = 12)
  layout(matrix(c(1:2), ncol = 2, byrow = T))
  par(mar = rep(0.5, 4))
  doplot(map_early, xlim, ylim)
  doplot(procrustesMap(map_conv, map_early, sera = FALSE, translation = FALSE), xlim, ylim)
  dev.off()
  
  
  png(file.path(figure_dir, paste0(map_specifics, "_merged_comparison.png")), 13, 2.5, units = 'in', res=300, pointsize = 12)
  layout(matrix(c(1:4), ncol = 4, byrow = T))
  par(mar = rep(0.5, 4))
  doplot(map_early, xlim, ylim)
  doplot(procrustesMap(map_conv, map_early, sera = FALSE, translation = FALSE), xlim, ylim)
  doplot(procrustesMap(muehlemann_map, map_early, sera = FALSE, translation = TRUE, scaling = FALSE), xlim+2, ylim+1, plot_labels = TRUE, plot_target_labels = agNames(map_early))
  doplot(procrustesMap(muehlemann_map, map_conv, sera = FALSE, translation = TRUE, scaling = FALSE), xlim+2, ylim+1, plot_labels = TRUE, plot_target_labels = agNames(map_conv))
  dev.off()
  
  png(file.path(figure_dir, paste0(map_specifics, "_wilks_roessler_comparison.png")), 14, 3, units = 'in', res=300, pointsize = 12)
  layout(matrix(c(1:4), ncol = 4, byrow = T))
  par(mar = rep(0.5, 4))
  doplot(map_early, xlim, ylim)
  doplot(procrustesMap(map_conv, map_early, sera = FALSE, translation = FALSE), xlim, ylim)
  doplot(procrustesMap(wilks_map, map_conv, sera = FALSE, translation = TRUE, scaling = FALSE), xlim+1, ylim, plot_labels = TRUE, plot_target_labels = agNames(map_conv)[agNames(map_conv) != "WT"], alpha_right = TRUE)
  doplot(procrustesMap(roessler_map, map_conv, sera = FALSE, translation = TRUE, scaling = FALSE), xlim+1, ylim, plot_labels = TRUE, plot_target_labels = agNames(map_conv)[agNames(map_conv) != "WT"], alpha_right = TRUE)
  dev.off()
  
  
  
  # do procrustes of PV, LV maps
  png(file.path(figure_dir,"validation", "procrustes", paste0(map_specifics, "_procrustes_assay.png")), 10, 6, units = 'in', res=300, pointsize = 12)
  layout(matrix(c(1:6), ncol = 3, byrow = T))
  par(mar = rep(0.5, 4))
  doplot(map_lv, xlim, ylim, plot_labels = FALSE)
  text(xlim[1]+0.4, ylim[2]-0.4, "A", cex = 2)
  doplot(map_pv, xlim, ylim, plot_labels = FALSE)
  text(xlim[1]+0.4, ylim[2]-0.4, "B", cex = 2)
  doplot(procrustesMap(map_lv, map_pv), xlim, ylim, plot_labels = FALSE)
  text(xlim[1]+0.4, ylim[2]-0.4, "C", cex = 2)
  doplot(procrustesMap(map_lv, map_conv, sera = FALSE), xlim, ylim, plot_labels = FALSE)
  text(xlim[1]+0.4, ylim[2]-0.4, "D", cex = 2)
  doplot(procrustesMap(map_pv, map_conv, sera = FALSE), xlim, ylim, plot_labels = FALSE)
  text(xlim[1]+0.4, ylim[2]-0.4, "E", cex = 2)
  dev.off()
  
  
  
  
}


