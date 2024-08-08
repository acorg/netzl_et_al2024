
# determine plot heights for boxplot
boxplot_widths_heights <- function(data, x_var, facet_var) {
  
  max_width <- 8
  row_height <- 3
  
  points_per_row <- 40
  
  point_width <- max_width/points_per_row
  
  total <- length(unique(data[,x_var])) * length(unique(data[,facet_var]))
  nrow <- ceiling(total/points_per_row)
  
  plot_height <- row_height * nrow
  plot_width <- max_width * (total/(points_per_row*nrow))
  
  return(list("height" = plot_height, "width" = plot_width, "rows" = nrow))
  
}

# boxplot of fold drops
do_boxplot_fold_change <- function(data, comp_antigen, x_var, facet_var, x_label, nrows = 2){

  data$standardise_encounters <- factor(as.character(data$standardise_encounters), 
                                        levels = rev(levels(data$standardise_encounters)))
  data$standardise_encounters <- gsub("WT", "First wave", data$standardise_encounters)
  
  xmax <- floor(min(data$log_fold_change, na.rm = T))
  xmin <- ceiling(max(data$log_fold_change, na.rm = T))
 
  data$arrow_length <- as.character(data$arrow_length)
  
  
  colors <- c(colors, "1"= "black","0.5"= "grey40")

  data %>%
    ggplot(aes_string(y = "log_fold_change", x = x_var, color = "OmicronVariant")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_jitter(aes(color = ifelse(is.na(arrow_length), OmicronVariant, arrow_length)), alpha = 0.4) +
    geom_boxplot() +
    scale_color_manual(values = colors) +
    facet_wrap(
      facets = as.formula(paste("~", facet_var)),
      nrow = nrows
    ) +
    xlab(x_label) + 
    scale_y_reverse(
      name = paste0("Fold change from ", gsub("WT", "Prototype", comp_antigen)), 
      labels = function(x) ifelse(x > 0, paste0("-",2^x),  2^abs(x)), 
      breaks = c(xmin:xmax)
    ) +
    theme_bw() +
    theme(legend.position = "none", 
          strip.background =element_rect(fill="white"),
          axis.text.x = element_text(hjust=1, angle=45)) -> bp
  
  return(bp)
}

# boxplot of titers
do_boxplot_titer <- function(data, comp_antigen, x_var, facet_var, x_label, nrows = 2){
  
  
  data$standardise_encounters <- factor(as.character(data$standardise_encounters), 
                                        levels = rev(levels(data$standardise_encounters)))
  data$standardise_encounters <- gsub("WT", "First wave", data$standardise_encounters)
  
  xmin <- floor(min(data$Log2Omi, na.rm = T))
  xmax <- ceiling(max(data$Log2Omi, na.rm = T))
  
  colors <- c(colors, "1"= "black","0.5"= "grey40")
  data %>%
    ggplot(aes_string(y = "Log2Omi", x = x_var, color = "OmicronVariant")) +
    geom_jitter(aes(color = ifelse(is.na(arrow_length), OmicronVariant, arrow_length)), alpha = 0.4) +
    geom_boxplot() +
    scale_color_manual(values = colors) +
    facet_wrap(
      facets = as.formula(paste("~", facet_var)),
      nrow = nrows
    ) +
    xlab(x_label) + 
    scale_y_continuous(
      name = "Titer", 
      labels = function(x) 2^x*10, 
      breaks = c(xmin:xmax)
    ) +
    theme_bw() +
    theme(legend.position = "none", 
          strip.background =element_rect(fill="white"),
          axis.text.x = element_text(hjust=1, angle=45)) -> bp
  
  return(bp)
}


save_boxplots <- function(data, comp_antigen, x_var, facet_var, x_label, to_save, plot_prefix = "", do_titer_plot = F) {
  
  data$`Comparator antigen`[data$`Comparator antigen` %in% c("WT", "D614G")] <- "WT"
  
  data <- data %>% filter(`Comparator antigen` == comp_antigen)
  
  fig_paras <- boxplot_widths_heights(data, x_var, facet_var)
  
  plot_name <- paste0(plot_prefix, x_var, "_by_", facet_var, ".png")
  
  fc_plot <- do_boxplot_fold_change(data, comp_antigen, x_var, facet_var, x_label, nrows = fig_paras$rows)
  titer_plot <- do_boxplot_titer(data, comp_antigen, x_var, facet_var, x_label, nrows = fig_paras$rows)
  
  ggsave(filename = file.path(to_save, paste0("fold_change-", comp_antigen, "-", plot_name)), 
        plot = fc_plot, dpi = 300, width = fig_paras$width, height = fig_paras$height)
  
  if(do_titer_plot) {
    ggsave(filename = file.path(to_save, paste0("titers-", plot_name)), 
           plot = titer_plot, dpi = 300, width = fig_paras$width, height = fig_paras$height)
  }

}
  
  
  
  
  
  
  
  
  