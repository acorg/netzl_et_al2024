get_plot_heights <- function(data, group_by, max_plot_height = 14){
  
  min_plot_height <- 2
  max_plot_height <- max_plot_height
  
  data$`Comparator antigen`[data$`Comparator antigen` == "D614G"] <- "WT"
  
  # take ba.1 WT standardise encounters as standard height
  data %>% filter(OmicronVariant == "BA.1") %>%
    filter(`Comparator antigen` %in% c("WT", "D614G")) %>%
    count() -> base_count
  
  base_count <- base_count + length(unique(data$standardise_encounters))
  
  step_by <- (max_plot_height - min_plot_height)/base_count
  
  n_data_points <- data %>% group_by(OmicronVariant, `Comparator antigen`) %>% 
    summarize(n = n()) %>%
    ungroup()
  
  for(row in 1:nrow(n_data_points)) {

    n_unique <- 0
    for(grouping in group_by) {
      omicron_var <- n_data_points[[row, "OmicronVariant"]]
      comp_antigen <- n_data_points[[row, "Comparator antigen"]]
      
      n_unique <- n_unique + length(unique(data[data$OmicronVariant == omicron_var & data$`Comparator antigen` == comp_antigen,grouping])) -1
    }
    n_data_points$n[row] <- n_data_points$n[row] + n_unique
  }
  
  
  n_data_points$figure_height <- unlist(lapply(n_data_points$n, function(x) x*step_by + min_plot_height))
  
  return(n_data_points)                                              
   
}

plot_forest_plot_by_sublineage <- function(forest_data, sublineage = c("BA.1"), group_by = c("standardise_encounters"),
                                           xmin = -8, xmax = 4, point_size = 1, axis_text_size = 2.5, path_save = "",
                                           max_plot_height = 14) {
  
  plot_name <- paste0(paste0(sublineage, collapse = "-"), "_", paste0(group_by, collapse = "-"))
  plot_name <- gsub("BA.4\\/5", "BA.4.5", plot_name)
  plot_name_mean <- paste0(plot_name, "_mean")
  
  
  plot_data <- forest_data %>% filter(OmicronVariant %in% sublineage) %>% unique() %>%
    character_rowlabel()
  
  # get height for subplots from forest data
  plot_heights <- get_plot_heights(forest_data, group_by, max_plot_height) %>%
    filter(OmicronVariant == sublineage) %>%
    column_to_rownames("Comparator antigen") %>%
    select(figure_height)
  
  
  # order by comparator antigen
  plot_data <- reorder_data(plot_data[order(plot_data$`Comparator antigen`),], rev = TRUE)
  
  #  order by overall drop
  plot_data <- reorder_data(plot_data[order(plot_data$log_fold_change),], rev = TRUE)
  
  # order by fold drop within serum groups
  plot_data <- reorder_data(plot_data[order(plot_data$standardise_encounters),], rev = FALSE)
  
  save_in_scaled_format(plot_data, which_plot = "titer_drop", to_save = plot_name,
                        hline_by = group_by, show_mean = FALSE, single_plots = FALSE,
                        axis_text_size = axis_text_size,
                        height_wt = plot_heights["WT", "figure_height"],
                        height_alpha =plot_heights["Alpha", "figure_height"],
                        height_beta = plot_heights["Beta", "figure_height"],
                        height_delta = plot_heights["Delta", "figure_height"],
                        path_to_save = file.path(path_save, "fold_drops"))

  save_in_scaled_format(plot_data, which_plot = "titer_drop", to_save = plot_name_mean,
                        hline_by = group_by, show_mean = TRUE, single_plots = FALSE,
                        axis_text_size =axis_text_size,
                        height_wt = plot_heights["WT", "figure_height"],
                        height_alpha =plot_heights["Alpha", "figure_height"],
                        height_beta = plot_heights["Beta", "figure_height"],
                        height_delta = plot_heights["Delta", "figure_height"],
                        path_to_save = file.path(path_save, "fold_drops"))

  # do titer forest plots here
  plot_data_titer <- plot_data %>% filter(!is.na(TitersHAg))
  
  # remove here where Webplotdigitzer was used for titers
  if(grepl("lessWD", path_save)){
    plot_data_titer %>% 
      filter(Webplotdigitizer == "n") -> plot_data_titer
  }
  
  # increase plot heights
  plot_heights$figure_height <- plot_heights$figure_height + 0.5
  
  save_in_scaled_format(plot_data_titer, which_plot = "not_titer_drop", to_save = plot_name, 
                        hline_by = group_by, show_mean = FALSE, single_plots = FALSE, 
                        axis_text_size = axis_text_size, 
                        height_wt = plot_heights["WT", "figure_height"],
                        height_alpha =plot_heights["Alpha", "figure_height"],
                        height_beta = plot_heights["Beta", "figure_height"],
                        height_delta = plot_heights["Delta", "figure_height"],
                        path_to_save = file.path(path_save, "titers"),
                        order_by_omicron = F)
  
  
  plot_name_mean <- paste0(plot_name, "_mean")
  
  save_in_scaled_format(plot_data_titer, which_plot = "not_titer_drop", to_save = plot_name_mean,
                        hline_by = group_by, show_mean = TRUE, single_plots = FALSE,
                        axis_text_size = axis_text_size,
                        height_wt = plot_heights["WT", "figure_height"],
                        height_alpha =plot_heights["Alpha", "figure_height"],
                        height_beta = plot_heights["Beta", "figure_height"],
                        height_delta = plot_heights["Delta", "figure_height"],
                        path_to_save = file.path(path_save, "titers"),
                        order_by_omicron = F)
  
}