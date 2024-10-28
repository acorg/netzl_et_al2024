ag_colors <- read.csv(file.path("data", "metadata", "ag_colors.csv"))

standardise_time_format <- function(forest_data, group_by_sr_group = FALSE, group_not_by_comparator= FALSE, group_by_omicron = TRUE, add_time_since_first = TRUE){
 
   for(r in 1:nrow(forest_data)){
    test <- forest_data[r,]
    temp_year <- test$year
    if(grepl("\\.", test$year)){
      temp_year <- strsplit(test$year, "\\.")[[1]][1]
    } 
    
    temp <- str_locate(test$Sourcelink,temp_year)
    
    if(!is.na(temp[1])){
      forest_data$year[r] <- substr(test$Sourcelink, temp[1], temp[1]+9)  
    }
    
  }
  
  forest_data$standardise_time <- as.Date(forest_data$year, format = "%Y.%m.%d")
  
  if(!add_time_since_first){
    return(forest_data)
  }
  
  if(group_by_sr_group){
    
    if(group_by_omicron){
      
      if(group_not_by_comparator){
        
        forest_data %>%
          group_by(standardise_encounters, OmicronVariant) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data
        
      } else {
        forest_data %>%
          group_by(OmicronVariant, standardise_encounters, `Comparator antigen`) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data      
      }
      
    } else {
      
      if(group_not_by_comparator){
        
        forest_data %>%
          group_by(standardise_encounters) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data
        
      } else {
        forest_data %>%
          group_by(`Comparator antigen`, standardise_encounters) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data      
      }
    }
    
   

    
  } else {
    
    if(group_by_omicron){
      
      if(group_not_by_comparator){
        forest_data %>%
          group_by(OmicronVariant) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data  
      } else {
        
        forest_data %>%
          group_by(OmicronVariant, `Comparator antigen`) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data  
      }
    } else {
      
      if(group_not_by_comparator){
        forest_data %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data  
      } else {
        
        forest_data %>%
          group_by(`Comparator antigen`) %>%
          mutate(first_record = min(standardise_time, na.rm = TRUE),
                 time_since_first = standardise_time - first_record, 
                 time_since_first = as.numeric(time_since_first)) %>%
          ungroup() -> forest_data  
      }
      
    }
    
    
  }
 
  
  return(forest_data)
}

add_time_group <- function(forest_data_sub, dates, since_first = FALSE){
  
  forest_data_sub$date_group <- NA
  
  if(since_first){
    for(r in 1:nrow(forest_data_sub)){
      forest_data_sub$date_group[r] <- max(dates[dates <= forest_data_sub$time_since_first[r]])
    }
    
  
  } else {
   
    for(r in 1:nrow(forest_data_sub)){
      forest_data_sub$date_group[r] <- max(dates[dates < forest_data_sub$standardise_time[r]])
    }
    
    forest_data_sub$date_group<- sapply(forest_data_sub$standardise_time, function(x){
      max_d <- max(dates[dates < x])
      unique(grep(max_d, dates))
    })
  }
  
  return(forest_data_sub)
}

calculate_mean_fold_change_over_time <- function(forest_data_sub, dates){
  
  forest_data_sub %>%
    group_by(date_group, OmicronVariant, standardise_encounters, `Comparator antigen`) %>%
    summarize(n = length(log_fold_change),
              lower = Rmisc::CI(log_fold_change)["lower"],
              upper = Rmisc::CI(log_fold_change)["upper"],
              log_fold_change = mean(log_fold_change, na.rm = T)) %>%
    mutate(date_group = as.numeric(date_group),
           standardise_time = dates[date_group],
           shape = "GMT",
           y = 7)-> forest_data_plot
  
  return(forest_data_plot)
}

calculate_mean_over_time <- function(forest_data_sub, dates, target_var = "log_fold_change", since_first = FALSE){
  
  var <- sym(target_var)
  
  forest_data_sub %>%
    group_by(date_group, OmicronVariant, standardise_encounters, `Comparator antigen`) %>%
    summarize(n = length(!!var),
              lower = Rmisc::CI(!!var)["lower"],
              upper = Rmisc::CI(!!var)["upper"],
              log_fold_change = mean(!!var, na.rm = T)) %>%
    mutate(date_group = as.numeric(date_group),
           time_since_first = date_group,
           shape = "GMT",
           y = 7)-> forest_data_plot
  
  if(!since_first){
    forest_data_plot <- forest_data_plot %>%
      mutate(standardise_time = dates[date_group])
  }
  return(forest_data_plot)
}

calculate_cumulative_mean <- function(forest_data_sub, dates, target_var = "log_fold_change") {
  
  var <- sym(target_var)
  # now make it cumulative mean
  cumulative_means <- list()
  for(date in 1:length(dates)){
    
    forest_data_sub %>%
      filter(date_group <= dates[date]) %>%
      group_by(OmicronVariant, standardise_encounters, `Comparator antigen`) %>%
      summarize(n = length(!!var),
                lower = Rmisc::CI(!!var)["lower"],
                upper = Rmisc::CI(!!var)["upper"],
                log_fold_change = mean(!!var, na.rm = T)) %>%
      mutate(standardise_time = dates[date],
             shape = "GMT",
             y = 7,
             date_group = date,
             time_since_first = standardise_time)-> temp_plot
    
    cumulative_means[[date]] <- temp_plot
    
  }
  
  cum_means <- do.call(rbind, cumulative_means)
  
  return(cum_means)
}

plot_over_time_stepwise_mean <- function(forest_data_sub, forest_data_plot, dates, comp_antigen, sr_groups, ymin = -8, ymax = 10, stepwise = TRUE,
                                         var = "log_fold_change", plot_counts = TRUE){
  
  forest_data_sub <- forest_data_sub %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% sr_groups)
  
  forest_data_plot <- forest_data_plot %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% sr_groups)
  
  cols <- ag_colors %>%
    filter(Antigen %in% forest_data_plot$OmicronVariant)
  
  plot_colors <- cols$Color
  names(plot_colors) <- cols$Antigen

  forest_data_plot %>%
    ungroup() %>%
    mutate(y = ymax - 0.5,
           y_offset = (((as.numeric(date_group)+1) %% 2)/2),
           y = y -y_offset) -> forest_data_plot
  
  # forest_data_plot$y <- rep(c(ymax - 0.5, ymax-1), nrow(forest_data_plot)/2)[1:nrow(forest_data_plot)]

  forest_data_sub %>%
    ggplot(aes(x = standardise_time, y = log_fold_change, color = OmicronVariant)) + 
    geom_point(aes(shape = shape),  alpha = 0.7, color = "grey20", size = 1) -> gp
  
  if(stepwise){
    gp + stat_smooth(method = "lm", se=FALSE, color="black") -> gp 
  }
  
    gp + 
    geom_pointrange(data = forest_data_plot, aes(ymin = lower, ymax = upper, shape = shape), alpha = 1) -> gp
    
  if(plot_counts){
    gp + geom_text(data = forest_data_plot, aes(y = y, label = n), colour = "black", size = 2) -> gp
  }
    

  gp +    
    facet_grid(standardise_encounters~factor(OmicronVariant,
               levels = c("BA.1", "BA.1.1", "BA.2", "BA.3", "BA.4/5", "BA.2.12.1", "BA.2.75"))) + 
    scale_x_continuous(breaks = dates[seq(1, length(dates), by=4)],
                       name = "Date") +
    scale_shape_manual(values = c("sample" = 20, "GMT" = 21)) +
    scale_color_manual(values = plot_colors) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.background.x = element_rect(fill = "white"),
          strip.background.y = element_rect(fill = "white")) -> gp
    
    if(!stepwise){
      gp + 
        geom_line(data = forest_data_plot) -> gp
    }
  
    if(var == "log_fold_change") {
      gp <- gp + 
        scale_y_continuous(limits = c(ymin,ymax),
                           name = paste0("Fold change from ", gsub("D614G", "614D/614G", comp_antigen)),
                           breaks = seq(ymin, ymax-1, by = 2),
                           labels = function(x) ifelse(x < 0, paste0("-", 2^abs(x), "x"), paste0(2^x, "x")))
        
    } else {
      
      gp <- gp + 
        scale_y_continuous(limits = c(ymin,ymax),
                           name = "Titer",
                           breaks = seq(ymin, ymax-1, by = 2),
                           labels = function(x) 2^x*10)
      
    }
    
  return(gp)
  
}


plot_over_time_stepwise_cum <- function(forest_data_sub, forest_data_plot, cum_means, dates, comp_antigen, sr_groups, ymin = -8, ymax = 10, stepwise = TRUE,
                                         var = "log_fold_change", plot_counts = TRUE){
  
  forest_data_sub <- forest_data_sub %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% sr_groups)
  
  forest_data_plot <- forest_data_plot %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% sr_groups)
  
  cum_means <- cum_means %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% sr_groups)
  
  cols <- ag_colors %>%
    filter(Antigen %in% forest_data_plot$OmicronVariant)
  
  plot_colors <- cols$Color
  names(plot_colors) <- cols$Antigen
  
  forest_data_plot %>%
    ungroup() %>%
    mutate(y = ymax - 0.5,
           y_offset = (((as.numeric(date_group)+1) %% 2)/2),
           y = y -y_offset) -> forest_data_plot
  
  # forest_data_plot$y <- rep(c(ymax - 0.5, ymax-1), nrow(forest_data_plot)/2)[1:nrow(forest_data_plot)]
  
  forest_data_sub %>%
    ggplot(aes(x = time_since_first, y = log_fold_change, color = OmicronVariant)) + 
    geom_point(aes(shape = shape),  alpha = 0.7, color = "grey20", size = 1) + 
    geom_pointrange(data = forest_data_plot, aes(ymin = lower, ymax = upper, shape = shape), alpha = 1) -> gp
  
  if(plot_counts){
    gp + geom_text(data = forest_data_plot, aes(y = y, label = n), colour = "black", size = 2) -> gp
  }
  
  
  gp +    
    facet_grid(standardise_encounters~factor(OmicronVariant,
               levels = c("BA.1", "BA.1.1", "BA.2", "BA.3", "BA.4/5", "BA.2.12.1", "BA.2.75"))) + 
    scale_x_continuous(name = "Time since first record (days)") +
    scale_shape_manual(values = c("sample" = 20, "GMT" = 21)) +
    scale_color_manual(values = plot_colors) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.background.x = element_rect(fill = "white"),
          strip.background.y = element_rect(fill = "white")) -> gp

    gp + 
      geom_line(data = cum_means, color = "red") + 
      geom_line(data = cum_means, aes(y = upper), linetype = "dotted", color = "red", linewidth = 0.7) +
      geom_line(data = cum_means, aes(y = lower), linetype = "dotted", color = "red", linewidth = 0.7) -> gp

  
  if(var == "log_fold_change") {
    gp <- gp + 
      scale_y_continuous(limits = c(ymin,ymax),
                         name = paste0("Fold change from ", gsub("D614G", "614D/614G", comp_antigen)),
                         breaks = seq(ymin, ymax-1, by = 2),
                         labels = function(x) ifelse(x < 0, paste0("-", 2^abs(x), "x"), paste0(2^x, "x")))
    
  } else {
    
    gp <- gp + 
      scale_y_continuous(limits = c(ymin,ymax),
                         name = "Titer",
                         breaks = seq(ymin, ymax-1, by = 2),
                         labels = function(x) 2^x*10)
    
  }
  
  return(gp)
  
}


plot_over_time_lowess <- function(forest_data_sub, comp_antigen, target_sr_groups, ymin = -8, ymax = 10, stepwise = TRUE,
                                         var = "log_fold_change", span = 0.75, cum_means = NULL, omicron_variant_levels = c("BA.1", "BA.1.1", "BA.2", "BA.3", "BA.4/5", "BA.2.12.1", "BA.2.75"),
                                  show_loess = TRUE,
                                  show_linear = FALSE,
                                  cum_means_long = NULL,
                                  titers = FALSE){
  
  forest_data_sub <- forest_data_sub %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% target_sr_groups) %>%
    filter(OmicronVariant %in% omicron_variant_levels) %>%
    mutate(standardise_encounters = factor(standardise_encounters, levels = target_sr_groups))
  
  cols <- ag_colors %>%
    filter(Antigen %in% forest_data_sub$OmicronVariant)
  
  plot_colors <- cols$Color
  names(plot_colors) <- cols$Antigen
  
  mean_colors <- c("Cumulative" = "#2297E6",
                   "Total" = "tomato1")
  
  plot_colors <- c(plot_colors, mean_colors)
  # forest_data_plot$y <- rep(c(ymax - 0.5, ymax-1), nrow(forest_data_plot)/2)[1:nrow(forest_data_plot)]
  
  forest_data_sub %>%
    ggplot(aes(x = time_since_first, y = log_fold_change)) + 
    geom_point(aes(shape = shape), alpha = 0.7, color = "black", size = 1) -> gp 
  
  if(show_loess){
    gp + stat_smooth(method = "loess", se=T,degree=1, span = span, color = "#2297E6", fill = "lightblue",
                     linewidth = 0.5) -> gp
  }
  
  if(show_linear){
    
    gp + stat_smooth(method = "lm", se=T,degree=1, span = span, color = "#a48ff5", fill = "#d2b7f5",
                     linewidth = 0.5) + 
      stat_regline_equation(color = "#a48ff5")-> gp
  }  
  
  
  if(!is.null(cum_means)){
    cum_means %>% 
      filter(`Comparator antigen` == comp_antigen) %>%
      filter(standardise_encounters %in% target_sr_groups) %>%
      filter(OmicronVariant %in% omicron_variant_levels) %>%
      mutate(standardise_encounters = factor(standardise_encounters, levels = target_sr_groups),
             Average = "Cumulative") -> cum_means
    
    gp <- gp + 
      geom_ribbon(data = cum_means, aes(ymin = lower, ymax = upper, fill = Average), alpha = 0.2) +
      geom_line(data = cum_means,aes(color = Average), linewidth = 0.5) -> gp 
    
  }
  
  if(!is.null(cum_means_long)){
    cum_means_long %>%
      group_by(standardise_encounters, `Comparator antigen`, OmicronVariant) %>%
      mutate(log_fold_change = unique(log_fold_change[time_since_first == max(time_since_first, na.rm = TRUE)]),
             lower = unique(lower[time_since_first == max(time_since_first, na.rm = TRUE)]),
             upper = unique(upper[time_since_first == max(time_since_first, na.rm = TRUE)])) %>%
      ungroup() %>%
      filter(`Comparator antigen` == comp_antigen) %>%
      filter(standardise_encounters %in% target_sr_groups) %>%
      filter(OmicronVariant %in% omicron_variant_levels) %>%
      mutate(standardise_encounters = factor(standardise_encounters, levels = target_sr_groups)) %>%
      filter(time_since_first <= max(cum_means$time_since_first)) %>%
      mutate(Average = "Total")-> cum_means_long
    
    # set all values to max tp
    cum_means_long %>%
      group_by(standardise_encounters, `Comparator antigen`, OmicronVariant) %>%
      mutate(log_fold_change = log_fold_change[])
    
    gp <- gp + 
      geom_ribbon(data = cum_means_long, aes(ymin = lower, ymax = upper, fill = Average), alpha = 0.2) +
      geom_line(data = cum_means_long, aes(color = Average), linewidth = 0.3) -> gp 
    
  }
  
  gp + 
    facet_grid(factor(OmicronVariant,
                      levels = omicron_variant_levels)~standardise_encounters) + #factor(standardise_encounters) , levels = c("WT conv", "2x Vax", "3x Vax", "Inf + Vax", "Vax + Inf", "Vax + BA.1", "Vax + BA.2")) + 
    scale_x_continuous(name = "Time since first available data (months)",
                       breaks = seq(0, max(forest_data_sub$time_since_first), by = 30),
                       labels = function(x) x/30) +
    scale_shape_manual(values = c("sample" = 20, "GMT" = 21)) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) + 
    theme_bw() + 
    guides(shape = "none") + 
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          strip.background.x = element_rect(fill = "white"),
          strip.background.y = element_rect(fill = "white")) -> gp
  
  # if(!stepwise){
  #   gp + 
  #     geom_line(data = forest_data_plot) -> gp
  # }
  
  if(var == "log_fold_change" & !titers) {
    gp <- gp + 
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey20") +
      scale_y_continuous(name ="Log2(fold change) from WT",  #paste0("Log2(fold change) from ", gsub("D614G", "614D/G", comp_antigen)),
                         #    labels = function(x) ifelse(x < 0, paste0("-", 2^abs(x), "x"), paste0(2^x, "x")),
                         breaks = seq(ymin, ymax, by = 1)) + 
      coord_cartesian(ylim = c(ymin, ymax))
    
  } else {
    
    gp <- gp + 
      scale_y_continuous(name = "Titer",
                         breaks = seq(ymin, ymax, by = 2),
                         labels = function(x) round(2^x*10)) +
      coord_cartesian(ylim = c(ymin, ymax))
    
  }
  
  return(gp)
  
}
