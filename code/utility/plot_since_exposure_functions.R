ag_colors <- read.csv(file.path("data", "metadata", "ag_colors.csv"))


fold_change_since_exposure <- function(data, variants = c("BA.1"), comp_antigen = c("D614G"), target_sr_groups = c("2x Vax", "3x Vax"), 
                                       ymin = -8, ymax = 1) {
  
  temp <- data %>%
    filter(OmicronVariant %in% variants) %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% target_sr_groups) %>%
    mutate(standardise_encounters = factor(standardise_encounters, levels = target_sr_groups))
  
  ggplot(data = temp, aes_string(x="mean_days", y = "log_fold_change")) +
    geom_point(alpha = 0.7, color = "black", size = 1) + 
    # geom_smooth(method='lm', color = "#a48ff5", fill = "#d2b7f5",
    #             linewidth = 0.5) + 
    # stat_regline_equation(color = "#a48ff5",
    #                       label.y.npc = "top")+
    facet_grid(OmicronVariant ~ standardise_encounters) +
    #  scale_y_continuous(name = names[x], labels = function(y) 2^y, breaks = c(0:10), limits = c(0,10)) +
    scale_x_continuous(name = "Time since exposure (months)", labels = function(d) d/30, breaks = seq(from = 0, to =360, by = 30), limits = c(0, 360)) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey20") +
    scale_y_continuous(name = paste0("Log2 fold change from ", gsub("D614G", "614D/G", comp_antigen)),
                      # labels = function(x) ifelse(x < 0, paste0("-", 2^abs(x), "x"), paste0(2^x, "x")),
                       breaks = seq(ymin, ymax, by = 1)) + 
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_bw() +
    theme(strip.background =element_rect(fill="white")) -> p
  
  return(p)
  
}

titer_since_exposure <- function(data, variants = c("BA.1"), comp_antigen = c("D614G"), target_sr_groups = c("2x Vax", "3x Vax"), 
                                       ymin = -8, ymax = 1) {
  
  temp <- data %>%
    filter(!is.na(Log2Omi)) %>%
    filter(OmicronVariant %in% variants) %>%
    filter(`Comparator antigen` %in% comp_antigen) %>%
    filter(standardise_encounters %in% target_sr_groups) %>%
    mutate(standardise_encounters = factor(standardise_encounters, levels = target_sr_groups))
  
  ggplot(data = temp, aes_string(x="mean_days", y = "Log2Omi")) +
    geom_point(alpha = 0.7, color = "black", size = 1) + 
    # geom_smooth(method='lm', color = "#a48ff5", fill = "#d2b7f5",
    #             linewidth = 0.5) + 
    # stat_regline_equation(color = "#a48ff5",
    #                       label.y.npc = "top") +
    facet_grid(OmicronVariant ~ standardise_encounters) +
    #  scale_y_continuous(name = names[x], labels = function(y) 2^y, breaks = c(0:10), limits = c(0,10)) +
    scale_x_continuous(name = "Time since exposure (months)", labels = function(d) d/30, breaks = seq(from = 0, to =360, by = 30), limits = c(0, 360)) + 
    scale_y_continuous(name = "Log2(GMT/10)",
                       breaks = seq(ymin, ymax, by = 2)) + 
    coord_cartesian(ylim = c(ymin, ymax)) +
    theme_bw() +
    theme(strip.background =element_rect(fill="white")) -> p
  
  return(p)
  
}
