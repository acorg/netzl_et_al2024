mapValues <- function(val_fn, name_fn, map) {
  function(map) {
    values <- val_fn(map)
    names(values) <- name_fn(map)
    values
  }
}

calc_fold_changes <- function(map, titer_table) {
  # Setup to store results
  all_group_results <- tibble(
    sr_group = character(0),
    ag = character(0),
    diff = numeric(0),
    diff_upper = numeric(0),
    diff_lower = numeric(0),
    homologous = logical(0)
  )
  
  # Append results for each group
  for (sr_group in names(homologous_ags)) {
    
    sr <- which(srGroups(map) == sr_group)
    homo_ag_name <- homologous_ags[sr_group]
    homo_ag <- match(homo_ag_name, agNames(map))
    sr_group_results <- tibble(
      sr_group = rep(sr_group, numAntigens(map)),
      ag = character(numAntigens(map)),
      diff = numeric(numAntigens(map)),
      diff_upper = numeric(numAntigens(map)),
      diff_lower = numeric(numAntigens(map)),
      homologous = logical(numAntigens(map))
    )
    
    for (ag in seq_len(numAntigens(map))) {
      
      homologous_titers <- titer_table[homo_ag, sr]
      ag_titers <- titer_table[ag, sr]
      titer_diff_est <- data.frame(
        mean_diff = Rmisc::CI(na.omit(log2(as.numeric(ag_titers)) - log2(as.numeric(homologous_titers))))[["mean"]],
        mean_diff_upper = Rmisc::CI(na.omit(log2(as.numeric(ag_titers))- log2(as.numeric(homologous_titers))))[["lower"]],
        mean_diff_lower = Rmisc::CI(na.omit(log2(as.numeric(ag_titers))- log2(as.numeric(homologous_titers))))[["upper"]]
      ) 
      
      
      # Populate results
      sr_group_results$ag[ag] <- agNames(map)[ag]
      if (ag == homo_ag) {
        sr_group_results$diff[ag] <- 0
        sr_group_results$diff_upper[ag] <- 0
        sr_group_results$diff_lower[ag] <- 0
        sr_group_results$homologous[ag] <- TRUE
      } else {
        sr_group_results$diff[ag] <- titer_diff_est$mean_diff
        sr_group_results$diff_upper[ag] <- titer_diff_est$mean_diff_upper
        sr_group_results$diff_lower[ag] <- titer_diff_est$mean_diff_lower
        sr_group_results$homologous[ag] <- FALSE
      }
      
    }
    
    # Append the results
    all_group_results <- bind_rows(all_group_results, sr_group_results)
    
  }
  
  
  # Remove NA diffs
  all_group_results %>% 
    filter(
      !is.na(diff)
    ) -> all_group_results
  
  return(all_group_results)
  
}

foldchange <- \(x) {
  
  xabs <- abs(x)
  foldchange <- 2^xabs
  foldchange[x < 0] <- -foldchange[x < 0]
  as.character(round(foldchange, 1))
  
}


# plot fold change from homologous per antigen assay type
do_fold_change_plot <- function(combo, sr_group_name) {
  combo %>%
    filter(
      sr_group == sr_group_name
    ) -> sr_group_results
  
  combined_order <- sr_group_results %>% filter(`Antigen type` == "Combined")
  temp_ag_levels <- unique(combined_order$ag[order(-combined_order$diff)])
  
  sr_group_results$x <- match(sr_group_results$ag, temp_ag_levels)
  
  sr_group_results %>%
    ggplot(
      aes(
        x = x,#ifelse(`Antigen type` == "Live-virus", x - 0.15, x + 0.15),
        y = diff,
        ymin = diff_lower,
        ymax = diff_upper,
        shape = `Antigen type`,
        group = `Antigen type`
      )
    ) + 
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) + 
    geom_line(aes(color = sr_group,
                  linetype = `Antigen type`),
              position = position_dodge(width = 0.35),) +
    geom_pointrange(
      aes(color = ag),
      show.legend = TRUE,
      position = position_dodge(width = 0.35),
      size = 0.5,
    ) +
    geom_point(
      color = "white",
      show.legend = FALSE,
      position = position_dodge(width = 0.35),
      size = 0.7,
    ) +
    scale_x_continuous(
      breaks = 1:length(temp_ag_levels),
      labels = temp_ag_levels
    ) + 
    scale_color_manual(
      values = c(agFillValues(map_lv), srGroupValues(map_lv))
    ) + 
    scale_y_continuous(
      breaks = 3:min(floor(combo$diff), na.rm=T),
      labels = foldchange
    ) +
    guides(colour = "none",
           linetype = guide_legend("Antigen type")
           ) + 
    coord_cartesian(
      ylim = c(min(floor(combo$diff), na.rm=T), 3.5)
    ) +
    labs(
      x = "",
      y = "Fold change from homologous",
      title = sr_group_name
    ) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.title.y = element_text(
        size = 8
      ),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) -> gp
  
  gp <- gp + 
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Combined"),
              aes(
                x = x,
                y = 3.4,
                label = paste0("C: ",foldchange(diff))
              ),
              size = 1.7,
              color = "grey20"
    ) +
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Live-virus"),
              aes(
                x = x,
                y = 2.8,
                label = paste0("LV: ",foldchange(diff))
              ),
              size = 1.7,
              color = "grey20"
    ) +
    geom_text(data = sr_group_results %>% filter(`Antigen type` == "Pseudovirus"),
              aes(
                x = x,
                y = 2.2,
                label = paste0("PV: ",foldchange(diff))
              ),
              size = 1.7,
              color = "grey20"
    )
  
  return(gp)
}

# plot ratio difference of LV vs. PV
do_ratio_plot <- function(combo, sr_group_name) {
  
  combo %>%
    filter(
      sr_group == sr_group_name
    ) -> sr_group_results
  
  combined_order <- sr_group_results %>% filter(`Antigen type` == "Combined")
  temp_ag_levels <- unique(combined_order$ag[order(-combined_order$diff)])
  
  sr_group_results %>%
    mutate(fc = as.numeric(foldchange(diff)),
           fc = as.numeric(ifelse(fc < 0, fc, 1/fc)))-> sr_group_results
  
 # combined_order <- sr_group_results %>% filter(`Antigen type` == "Combined")
 # temp_ag_levels <- unique(combined_order$ag[order(-combined_order$diff)])
  
  sr_group_results$x <- match(sr_group_results$ag, temp_ag_levels)
  
  sr_group_results %>%
    filter(!is.null(fc)) %>%
    select(sr_group, ag, `Antigen type`, fc, x) %>%
    pivot_wider(names_from = `Antigen type`, values_from = fc) %>%
    filter(`Live-virus` != "NULL" & Pseudovirus != "NULL") %>%
    mutate(ratio_full = `Live-virus`/Pseudovirus,
           ratio = log2(abs(ratio_full))) -> ratio_df
  # ratio = ifelse(ratio_full < 0, -ratio, ratio)) -> ratio_df
 
  ratio_df %>%
    ggplot(
      aes(
        x = x,#ifelse(`Antigen type` == "Live-virus", x - 0.15, x + 0.15),
        y = ratio
      )
    ) + 
    geom_line(aes(color = sr_group)) +
    geom_point(
      aes(color = ag),
      show.legend = FALSE,
      size = 2
    ) +
    scale_x_continuous(
      breaks = 1:length(temp_ag_levels),
      labels = temp_ag_levels,
      limits = c(1, length(temp_ag_levels))
    ) + 
    scale_color_manual(
      values = c(agFillValues(map_lv), srGroupValues(map_lv))
    ) + 
    scale_y_continuous(
      labels = function(x) round(2^x,1),
      breaks = seq(-3,3,1),
      limits = c(-3.5,3.3)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) + 
    labs(
      x = "",
      y = "Fold change LV/PV"
    ) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.title.y = element_text(
        size = 8
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) -> gp_ratio
  
  return(gp_ratio)
  
}
