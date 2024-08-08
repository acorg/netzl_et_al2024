# add Omicron Variants as Comparator antigen
melt_omicron_data <- function(forest_data) {
  
  omicron_data <- forest_data %>% filter(!(is.na(TitersOmicron)))
  omicron_data$TitersHAg <- omicron_data$TitersOmicron
  omicron_data$`Comparator antigen` <- omicron_data$OmicronVariant
  omicron_data$log_fold_change <- NA
  
  
  forest_data <- rbind(forest_data, omicron_data)
  
  return(forest_data)
}

# match specific sera in data that were titrated against multiple variants
match_sera <- function(full_data, variant) {
  
  variant_sera <- forest_data %>%
    filter(!is.na(log_fold_change)) %>%
    filter(OmicronVariant %in% variant) %>%
    select(Study, `Sera details long`, OmicronVariant) %>%
    unique()
  
  data <- full_data %>%
    filter(!is.na(log_fold_change)) %>%
    filter(Study %in% variant_sera$Study) %>%
    filter(`Sera details long` %in% variant_sera$`Sera details long`)
  
  return(data)
}



mean_fc_gmt_table <- function(data, variant, group_var = "standardise_encounters", comp_antigen = "WT", match_sera = TRUE) {
  
  data$`Comparator antigen`[data$`Comparator antigen` %in% c("WT", "D614G")] <- "WT"
  # match studies by sera
  if(match_sera) {
    data <- match_sera(data, variant)
  }
  
  # calculate mean fold change and GMTs
  calc_mean_per_grouping(data, group_var)%>% filter(`Comparator antigen` %in% comp_antigen) -> means
  
  groupings <- c(group_var, "OmicronVariant")
  
  # get the number of data points
  means %>%
    group_by(across(all_of(groupings))) %>%
    mutate(n_fc = length(!(is.na(log_fold_change)))-1,
           n_hag = length(!(is.na(Log2HAg)))-1,
           n_omi = length(!(is.na(Log2Omi)))-1
                         ) -> means
  
  # go from long to wide format
  # and add column labels
  unique(means[grepl("Mean", means$rowlabel),]) %>% 
    mutate(fc_label = paste0(round(2^log_fold_change, 1),
                             " (", round(2^fc_lower, 1), "; ", round(2^fc_upper, 1), ")\n",
                             "n=" ,n_fc),
           gmt_omi_label = paste0(round(2^Log2Omi*10, 0),
                              " (", round(2^lower_omi*10, 0), "; ", round(2^higher_omi*10, 0), ")\n",
                              "n=" ,n_omi),
           gmt_hag_label = paste0(round(2^Log2HAg*10, 0),
                                  " (", round(2^lower_hag*10, 0), "; ", round(2^higher_hag*10, 0), ")\n",
                                  "n=" ,n_hag)) %>%
    select("OmicronVariant", group_var, "fc_label", "gmt_omi_label", "gmt_hag_label") %>%
    mutate("Lineage" = paste0(variant, collapse = ", "),
           "Comparator antigen" = paste0(comp_antigen, collapse = ", ")) -> mean_table
  
  mean_table <- apply(mean_table, c(1,2), function(x) gsub("\\(NA; NA\\)", "", x))

  mean_table <- as.data.frame(mean_table)
  mean_table$standardise_encounters <- factor(mean_table$standardise_encounters, levels = rev(levels(data$standardise_encounters)))

  mean_table <- mean_table[,c("Comparator antigen", group_var, "OmicronVariant", "Lineage", "fc_label", "gmt_omi_label", "gmt_hag_label")]
  fc_table <- mean_table %>% select(!gmt_omi_label:gmt_hag_label) %>%
    pivot_wider(names_from = OmicronVariant, values_from = fc_label)
  fc_table[is.na(fc_table)] <- "-"
  
  titer_table <- mean_table %>% select(!fc_label) %>%
    pivot_wider(names_from = OmicronVariant, values_from = c("gmt_hag_label", "gmt_omi_label"))
  titer_table[is.na(titer_table)] <- "-"
                
  
  table_list <- list("full_table" = mean_table,
                     "fc_table" = fc_table,
                     "titer_table" = titer_table
                     )
  return(table_list)
}

get_fc_table <- function(data, variant, group_var = "standardise_encounters", comp_antigen = "WT", match_sera = TRUE) {
  mean_fc_gmt_table(data, variant, group_var, comp_antigen, match_sera)$fc_table
}

get_titer_table <- function(data, variant, group_var = "standardise_encounters", comp_antigen = "WT", match_sera = TRUE) {
  mean_fc_gmt_table(data, variant, group_var, comp_antigen, match_sera)$titer_table
}

wider_fc_variant_tables <- function(data, variants = c("BA.1", "BA.2", "BA.4/5", "BA.1.1", "BA.2.12.1", "BA.3"),
                                    group_var = "standardise_encounters", comp_antigen = "WT", match_sera = TRUE) {
  
  data <- data %>%
    filter(OmicronVariant %in% variants) %>%
    filter(`Comparator antigen` %in% comp_antigen)
  
  matched_fc <- get_fc_table(data, variants[1], group_var, comp_antigen, match_sera)
  for(var in variants[2:length(variants)]) {
    matched_fc <- matched_fc %>%
      rbind(., get_fc_table(data, var, group_var, comp_antigen, match_sera))
  }
  
  matched_fc <- matched_fc[order(matched_fc %>% pull(group_var)),]
  colnames(matched_fc) <- gsub("standardise_encounters", "Serum group", colnames(matched_fc))
  
  return(matched_fc)
    
}

wider_titer_variant_tables <- function(data, variants = c("BA.1", "BA.2", "BA.4/5", "BA.1.1", "BA.2.12.1", "BA.3"),
                                       group_var = "standardise_encounters", comp_antigen = "WT", match_sera = TRUE) {
  
  data <- data %>%
    filter(OmicronVariant %in% variants) %>%
    filter(`Comparator antigen` %in% comp_antigen)
  
  matched_titer <- get_titer_table(data, variants[1], group_var, comp_antigen, match_sera) 
  
  for(var in variants[2:length(variants)]) {
    matched_titer <- matched_titer %>%
      rbind(., get_titer_table(data, var, group_var, comp_antigen, match_sera))
  }
  
  matched_titer <- matched_titer[order(matched_titer %>% pull(group_var)),]
  matched_titer[comp_antigen] <- matched_titer[,grep("gmt_hag", colnames(matched_titer))[[1]]]
  colnames(matched_titer) <- gsub("standardise_encounters", "Serum group", colnames(matched_titer))
  
  matched_titer <- matched_titer[,!grepl("gmt_hag", colnames(matched_titer))]
  
  colnames(matched_titer) <- gsub("gmt_omi_label_", "", colnames(matched_titer))
  
  return(matched_titer)
}

huxtable_row_borders <- function(ht, column_name, base){
  
  start <- base+1
  
  # set row format
  col_number <- grep(column_name, colnames(ht))[[1]]
  
  nr_rows <- ht[start:nrow(ht),] %>% 
    group_by(across(all_of(column_name))) %>% count() %>%
    column_to_rownames(column_name)
  
  order <- unique(ht[start:nrow(ht),] %>% pull(column_name))
  
  nr_rows <- nr_rows[order,]
  
  nr_rows <- nr_rows[nr_rows > 0]
  
  nr_rows <- c(0, nr_rows)
  
  for(i in 2:(length(nr_rows))) {
    start_pos <-start+sum(nr_rows[(i-1):1])
    
    end_pos <-(start-1)+sum(nr_rows[i:1])
    
    ht <- ht %>%
      merge_cells(.,row = start_pos:end_pos, col = col_number) %>%
      set_bottom_border(row =end_pos, col = col_number:ncol(ht), value = 0.2)
    
  }
  
  return(ht)
}

format_huxtable_variant_table <- function(data, target_variants, column_names = c("Comparator antigen", "Serum group"), font_size = 6) {
  
  base <- 1

  hux_t <- hux(data) 
  
  if(TRUE %in% grepl("gmt_", colnames(data))) {
    print("in here")
    base <- base + 1
    
    col_names <- colnames(data) %>%
      gsub("gmt_hag_label_", "", .) %>%
      gsub("gmt_omi_label_", "", .)
    
    hux_t <- hux_t %>%
      insert_row(col_names, after = 1)
    
    hux_t[1,grepl("hag", colnames(hux_t))] <- "GMT Comparator antigen"
    hux_t[1,grepl("omi", colnames(hux_t))] <- "GMT Omicron"
    
    hux_t[1, colnames(hux_t) == col_names] <- ""
    hux_t <- hux_t %>%
      merge_cells(row = 1, col = grepl("hag", colnames(hux_t))) %>%
      merge_cells(row = 1, col = grepl("omi", colnames(hux_t)))
    }
  
  
  hux_t <- hux_t %>%
    set_align(everywhere, everywhere, "center") %>%
    set_all_padding(1) %>% 
    set_outer_padding(0.1) %>% 
    set_bold(row = 1, col = everywhere) %>% 
    set_bold(row = everywhere, col = 1:length(column_names)) %>% 
    set_bottom_border(row = base, col = everywhere) %>% 
    set_width(1) %>% 
    set_font_size(font_size) %>% 
    set_font_size(row = 1:base, col = everywhere, value = font_size + 1) %>% 
    theme_article()
  
  non_variant_cols <- length(column_names) + 1
  
  # make cells bold
  for(n in (base+1):nrow(hux_t)) {
    for(c in non_variant_cols:ncol(hux_t)) {
     
      if(hux_t$Lineage[n] == hux_t[base,c]) {
        hux_t <- hux_t %>% 
          set_bold(row = n, col = c)
      } 
    }
  }
  
  for(col in column_names) {
    hux_t <- huxtable_row_borders(hux_t, col, base)
  }

  hux_t <- hux_t %>%
    set_right_border(row = everywhere, col = colnames(hux_t)[!grepl(paste(target_variants, collapse = "|"), colnames(hux_t))], value = 0.2) %>%
    set_outer_borders(row = everywhere, col = everywhere)
  
  if(TRUE %in% grepl("gmt_", colnames(data))) {
    hux_t <- hux_t %>%
      set_right_border(row = everywhere, col = (ncol(hux_t) - length(target_variants)), value = 0.2) %>%
      set_width(1.2)
  }
  
  return(hux_t)
}

save_summary_table <- function(table, file_name, table_dir, sheet_name){
  
  if(TRUE %in% grepl("titer", file_name)) {
    table <- table[,2:ncol(table)] %>%
      set_outer_borders(row = everywhere, col = everywhere)
  }
  
  ft <- as_flextable(table)
  
  my_doc <- officer::read_docx()
  my_doc <- flextable::body_add_flextable(
      my_doc, ft)
  print(my_doc, target =
            file.path(table_dir, sheet_name, file_name))
  
}


mean_ecounter_table_per_ag <- function(mean_data_encounter, target_ag, word_doc, table_counter = 1, only_omicron = FALSE){
  
  comp_ags <- NULL
  if(target_ag == "BA.1"){
    mean_data_encounter %>%  select(!gmt_Omic) %>% 
      filter(OmicronVariant == target_ag) %>%
      pivot_wider(names_from = c("Comparator antigen"), values_from = c("gmt_hAG", "mean_fold_drop", "mean_fold_drop_to_WT")) %>%
      select(!c("OmicronVariant", "mean_fold_drop_to_WT_WT", paste0("mean_fold_drop_", target_ag), paste0("mean_fold_drop_to_WT_", target_ag)))-> mean_encounters_wide
    
  } else {
    
    mean_data_encounter %>% 
      filter(OmicronVariant == target_ag) %>%
      select(!mean_fold_drop_to_WT) -> mean_data_encounter
    
    comp_ags <- unique(mean_data_encounter$`Comparator antigen`)
    comp_ags <- comp_ags[comp_ags %in%c("WT", "Alpha", "Beta", "Gamma", "Delta")]
    
    mean_data_encounter %>% 
      pivot_wider(names_from = c("Comparator antigen"), values_from = c("gmt_hAG","gmt_Omic", "mean_fold_drop")) %>%
      select(c("standardise_encounters", paste0("gmt_Omic_", target_ag), paste0("mean_fold_drop_", comp_ags))) -> mean_encounters_wide
    
  }
  
  mean_encounters_wide <- mean_encounters_wide[rev(order(mean_encounters_wide$standardise_encounters)),]
  mean_encounters_wide <- clear_NA_from_table(mean_encounters_wide)
  
  # create pretty table
  ht <- hux(mean_encounters_wide) %>%
    insert_row(colnames(mean_encounters_wide))
  
  ht <- format_huxtable_cols_standard(ht, target_ag = target_ag, only_omicron = only_omicron, comp_ags = comp_ags)
  
  ft <- as_flextable(ht)
  ft <- flextable::set_table_properties(
    ft,
    layout = "autofit",
    width = 1,
    align = "left"
  )
  
  if(table_counter > 1){
    word_doc <- word_doc %>%
      officer::body_add_break()
  }
  word_doc <- word_doc %>%
    officer::body_add_par(paste0("Supplementary Table ", table_counter, ": ",target_ag, " Geometric Mean Titer (GMT) and mean fold drops per serum group (conv = convalescent). Mean fold drops were calculated from fold drops per study, not from GMT fold drops. Studies reporting fold drops but not GMTs result in a discrepancy between GMT based mean fold drop and individual study based mean fold drop. The 95%CI is given in parentheses, the number of data points n in the next line. WT summarizes wild-type like strains (e.g: 614D, 614G)."))
  word_doc <- flextable::body_add_flextable(
    word_doc, ft)
  
  return(word_doc)
}

format_stat_table_flex <- function(tab){
  
  ht <- hux(tab)
  
  ht %>% set_align(1, everywhere, "center") %>%
    set_outer_padding(0.1) %>% 
    set_width(0.9) %>% 
    set_font_size(7) %>%
    theme_article() -> ht
  
  ft <- as_flextable(ht)
  ft <- flextable::set_table_properties(
    ft,
    layout = "autofit",
    width = 1,
    align = "left"
  ) 
  
  return(ft)
}

# function to remove NA values from table
clear_NA_from_table <- function(table) {
  levels_encounter <- levels(table$standardise_encounters)
  table <- as.data.frame(table)
  for(col in colnames(table)) {
    table[,col] <- gsub("NaN", "NA", table[,col])
    table[,col] <- gsub("NA\\\n\\(NA\\; NA\\)", "", table[,col])
    table[,col] <- gsub("\\(NA\\; NA\\)", "", table[,col])
    table[,col] <- gsub("NA \\(n=1\\)", "", table[,col])
    
  }
  
  table$standardise_encounters <- factor(table$standardise_encounters, levels = levels_encounter)
  return(table)
}

add_ba.1.1_stats <- function(table1, table2) {
  
  table1_rows <- match(interaction(table2$standardise_encounters, table2$OmicronVariant),
                       interaction(table1$standardise_encounters, table1$OmicronVariant))
  
  table2_rows <- 1:length(table1_rows)
  
  for(col in c("mean_fold_drop", "gmt_hAG", "gmt_Omic"))
    for(x in 1:length(table1_rows)) {
      table1[[table1_rows[x], col]] <- paste0(table1[[table1_rows[x], col]], "\n\n", table2[[table2_rows[x], col]])
    }
  
  return(table1)
}

