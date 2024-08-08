# preparation for forest plots, tables
format_table_forest_plots <- function(forest_data, remove_4xVax_ba45bt_samples = TRUE) {
  colnames(forest_data) <- gsub("\\.", " ", colnames(forest_data))
  
  forest_data$Sourcelink <- gsub(" ", "", forest_data$Sourcelink)
  
  forest_data <- forest_data %>%
    mutate(log_fold_change = conv_to_logfold(`numerical Titre drop`)) 
  
  forest_data <- as.data.frame(forest_data)
  # ----------------------------------------------- prepare data -----------------------------------------------
  forest_data <- clean_up_data(forest_data)
  
  # factorise variables and check for NAs
  forest_data <- factorise(forest_data)

  # get proper rowlabel with monospacing
  time_forest_data <- standardise_time_to_mean_days(forest_data)
  time_forest_data <- character_rowlabel_standard_time(time_forest_data)
  
  # add time rowlabel
  forest_data <- character_rowlabel(forest_data)
  forest_data$rowlabel <- time_forest_data$rowlabel
  
  
  forest_data$arrow_length <- unlist(lapply(forest_data$`Titre drop`, function(x) {
    if(length(grep(">>|large", x))>0) {
      1
    } else if(length(grep(">", x))>0) {
      0.5
    } else {
      NA
    }
  }))
  
  # add log2 titers and log2 omicron
  forest_data$Log2HAg[is.na(forest_data$Log2HAg)] <- log2(forest_data$TitersHAg[is.na(forest_data$Log2HAg)]/10)
  
  # get omicron titer by subtracting log fold change from log2 Titer of comparator antigen
  forest_data$Log2Omi <- forest_data$Log2HAg - forest_data$log_fold_change
  forest_data$TitersOmicron[is.na(forest_data$TitersOmicron)] <- 2^(forest_data$Log2Omi[is.na(forest_data$TitersOmicron)])*10
  
  forest_data <- pretty_plot_names(forest_data)
  
  # remove files without source link
  forest_data <- forest_data[grepl("http|DOI|doi", forest_data$Sourcelink), ]
  
  # remove sr groups with just 1 sample
  if(remove_4xVax_ba45bt_samples){
    forest_data %>%
      filter(!(standardise_encounters %in% c("Vax + BA.4/5", "4x Vax", "BA.4/5 conv"))) -> forest_data
  }
  # add variant data
  # melt variant titer data into comparator antigen
  forest_data <- melt_variant_comparator(forest_data, "Alpha")
  forest_data <- melt_variant_comparator(forest_data, "Beta")
  forest_data <- melt_variant_comparator(forest_data, "Gamma")
  forest_data <- melt_variant_comparator(forest_data, "Delta")
  
  # set titers below 10 to 5 (log2 of below 0 to -1)
  forest_data$Log2HAg[forest_data$Log2HAg < 0] <- -1
  forest_data$Log2Omi[forest_data$Log2Omi < 0] <- -1
  
  # add row long for plotting purposes
  forest_data["Row_long"] <- 1:nrow(forest_data)
  
  # remove single nucleotide subst omicron variants
  forest_data <- forest_data %>%
    filter(OmicronVariant %in% c("BA.1","BA.1.1", "BA.2","BA.2.12.1", "BA.3", "BA.4/5", "BA.2.75"))
  
  
  return(forest_data)
}