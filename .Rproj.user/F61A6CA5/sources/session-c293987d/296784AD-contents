# netzl_et_al2024

This repository contains the code and output for the paper by Netzl et al.
Please cite the original publication if any data or code from this repository is used. 

The repository's DOI was created with Zenodo (https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)

All data and metadata can be found in the `data` directory. The raw data is a local copy of [this google sheet](https://docs.google.com/spreadsheets/d/1IvUwoWMAtJULnN-pohUoPio0RmHZIRj8OTtv9r4RFl4/edit?gid=0#gid=0) and stored in  the `data/google_sheet_tables` directory. 

The forest plots in Figure 1 and supplementary Figure is created with the `code/plotting/forest_plots.R` script.

The `code/plotting/fold_change_over_time_cumulative.R` creates Figure 2. This script also creates Table 1, stored in `data/summary_tables`.

The `code/plotting/data_since_exposure.R` creates the supplementary Figure showing titers and fold changes over time since exposure. 

The `code/cartography` directory contains the scripts to create and plot antigenic maps, antibody landscapes and perform map diagnostics. The first script to run is `create_titer_table.R`, which stores titer tables in `data/titer_tables`, followed by `map_creation.R` to create the base maps. The maps are save in `data/maps`. The comparison map  `data/maps/roessler_et_al2023.ace` was originally published in Rössler, Netzl, et al. Nat Commun 14, 5224 (2023) (https://doi.org/10.1038/s41467-023-41049-4)
and the `Wilks_et_al_map_ndsubset_no_outliers_slope_adjusted.ace` map in Wilks, Mühlemann, et al., Science (2023) (https://doi.org/10.1126/science.adj0070).

Table 2 is created in the `code/summary_tables/gmt_fc_table_updated.R` script. This directory contains other scripts to generate the supplementary tables. All tables are stored in the `data/summary_tables` directory.

In general, output figures are stored in the `figure` directory, tables in `data/summary_tables`.

All analyses were performed in R version 4.2.2 (2022-10-31).
R Core Team (2022). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna,
  Austria. URL https://www.R-project.org/.
  
Antigenic maps were constructed using the Racmacs package, Version 1.1.35:
Wilks S (2022). _Racmacs: R Antigenic Cartography Macros_. https://acorg.github.io/Racmacs,
  https://github.com/acorg/Racmacs.
  
Antibody landscapes were constructed using the ablandscapes package, Version 1.1.0: 
Wilks S (2021). _ablandscapes: Making Antibody landscapes Using R_. R package
  version 1.1.0, <https://github.com/acorg/ablandscapes>.

