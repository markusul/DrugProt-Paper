# Load required packages
library(dplyr)
library(tidyr)

# Load prepared data
load("data/prepData.RData")
agg_input <- data   # keep original under a clearer name

# Number of protein measurement columns (first n_protein columns of raw_data)
n_protein <- 5519
names(agg_input)[n_protein + 0:2]
names(agg_input)[1]

prot_names <- names(agg_input)[1:n_protein]

# Determine additional perturbation-related columns that follow the protein columns
n_pert_single <- length(unique(agg_input$pertLabel[agg_input$type == 'singleDrug']))
pert_names <- names(agg_input)[(1 + n_protein):(n_protein + n_pert_single)]

# treat all no drug experiments equally (singleDrug/drugCombination experiments)
agg_input$pertLabel[agg_input$pertLabel == 'no no'] <- 'no'

# Drop columns that are not needed for aggregation of measurements
drop_cols <- c('type', 'Sample_ID', 'BioRep', 'machine', 'NY')
agg_input <- agg_input %>% select(-any_of(drop_cols))

# Aggregate measurements to the mean per group (protein_plate, pertLabel, Anchor_dose, Library_dose, pert_time)
# Aggregate protein columns, perturbation indicator columns and IC50 by mean (na.rm = TRUE)
agg_data <- agg_input %>%
    group_by(protein_plate, pertLabel, Anchor_dose, Library_dose, pert_time) %>%
    summarise(across(all_of(c(prot_names, pert_names, "IC50")), mean, na.rm = TRUE),
                        .groups = "drop")

# Extract baseline protein measurements (pertLabel == 'no')
baseline_prot <- agg_data %>% filter(pertLabel == 'no')

# Keep only perturbation rows (exclude baseline rows) for constructing perturbation records
agg_pert <- agg_data %>% filter(pertLabel != 'no')

# Use only cell lines that have baseline data
valid_plates <- intersect(unique(agg_pert$protein_plate), unique(baseline_prot$protein_plate))
agg_pert <- agg_pert %>% filter(protein_plate %in% valid_plates)

# Prepare baseline rows to attach to each perturbation record for the same protein_plate
pert_info_cols <- c('protein_plate', 'pertLabel', 'Anchor_dose', 'Library_dose', 'IC50', pert_names)
baseline_rows <- agg_pert %>%
    select(all_of(pert_info_cols)) %>%
    distinct() %>%
    mutate(pert_time = 0)

# Fill protein measurements in baseline_rows from baseline_prot by protein_plate
# For each protein_plate in baseline_rows, copy the prot_names values from baseline_prot
baseline_rows <- baseline_rows %>%
    left_join(baseline_prot %>% select(protein_plate, all_of(prot_names)),
                        by = "protein_plate")

# Combine perturbation aggregated data with the baseline rows (so each perturbation has a time 0 baseline)
combined <- bind_rows(agg_pert, baseline_rows)

# Pivot protein measurements wide so that each pert_time becomes its own set of protein columns
# Values are the protein columns; non-value columns (grouping/meta) are kept as is.
combined_wide <- combined %>%
    pivot_wider(names_from = pert_time, values_from = all_of(prot_names),
                            names_sep = "_")  # add suffix so column names are unique

# Remove rows that have any missing protein measurements across the time points
combined_wide <- combined_wide %>% drop_na()

# Add a type column: 'singleDrug' by default, 'drugCombination' when Anchor_dose != 0
agg_data <- combined_wide %>%
    mutate(type = if_else(Anchor_dose != 0, "drugCombination", "singleDrug"))

# Save the processed aggregated data
save(agg_data, file = "data/aggData.RData")