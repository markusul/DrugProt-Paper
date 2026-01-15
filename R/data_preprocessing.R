print("start data preprocessing")

library(dplyr)
library(fastDummies)

# Impute missing values with 80% of minimum
# as we assume missing values are below detection limit
na_imputation <- function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) * 0.8
  x
}

# Load data
dat <- read.csv("data/ProteinMatrix_sampleID_MapEC50_20240229.csv")

# Create drug lookup
# What drug number corresponds to which drug name
drugnames <- unique(dat[, c('pert_iname', 'pert_id')])
drug_lookup <- setNames(drugnames$pert_iname, drugnames$pert_id)
save(drug_lookup, file = "data/drugLookup.RData")

# ---- Protein data preprocessing ----
p_idx <- 2:5586
data_protein <- dat[, p_idx]

# Remove constant columns and keep only HUMAN proteins
data_protein <- data_protein[, apply(data_protein, 2, function(x) length(unique(x))) > 1]
data_protein <- data_protein[, grepl("HUMAN", names(data_protein))]

# Analyze missing values
na_count <- colMeans(is.na(data_protein))
save(na_count, file = "data/na_count.RData")

# Impute and log transform
data_protein <- apply(data_protein, 2, na_imputation)
data_protein <- log(data_protein)

# ---- Drug data preprocessing ----
drugs <- dat$pert_id
data_drugs <- dummy_cols(data.frame(drug = drugs), remove_first_dummy = FALSE)[, -1]
data_drugs <- data_drugs * 10

combination_idx <- which(drugs == "")

# Update drug combinations with actual doses
for (i in combination_idx) {
  drugAB <- paste0("drug_", strsplit(dat[i, 'drugIdAB'], " ")[[1]])
  if (all(drugAB %in% colnames(data_drugs))) {
    data_drugs[i, drugAB[1]] <- dat[i, 'Anchor_dose']
    data_drugs[i, drugAB[2]] <- dat[i, 'Library_dose']
  }
}

# Remove reference columns
data_drugs <- data_drugs[, !colnames(data_drugs) %in% c("drug_", "drug_no")]

# ---- Additional metadata ----
data_additional <- dat[, c('pert_time', 'protein_plate', 'machine', 'BioRep', 'Sample_ID', 'Anchor_dose', 'Library_dose')]
data_additional$Anchor_dose[is.na(data_additional$Anchor_dose)] <- 0
data_additional$Library_dose[is.na(data_additional$Library_dose)] <- 0

# ---- Response variable ----
data_response <- dat$EC50
data_response[combination_idx] <- dat$Combo.IC50[combination_idx]

# ---- Combine all data ----
type <- rep('singleDrug', nrow(data_protein))
type[combination_idx] <- 'drugCombination'
type[which(drugs == 'no')] <- 'noDrug'

pertLabel <- dat$pert_id
pertLabel[combination_idx] <- dat$drugIdAB[combination_idx]
data_response[drugs == 'no'] <- Inf

data <- cbind(data_protein, data_drugs, data_additional, type, pertLabel, IC50 = data_response, NY = dat$NY)
data[data$type == 'singleDrug', 'IC50'] <- log2(data[data$type == 'singleDrug', 'IC50'])

save(data, file = 'data/prepData.RData')