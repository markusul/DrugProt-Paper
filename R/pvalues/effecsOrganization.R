load('data/order.RData')
load("data/laggedData.RData")
source("R/pvalues/drugDesign.R")

tp <- expTimes[2]
D <- datI[datI$pert_time == tp, pert_names]

## prepare design matrix with interactions
drug_design <- model.matrix(~ -1 + .^2, data = D)

# remove treatments without data
with_data <- apply(drug_design, 2, function(x) length(unique(x))) != 1
treatNames <- colnames(drug_design)[with_data]

load('data/drugLookup.RData')

replace_drug_ids <- function(x) {
  # 1. Split the string by backticks and underscores
  # E.g., "`drug_#123`" becomes c("", "drug", "#123", "")
  # E.g., "`drug_EV`"   becomes c("", "drug", "EV", "")
  ids <- strsplit(x, "_|`")[[1]]
  
  # 2. Filter out empty strings and the "drug" prefix
  ids <- ids[ids != "" & ids != "drug" & ids != ":"]
  
  # If there's nothing left, return the original input
  if(length(ids) == 0) return(x)
  
  # 3. Process each remaining part
  names <- sapply(ids, function(id) {
    # Check if this part is a "number" (ID). 
    # We remove the '#' for the numeric check to handle both '123' and '#123'.
    is_numeric_id <- grepl("#", id) || !is.na(as.numeric(gsub("#", "", id)))
    
    if (is_numeric_id) {
      # It's an ID: Try to get the name from drug_lookup
      val <- drug_lookup[id]
      # If it's found in the table, return the name; otherwise return the ID itself
      if (!is.na(val)) return(unname(val))
      return(id)
    } else {
      # It's already a name (like "EV", "LA", "MK"): Return it as is
      return(id)
    }
  })
  
  # 4. Join parts back with ":" (useful for drug combinations)
  paste(names, collapse = ":")
}


treatDNames <- sapply(treatNames, replace_drug_ids)
names(treatNames) <- treatDNames
save(treatNames, file = "results/Coef/drugs/treatNames.RData")


for(effect_type in c("direct", "total")){
# 6 h has no lagged proteins, so direct and total effect coincide there
res <- lapply(expTimes, function(tp){
suffix <- if(effect_type == "total" && tp != expTimes[1]) "_mode2" else ""
  
D <- datI[datI$pert_time == tp, pert_names]

dd <- drugDesign(D)
drug_design <- dd$design
dLabels_measured <- dd$labels
dlabels_model <- dd$labels_model


res <- sapply(1:length(prot_names), function(protein){
  load(file = paste0('results/DrugEffects/', protein , '_', tp, suffix, '.RData'))
  sapply(effects.drugs.debiased, function(dEffect){
    total_effects <- drug_design[, names(dEffect)] %*% dEffect
    meanCoef <- mean(total_effects[total_effects != 0])
    meanCoef
  })
})
res
})


res2 <- array(unlist(res), dim = c(nTreatment, length(prot_names), length(expTimes)))
out_dir <- if(effect_type == "direct") "results/Coef/drugs/" else "results/Coef/drugs_total/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for(drug in 1:nTreatment){
  dEff <- res2[drug, , ]
  save(dEff, file = paste0(out_dir, drug, '.RData'))
}
}

# effects of proteins
for(tp in expTimes[-1]){
  for(protein in 1:length(prot_names)){
    load(file = paste0('results/ProteinEffects/', protein , '_', tp, '.RData'))
    save(bhat, file = paste0('results/Coef/proteins/', protein , '_', tp, '.RData'))
  }
}


