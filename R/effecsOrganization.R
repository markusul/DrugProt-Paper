load('data/order.RData')
load("data/laggedData.RData")

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


res <- lapply(expTimes, function(tp){
  
D <- datI[datI$pert_time == tp, pert_names]

## prepare design matrix with interactions
drug_design <- model.matrix(~ -1 + .^2, data = D)
dLabels <- c(colnames(drug_design), colnames(drug_design))

# remove treatments without data
with_data <- apply(drug_design, 2, function(x) length(unique(x))) != 1
drug_design <- drug_design[, with_data]
dLabels_measured <- colnames(drug_design)
dlabels_model <- c(dLabels_measured, dLabels_measured)

# add intercept to each treatment
drug_intercept <- drug_design
drug_intercept[drug_design != 0] <- 1
colnames(drug_intercept) <- paste0(colnames(drug_design), "_intercept")

# combine intercept and drug design
drug_design <- cbind(drug_design, drug_intercept)


res <- sapply(1:length(prot_names), function(protein){
  load(file = paste0('results/DrugEffects/', protein , '_', tp, '.RData'))
  sapply(effects.drugs.debiased, function(dEffect){
    total_effects <- drug_design[, names(dEffect)] %*% dEffect
    meanCoef <- mean(total_effects[total_effects != 0])
    meanCoef
  })
})
res
})


res2 <- array(unlist(res), dim = c(nTreatment, length(prot_names), length(expTimes)))

for(drug in 1:nTreatment){
  dEff <- res2[drug, , ]
  save(dEff, file = paste0('results/Coef/drugs/', drug, '.RData'))
}



for(tp in expTimes[-1]){
  for(protein in 1:length(prot_names)){
    load(file = paste0('results/ProteinEffects/', protein , '_', tp, '.RData'))
    save(bhat, file = paste0('results/Coef/proteins/', protein , '_', tp, '.RData'))
  }
}


