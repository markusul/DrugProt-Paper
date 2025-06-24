library(hdi)
set.seed(22)

load("data/laggedData.RData")
expTimes <- c(48)

# protein of interest
P <- prot_names[1]

for(t in expTimes){
  print(P)
  print(t)
  
  # Data for model
  Y <- datI[datI$pert_time == t, P] - datI[datI$pert_time == t, paste0(P, "_0")]
  D <- datI[datI$pert_time == t, pert_names]
  
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
  design <- drug_design
  
  #protein design
  laggedTime <- which(expTimes == t) - 1
  
  if(laggedTime > 0){
    max(datI[datI$pert_time == t, 'label'])
    sum(is.na(aggData[[laggedTime]]))
    sum(rowSums(is.na(aggData[[laggedTime]])) > 0)
    protein_design <- aggData[[laggedTime]][datI[datI$pert_time == t, 'label'], ]
    
    design <- cbind(design, protein_design)
    
    # remove samples without lagged protein measurements
    noLagged <- rowSums(is.na(design)) > 0
    design <- design[!noLagged, ]
    Y <- Y[!noLagged]
  }
  
  #labels
  single_effects <- rep(c(colnames(D), rep(NA, choose(ncol(D), 2))), 2)
  colnames(drug_design)
  
  #hdi fit
  fit <- lasso.proj(x = design, y = Y, return.Z = T, 
                    suppress.grouptesting = T, 
                    parallel = T, ncores = 100)
  #save Z
  Z <- fit$Z
  save(file = paste0('Z/', t, '.RData'), Z)
}

print("finished!")