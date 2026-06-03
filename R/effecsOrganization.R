load('data/order.RData')
load("data/laggedData.RData")

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
table(drug_design[, 5])
apply(drug_design, 2, table)


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


