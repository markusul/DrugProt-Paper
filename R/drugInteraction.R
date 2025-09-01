start_time <- Sys.time()

library(hdi)
RNGkind("L'Ecuyer-CMRG")
set.seed(22)

load("data/laggedData.RData")
expTimes <- c(6, 24, 48)

#res <- lapply(prot_names, function(P){
res <- parallel::mclapply(prot_names, mc.cores = 100, function(P){
lapply(expTimes, function(t){
#parallel::mclapply(expTimes, mc.cores = 3, function(t){
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

    # differential expression to baseline
    protein_design <- protein_design - datI[datI$pert_time == t, paste0(prot_names, "_0")]
    colnames(protein_design) <- prot_names
    
    design <- cbind(design, protein_design)
    
    # remove samples without lagged protein measurements
    noLagged <- rowSums(is.na(design)) > 0
    design <- design[!noLagged, ]
    Y <- Y[!noLagged]
  }
  
  #labels
  single_effects <- rep(c(colnames(D), rep(NA, choose(ncol(D), 2))), 2)
  colnames(drug_design)
  
  # load projections
  load(paste0('Z/', t, '.RData'))
  
  #hdi fit with robustness against model misspecifications
  fit <- lasso.proj(x = design, y = Y, Z = Z, robust = FALSE)
  
  # apply group testing for each treatments (intercept and effect)
  nDrugs <- ncol(D)
  pval.drugs <- sapply(dLabels_measured, function(l){fit$groupTest(which(dlabels_model == l), conservative = FALSE)})
  save(file = paste0('results/DrugEffects/', which(prot_names == P) , '_', t, '.RData'), 
       pval.drugs, dLabels_measured, dlabels_model, nDrugs, P, t)
  
  # collect p values for protein effects
  pval <- NULL
  bhat <- NULL
  betahat <- NULL
  
  if(laggedTime > 0){ # return p values for protein effects
    pval <- fit$pval
    pval <- pval[(length(dlabels_model)+1):length(pval)]

    bhat <- fit$bhat
    bhat <- bhat[(length(dlabels_model)+1):length(bhat)]

    betahat <- fit$betahat
    betahat <- betahat[(length(dlabels_model)+1):length(betahat)]
  }
  save(file = paste0('results/ProteinEffects/', which(prot_names == P) , '_', t, '.RData'), 
       pval, prot_names, P, t, bhat, betahat)
})
})

print("finished!")
end_time <- Sys.time()
print(end_time - start_time)
