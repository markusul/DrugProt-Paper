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
  
  #save model
  save(file = paste0('models/', which(prot_names == P) , '_', t, '.RData'), fit)
  
  # apply group testing for each treatments (intercept and effect)
  pMat <- matrix(NA, nrow = ncol(D), ncol = ncol(D))
  rownames(pMat) <- colnames(pMat) <- dLabels[1:ncol(D)]
  for(l in dLabels_measured){
    # use conservative = FALSE because this would correct also using the number of protein effects
    pval <- fit$groupTest(which(dlabels_model == l), conservative = TRUE)
    pval <- min(1, pval * length(dLabels_measured)) # bonferroni correction with the number of groups
    drugs <- strsplit(l, ":")[[1]]
    if(length(drugs) == 1) drugs <- c(drugs, drugs)
    pMat[drugs[1], drugs[2]] <- pval
    pMat[drugs[2], drugs[1]] <- pval
  }
  
  # collect p values for protein effects
  pval.corr <- NULL
  pval <- NULL
  
  if(laggedTime > 0){ # return p values for protein effects
    pval.corr <- fit$pval.corr
    pval.corr <- pval.corr[(length(dlabels_model)+1):length(pval.corr)]
    
    pval <- fit$pval
    pval <- pval[(length(dlabels_model)+1):length(pval)]
  }
  save(file = paste0('pvals/', which(prot_names == P) , '_', t, '.RData'), P, 
       pMat, pval, pval.corr, prot_names)
})
})

print("finished!")

