args = commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

source("R/pvalues/drugDesign.R")
library(hdi)
RNGkind("L'Ecuyer-CMRG")
set.seed(22)

load("data/laggedData.RData")
expTimes <- c(6, 24, 48)

# protein of interest
P <- prot_names[1]
tp <- expTimes[as.numeric(args[1])]

print(P)
print(tp)

# Data for model
Y_full <- datI[datI$pert_time == tp, P] - datI[datI$pert_time == tp, paste0(P, "_0")]
D <- datI[datI$pert_time == tp, pert_names]

dd <- drugDesign(D)
drug_design <- dd$design
dLabels_measured <- dd$labels
dlabels_model <- dd$labels_model
design <- drug_design

#protein design
laggedTime <- which(expTimes == tp) - 1
Y <- Y_full

if(laggedTime > 0){
  protein_design <- aggData[[laggedTime]][datI[datI$pert_time == tp, 'label'], ]

  # differential expression to baseline
  protein_design <- protein_design - datI[datI$pert_time == tp, paste0(prot_names, "_0")]
  colnames(protein_design) <- prot_names
  design <- cbind(design, protein_design)
  
  # remove samples without lagged protein measurements
  noLagged <- rowSums(is.na(design)) > 0
  design <- design[!noLagged, ]
  Y <- Y[!noLagged]
}

#labels
single_effects <- rep(c(colnames(D), rep(NA, choose(ncol(D), 2))), 2)

if(laggedTime > 0){
  #hdi fit
  fit <- lasso.proj(x = drug_design, y = Y_full, return.Z = T, suppress.grouptesting = T)
  #save Z
  Z <- fit$Z
  print(dim(Z))
  save(Z, file = paste0('Z/', tp, '_mode2.RData'))
  print("mode 2 saved")
}

#hdi fit
fit <- lasso.proj(x = design, y = Y, return.Z = T, suppress.grouptesting = T)
#save Z
Z <- fit$Z
print(dim(Z))
save(Z, file = paste0('Z/', tp, '.RData'))


print("finished!")
end_time <- Sys.time()
print(end_time - start_time)
