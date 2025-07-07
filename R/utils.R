prepAggData <- function(path = "data/aggData.RData"){
  load(path)
  load("data/protNames.RData")
  
  perturbations <- unique(agg_data$pertLabel[agg_data$type == 'drugCombination'])
  
  dat <- agg_data
  p_names <- c(prot_names_6, prot_names_24, prot_names_48)
  
  X <- dat[, p_names]
  X[, prot_names_6] <- X[, prot_names_6] - dat[, prot_names_0]
  X[, prot_names_24] <- X[, prot_names_24] - dat[, prot_names_0]
  X[, prot_names_48] <- X[, prot_names_48] - dat[, prot_names_0]
  
  Y <- dat$IC50
  
  A <- dat[, pert_names]
  
  return(list(X = X, Y = Y, A = A, p_names = p_names, perturbations = perturbations, pertLabel = dat$pertLabel, prot_names = prot_names))
}