# helper function to load and preprocess aggregated data for AnchorForest
prepAggData <- function(){
  # Load data
  load("data/aggData.RData")
  load("data/protNames.RData")
  
  # collect drug combination perturbations
  perturbations <- unique(agg_data$pertLabel[agg_data$type == 'drugCombination'])
  dat <- agg_data

  # organize covariate names as protein names at different time points
  p_names <- c(prot_names_6, prot_names_24, prot_names_48)
  
  # select aggregated protein expressions as covariates
  X <- dat[, p_names]

  # center protein expressions by baseline (0 hr) levels
  X[, prot_names_6] <- X[, prot_names_6] - dat[, prot_names_0]
  X[, prot_names_24] <- X[, prot_names_24] - dat[, prot_names_0]
  X[, prot_names_48] <- X[, prot_names_48] - dat[, prot_names_0]
  
  # select health outcome
  Y <- dat$IC50

  # Drug concentrations as Anchor variables
  A <- dat[, pert_names]
  
  return(list(X = X, Y = Y, A = A, p_names = p_names, 
              perturbations = perturbations, 
              pertLabel = dat$pertLabel, prot_names = prot_names))
}