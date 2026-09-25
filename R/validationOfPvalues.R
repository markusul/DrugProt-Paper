start_time <- Sys.time()

library(hdi)
RNGkind("L'Ecuyer-CMRG")
set.seed(22)

load("data/laggedData.RData")
expTimes <- c(6, 24, 48)

# confidence level to analyze
alpha <- 0.05

#Number of true effects sampled
nEffects <- 100
# take the estimated effects of a random protein
sampleP_vec <- sample(1:length(prot_names), nEffects, replace = F)

#number of repeated drawings for p value validation
nRep <- 200

res <- lapply(expTimes, function(tp){

# load projections
load(paste0('Z/', tp, '.RData'))

# Data for model
D <- datI[datI$pert_time == tp, pert_names]

## prepare design matrix with interactions
drug_design <- model.matrix(~ -1 + .^2, data = D)

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
laggedTime <- which(expTimes == tp) - 1

if(laggedTime > 0){
  protein_design <- aggData[[laggedTime]][datI[datI$pert_time == tp, 'label'], ]
  
  # differential expression to baseline
  protein_design <- protein_design - datI[datI$pert_time == tp, paste0(prot_names, "_0")]
  colnames(protein_design) <- prot_names
  design <- cbind(design, as.matrix(protein_design))
  
  # remove samples without lagged protein measurements
  noLagged <- rowSums(is.na(design)) > 0
  design <- design[!noLagged, ]
  datI <- datI[datI$pert_time == tp, ][!noLagged, ]
}

res <- lapply(sampleP_vec, function(sampleP){
  ####### simulate Y ######
  file <- paste0("results/DrugEffects/", sampleP, "_", tp, ".RData")
  load(file)
  drugEffects <- do.call(rbind, effects.drugs)
  
  sampledEffects <- c(drugEffects[, 1], drugEffects[, 2])
  
  if(laggedTime > 0){
    file <- paste0("results/ProteinEffects/", sampleP, "_", tp, ".RData")
    load(file)
    sampledEffects <- c(sampledEffects, betahat)
  }
  
  Y_true <- datI[datI$pert_time == tp, prot_names[sampleP]] - datI[datI$pert_time == tp, paste0(prot_names[sampleP], "_0")]
  
  Y_mu <- design %*% sampledEffects
  sigma <- sd(Y_true - Y_mu)
  
  nOut <- length(dLabels_measured) + if(laggedTime > 0) length(prot_names) else 0
  res <- parallel::mclapply(1:nRep, mc.cores = 100, function(i){
    tryCatch({
      Y <- Y_mu + rnorm(length(Y_mu), 0, sigma)
      fit <- lasso.proj(x = design, y = Y, Z = Z, robust = FALSE)
      
      pval.drugs <- sapply(dLabels_measured, function(l){
        fit$groupTest(which(dlabels_model == l), conservative = FALSE)
      })
      estim_effects <- pval.drugs < alpha
      
      if(laggedTime > 0){
        pval <- fit$pval[(length(dlabels_model) + 1):length(fit$pval)]
        estim_effects <- c(estim_effects, pval < alpha)
      }
      estim_effects
    }, error = function(e){
      message("tp ", tp, ", protein ", sampleP, ", rep ", i, ": ", conditionMessage(e))
      rep(NA, nOut)
    })
  })
  
  # safety net: anything that still isn't a proper vector (e.g. a killed fork)
  bad <- vapply(res, function(r) inherits(r, "try-error") || length(r) != nOut, logical(1))
  res[bad] <- list(rep(NA, nOut))
  
  resMat <- do.call(cbind, res)
  nFailed <- sum(colSums(is.na(resMat)) == nOut)
  ratioOfEffects <- rowMeans(resMat, na.rm = TRUE)
  
  true_effects <- apply(drugEffects, 1, function(x) any(x != 0))
  if(laggedTime > 0){
    true_effects <- c(true_effects, betahat != 0)
  }
  
  group <- c(rep("single", 63), rep("double", 59))
  if(laggedTime > 0) group <- c(group, rep("protein", length(prot_names)))
  
  type1 <- tapply(ratioOfEffects[!true_effects], group[!true_effects], mean)
  power <- tapply(ratioOfEffects[true_effects], group[true_effects], mean)
  list(type1 = type1, power = power, nFailed = nFailed)
})
do.call(rbind, res)
})

save(res, file = "results/PvalAnalysis.RData")


print("finished!")
end_time <- Sys.time()
print(end_time - start_time)

