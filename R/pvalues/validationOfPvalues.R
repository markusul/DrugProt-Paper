start_time <- Sys.time()

source("R/pvalues/drugDesign.R")
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
print(paste0("Analyzing time point ", tp, "h"))

# load projections
load(paste0('Z/', tp, '.RData'))

# Data for model
D <- datI[datI$pert_time == tp, pert_names]

dd <- drugDesign(D)
drug_design <- dd$design
dLabels_measured <- dd$labels
dlabels_model <- dd$labels_model
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
}

res <- lapply(sampleP_vec, function(sampleP){
  print(paste0("Analyzing protein ", prot_names[sampleP]))

  ####### simulate Y ######
  file <- paste0("results/DrugEffects/", sampleP, "_", tp, ".RData")
  load(file)

  sampledEffects <- numeric(length(dlabels_model))
  for(l in dLabels_measured) sampledEffects[which(dlabels_model == l)] <- effects.drugs[[l]]
  
  if(laggedTime > 0){
    file <- paste0("results/ProteinEffects/", sampleP, "_", tp, ".RData")
    load(file)
    sampledEffects <- c(sampledEffects, betahat)
  }
  
  # noise level estimated in the main fit (cv lasso residuals, degrees of freedom corrected)
  stopifnot(!is.null(sigmahat))
  Y_mu <- design %*% sampledEffects
  sigma <- sigmahat
  
  nOut <- length(dLabels_measured) + if(laggedTime > 0) length(prot_names) else 0
  res <- parallel::mclapply(1:nRep, mc.cores = 100, function(i){
    tryCatch({
      Y <- Y_mu + rnorm(length(Y_mu), 0, sigma)
      fit <- lasso.proj(x = design, y = Y, Z = Z, robust = FALSE, suppress.grouptesting = TRUE)
      
      # same group test as in drugInteraction.R
      pval.drugs <- sapply(dLabels_measured, function(l){
        g <- which(dlabels_model == l)
        min(1, length(g) * min(fit$pval[g]))
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
  
  true_effects <- sapply(effects.drugs, function(x) any(x != 0))
  if(laggedTime > 0){
    true_effects <- c(true_effects, betahat != 0)
  }
  
  group <- c(rep("single", 63), rep("double", 59))
  if(laggedTime > 0) group <- c(group, rep("protein", length(prot_names)))
  
  list(coef = data.frame(group = group,
                         true_effect = true_effects,
                         rate = ratioOfEffects),
       nFailed = nFailed)
})
do.call(rbind, res)
})

save(res, file = "results/PvalAnalysis.RData")


print("finished!")
end_time <- Sys.time()
print(end_time - start_time)

