library(SDModels)

args = commandArgs(trailingOnly = TRUE)
set.seed(42)
start_time <- Sys.time()

# helper to load and preprocess data
source("R/utils.R")

# load and preprocess data
dat <- prepAggData()
X <- as.matrix(dat$X)
Y <- dat$Y
A <- as.matrix(dat$A)
p_names <- dat$p_names
perturbations <- dat$perturbations
pertLabel <- dat$pertLabel

# only use features with enough variability
X <- X[, apply(X, 2, function(x)length(unique(x))) >= 200]
n <- length(Y)

# define 8 drug combination groups based on equal anchor drug
anchor_dose <- sapply(strsplit(perturbations, ' '), function(x) x[1])
Groups <- unique(anchor_dose)

envs <- rep(0, n)
for(group in Groups){
  test_pert <- perturbations[anchor_dose == group]
  envs[which(pertLabel %in% test_pert)] <- group
}
envs <- as.factor(envs)

# fit 200 trees for each environment left out
trees_envs <- rep(200, length(levels(envs)))
names(trees_envs) <- levels(envs)
trees_envs['0'] <- 0

# gamma values to try Anchor Forest
gamma_seq <- c(0.0000001, 0.3, 0.6, seq(1, 4, 0.2))
gamma <- gamma_seq[as.numeric(args[1])]

# fit Anchor Forest
fit_anchor <- SDForest(x = X, y = Y, A = A, envs = envs, nTree_leave_out = trees_envs,
                       Q_type = "no_deconfounding", gamma = gamma ** 2, cp = 0, mc.cores = 20, gpu = F)

# collect out-of-environment residuals
pred <- fit_anchor$ooEnv_predictions
re <- (Y - pred)**2
save(re, gamma, file = paste0("results/anchorG/resid_", args[1], ".RData"))

end_time <- Sys.time()
cat("Time taken:", end_time - start_time, "\n")