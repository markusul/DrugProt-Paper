library(SDModels)

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

# optimal gamma from cross-validation
gamma <- 3.443544

# fit Anchor Forest
fit_anchor <- SDForest(x = X, y = Y, A = A, nTree = 10,
                       Q_type = "no_deconfounding", gamma = gamma, 
                       cp = 0, mc.cores = 20, gpu = F)
# save fitted model as list
fit_anchor <- toList(fit_anchor)
save(fit_anchor, gamma, file = paste0("results/anchorG_opt.RData"))

end_time <- Sys.time()
cat("Time taken:", end_time - start_time, "\n")