args = commandArgs(trailingOnly = TRUE)
set.seed(42)

source("R/utils.R")

dat <- prepAggData()
X <- as.matrix(dat$X)
Y <- dat$Y
A <- as.matrix(dat$A)
p_names <- dat$p_names
perturbations <- dat$perturbations
pertLabel <- dat$pertLabel

singleDrugs <- !pertLabel %in% perturbations

X <- X[singleDrugs,]
Y <- Y[singleDrugs]
A <- A[singleDrugs,]
pertLabel <- pertLabel[singleDrugs]

X <- X[, apply(X, 2, function(x)length(unique(x))) >= 200]
n <- length(Y)

Groups <- unique(pertLabel)
envs <- rep(0, n)

for(test_pert in Groups){
  envs[which(pertLabel %in% test_pert)] <- test_pert
}

envs <- as.factor(envs)

# ranomly put envs in group
envs <- as.factor(as.numeric(envs))
# sample level ids
levels(envs) <- sample(levels(envs))

envs <- as.character(envs)
envs <- as.numeric(envs) %% 20 + 1
envs <- as.factor(envs)

trees_envs <- rep(50, length(levels(envs)))
names(trees_envs) <- levels(envs)

library(SDModels)
gamma_seq <- c(0.0000001, exp(seq(-2,1,0.5)), 5)

gamma <- gamma_seq[as.numeric(args[1])]

fit_anchor <- SDForest(x = X, y = Y, A = A, envs = envs, nTree_leave_out = trees_envs,
                       Q_type = "no_deconfounding", gamma = gamma ** 2, cp = 0, mc.cores = 20, gpu = F)

pred <- fit_anchor$ooEnv_predictions
re <- (Y - pred)**2
save(re, gamma, file = paste0("results/anchorS/resid_", args[1], ".RData"))
print("done")