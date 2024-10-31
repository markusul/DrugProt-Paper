#library(ranger)
library(SDForest)
source("R/utils.R")

dat <- prepAggData()
X <- as.matrix(dat$X)
Y <- dat$Y
A <- dat$A
p_names <- dat$p_names
perturbations <- dat$perturbations
pertLabel <- dat$pertLabel

#hist(apply(X, 2, function(x)length(unique(x))), breaks = 1000)

X <- X[, apply(X, 2, function(x)length(unique(x))) >= 200]
n <- length(Y)

anchor_dose <- sapply(strsplit(perturbations, ' '), function(x) x[1])
Groups <- unique(anchor_dose)

envs <- rep(0, n)

for(group in Groups){
  test_pert <- perturbations[anchor_dose == group]
  envs[which(pertLabel %in% test_pert)] <- group
}
envs <- as.factor(envs)

start <- Sys.time()
fit_sdf <- SDForest(x = X, y = Y, envs = envs, nTree_leave_out = 10, mc.cores = 100)
fit_sdf <- toList(fit_sdf)
end <- Sys.time()
print(paste0("Training time SDF: ", end - start))

start <- Sys.time()
fit_plain <- SDForest(x = X, y = Y, envs = envs, nTree_leave_out = 10, mc.cores = 100,
                    Q_type = "no_deconfounding")
fit_plain <- toList(fit_plain)
end <- Sys.time()
print(paste0("Training time PLAIN: ", end - start))

save(fit_sdf, fit_plain, fitle = "results/fits.Rdata")