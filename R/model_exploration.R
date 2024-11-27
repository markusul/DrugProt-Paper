library(SDForest)
library(ranger)

source("R/utils.R")

dat <- prepAggData()
X <- as.matrix(dat$X)
Y <- dat$Y
A <- as.matrix(dat$A)
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
trees_envs <- rep(50, length(levels(envs)))
names(trees_envs) <- levels(envs)
trees_envs['0'] <- 0

fit_plain <- ranger(x = X, y = Y)


load('results/fit_sdf.Rdata')
fit_sdf <- fromList(fit_sdf)

paths <- regPath(fit_sdf)
plotOOB(paths)

plot(paths)
