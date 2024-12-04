library(ranger)
library(SDForest)
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

fit_full <- ranger(y = Y, x = X, importance = 'impurity')
fit_full
plot(sort(fit_full$variable.importance))

imp_full <- names(sort(fit_full$variable.importance, decreasing = T))
imp_full


i <- 1
perf <- sapply(seq(1, 1000, 10), function(i){
  fit_lim <- ranger(y = Y, x = X[, imp_full[-c(1:i)]])
  fit_lim$r.squared
})

plot(seq(1, 1000, 10), perf)

