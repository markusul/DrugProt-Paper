library(SDModels)

set.seed(42)

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

preds <- list.files("results/anchorG")

res <- matrix(nrow = n, ncol = length(preds))
gamma_vec <- rep(NA, length(preds))

for(i in 1:length(preds)){
  load(paste0("results/anchorG/", preds[i]))
  res[, i] <- re
  gamma_vec[i] <- gamma
}

gamma_vec


mse <- apply(res, 2, function(resj) tapply(resj, as.factor(envs), mean))
mse <- mse[-which(rownames(mse) == "0"), ]
dim(mse)

quantilevec <- 0.005*(1:199)


#qPerf <- apply(mse, 2, quantile, probs = quantilevec)
qPerf <- mse

dfPerf <- data.frame(t(qPerf))
#names(dfPerf) <- quantilevec
names(dfPerf) <- rownames(mse)
dfPerf$gamma <- gamma_vec


library(tidyr)
library(ggplot2)

str(dfPerf)

dfPerf_g <- gather(dfPerf, 'quantile', 'mse', -gamma)
ggplot(dfPerf_g, aes(x = gamma, y = mse, group = quantile)) + 
  theme_bw() + 
  geom_line(col = rgb(0.2,0.2,0.2,0.8)) + 
  geom_vline(xintercept = 1)

gamma_vec

zero_model <- tapply(Y**2, envs, mean)

mean_model <- sapply(unique(envs), function(env){
  mean((Y[envs == env] - mean(Y[envs != env]))**2)
})


zero_perf <- dfPerf
zero_perf[, 1:8] <- t(t(dfPerf[, 1:8]) / c(zero_model[names(dfPerf)[1:8]]))
zero_perf_g <- gather(zero_perf, 'quantile', 'mse', -gamma)

ggplot(zero_perf_g, aes(x = gamma, y = mse, group = quantile)) + 
  theme_bw() + 
  geom_line(col = rgb(0.2,0.2,0.2,0.8)) + 
  geom_vline(xintercept = 1) + 
  ggtitle('MSE/zero model')

mean_perf <- dfPerf
mean_perf[, 1:8] <- t(t(dfPerf[, 1:8]) / c(mean_model[names(dfPerf)[1:8]]))
mean_perf_g <- gather(mean_perf, 'quantile', 'mse', -gamma)

ggplot(mean_perf_g, aes(x = gamma, y = mse, group = quantile)) + 
  theme_bw() + 
  geom_line(col = rgb(0.2,0.2,0.2,0.8)) + 
  geom_vline(xintercept = 1) +
  ggtitle('MSE/mean model') + 
  xlim(c(1, 5))

