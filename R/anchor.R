library(tidyr)
library(ggplot2)

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
perturbations

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

mse <- apply(res, 2, function(resj) tapply(resj, as.factor(envs), mean))
mse <- mse[-which(rownames(mse) == "0"), ]
qPerf <- mse

dfPerf <- data.frame(t(qPerf))
names(dfPerf) <- rownames(mse)
dfPerf$gamma <- gamma_vec ** 2

xseq <- seq(min(dfPerf$gamma), max(dfPerf$gamma), length.out = 10000)
dfPerf_g <- gather(dfPerf, 'intervention', 'mse', -gamma)
#dfPerf_g <- dfPerf_g[dfPerf_g$intervention == "#20", ]

fit <- stats::loess(`#20` ~ gamma, dfPerf)
fitted <- predict(fit, newdata = data.frame(gamma = xseq))
xseq[which.min(fitted)]
data2 <- data.frame(mse = fitted, gamma = xseq)

ggplot(dfPerf_g, aes(x = gamma, y = mse, group = intervention)) + 
  theme_bw() + 
  geom_smooth(xseq = xseq, se = FALSE, aes(linetype = intervention, col = intervention)) + 
  geom_point(aes(col = intervention, shape = intervention)) + 
  scale_shape_manual(values=0:7) +
  #geom_line(data = data2, colour = "blue", group = "3") + 
  geom_vline(xintercept = 1) +
  ylab("OOB MSE") + xlab(expression(gamma))
  # + coord_cartesian(ylim = c(6, 1), xlim = c(1, 4))

mean_perf <- sapply(unique(envs), function(env) mean((Y[envs == env] - mean(Y[envs != env]))**2))
mean_perf


load("results/anchorG_opt.RData")
plot(log(fit_anchor$var_importance))

fit_anchor <- fromList(fit_anchor)
plot(fit_anchor)

path <- regPath(fit_anchor)
plot(path)


