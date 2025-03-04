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


load('results/impProt.RData')
length(names(important_prot))
dim(X[, names(important_prot)])

perf <- sapply(Groups, function(group){
  fit_all <- ranger(x = X[envs != group, ], y = Y[envs != group])
  pred <- predict(fit_all, X[envs == group, ])
  1 - mean((pred$predictions - Y[envs == group])**2) / var(Y[envs == group])
})

perf

perf_imp <- sapply(Groups, function(group){
  fit_all <- ranger(x = X[envs != group, names(important_prot)], y = Y[envs != group])
  pred <- predict(fit_all, X[envs == group, names(important_prot)])
  1 - mean((pred$predictions - Y[envs == group])**2) / var(Y[envs == group])
})

perf_imp

perf_tree <- sapply(Groups, function(group){
  impData <- data.frame(Y = Y, X[, names(important_prot)])
  fit_simple <- rpart(Y ~ ., impData[envs != group, ], control = rpart.control(cp = 0, xval = 20))
  cp_min <- fit_simple$cptable[which.min(fit_simple$cptable[, 'xerror']), 'CP']
  pFit <- prune(fit_simple, cp_min)
  1 - mean((impData[envs == group, 'Y'] - predict(pFit, impData[envs == group, ]))**2) / var(impData[envs == group, 'Y'])
})

perf_tree
library(tidyr)

perf_df <- data.frame(full = perf, imp = perf_imp, tree = perf_tree)
perf_df$env <- as.factor(rownames(perf_df))
perf_df <- gather(perf_df, 'method', 'perf', -env)
perf_df$method <- as.factor(perf_df$method)

perf_df$method

plot(perf_df$method, perf_df$perf)

library(ggplot2)
ggplot(perf_df, aes(x = method, y = perf, group = env, col = env)) + 
  geom_line() + 
  geom_point() + 
  theme_bw()


library(rpart)
impData <- data.frame(Y = Y, X[, names(important_prot)])
names(impData) <- c('Y', 1:(ncol(impData)-1))

perf_tree_in <- sapply(1:nrow(impData), function(testInd){
  impData_test <- impData[testInd, ]
  impData_train <- impData[-testInd, ] 

  fit_simple <- rpart(Y ~ ., impData_trainl,m, control = rpart.control(cp = 0, xval = 20))
  cp_min <- fit_simple$cptable[which.min(fit_simple$cptable[, 'xerror']), 'CP']
  pFit <- rpart::prune(fit_simple, cp_min)
  (impData_test$Y - predict(pFit, impData_test))**2
})
perf_tree_in
1 - mean(perf_tree_in) / var(Y)


fit_simple <- rpart(Y ~ ., impData, control = rpart.control(cp = 0, xval = 20))
cp_min <- fit_simple$cptable[which.min(fit_simple$cptable[, 'xerror']), 'CP']
pFit <- rpart::prune(fit_simple, cp_min)
rattle::fancyRpartPlot(pFit, sub = '')


library(SDModels)
gamma_seq <- c(0.0000001, exp(seq(-2,1,0.5)), 5)
#gamma_seq <- c(0, 0.001, 0.5, 1, 1.5, 2, 5)

perf <- sapply(gamma_seq, function(gamma){
  fit_anchor <- SDForest(x = X[, names(important_prot)], y = Y, A = A, envs = envs, nTree_leave_out = trees_envs,
                         Q_type = "no_deconfounding", gamma = gamma ** 2, cp = 0.01, mc.cores = 8, gpu = T)

  pred <- fit_anchor$ooEnv_predictions
  re <- (Y - pred)**2
  tapply(re, envs, mean)
})

coll <- rgb(0.2,0.2,0.2,0.8)
matplot(gamma_seq, t(perf), type="l",col=coll,lty=1)

