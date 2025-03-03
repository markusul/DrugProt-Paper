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

fit_all <- ranger(x = X, y = Y)
fit_all

fit_imp <- ranger(x = X[, names(important_prot)], y = Y)
fit_imp

library(rpart)
impData <- data.frame(Y = Y, X[, names(important_prot)])
names(impData) <- c('Y', 1:(ncol(impData)-1))

testInd <- sample(1:nrow(impData), 300)
impData_test <- impData[testInd, ]
impData_train <- impData[-testInd, ] 

fit_simple <- rpart(Y ~ ., impData_train, control = rpart.control(cp = 0))
cp_min <- fit_simple$cptable[which.min(fit_simple$cptable[, 'xerror']), 'CP']
pFit <- prune(fit_simple, cp_min)

library(rattle)
fancyRpartPlot(pFit, sub = '')

1 - mean((impData_test$Y - predict(pFit, impData_test))**2) / var(impData_test$Y)

library(SDModels)
gamma_seq <- c(0.0000001, exp(seq(-2,1,0.5)))
#gamma_seq <- c(0, 0.001, 0.5, 1, 1.5, 2, 5)

perf <- sapply(gamma_seq, function(gamma){
  fit_anchor <- SDForest(x = X[, names(important_prot)], y = Y, A = A, envs = envs, nTree_leave_out = trees_envs,
                         Q_type = "no_deconfounding", gamma = gamma ** 2, cp = 0.001, max_size = 500)

  pred <- fit_anchor$ooEnv_predictions
  re <- (Y - pred)**2+
  tapply(re, envs, mean)
})

coll <- rgb(0.2,0.2,0.2,0.8)
matplot(gamma_seq, t(perf), type="l",col=coll,lty=1)

