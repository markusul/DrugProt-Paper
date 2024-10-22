library(ranger)
library(SDForest)
source("R/utils.R")

dat <- prepAggData()
X <- dat$X
Y <- dat$Y
A <- dat$A
p_names <- dat$p_names
perturbations <- dat$perturbations
pertLabel <- dat$pertLabel


anchor_dose <- sapply(strsplit(perturbations, ' '), function(x) x[1])
Groups <- unique(anchor_dose)

Groups
anchor_dose
perturbations
unique(pertLabel)
A
dim(X)

start <- Sys.time()
fit_sdf <- SDForest(x = X, y = Y, nTree = 100, return_data = F, multicore = F, 
                    mem_size = 4e+7, leave_out_ind = test_env, max_size = 100, 
                    envs = )
end <- Sys.time()
print(paste0("Training time: ", end - start))

print(paste0("MSE: ", fit_sdf$oob_loss))
print(paste0("SDE: ", fit_sdf$oob_SDloss))