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

# optimal gamma from cross-validation
gamma <- 7.965597

# fit Anchor Forest
options(future.globals.maxSize = 2.0 * 1e9) # increase max global size for parallel processing
fit_anchor <- SDForest(x = X, y = Y, A = A, nTree = 1000,
                       Q_type = "no_deconfounding", gamma = gamma, 
                       cp = 0, mc.cores = 100)
save(fit_anchor, gamma, file = paste0("results/anchorG_opt.RData"))

end_time <- Sys.time()
cat("Time taken:", end_time - start_time, "\n")