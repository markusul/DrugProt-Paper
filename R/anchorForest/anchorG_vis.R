library(tidyr)
library(ggplot2)
library(gridExtra)
theme_set(theme_bw(base_size = 14))
library(SDModels)
set.seed(42)

# load and preprocess data
source("R/anchorForest/utils.R")
dat <- prepAggData()

Y <- dat$Y
A <- as.matrix(dat$A)
p_names <- dat$p_names
perturbations <- dat$perturbations
pertLabel <- dat$pertLabel

# number of samples
n <- length(Y)
n

# collect results
A_Results <- c("Results from Anchor Forest", "-----------------------------------")
A_Results <- c(A_Results, paste0("Number of samples: ", n))

# drug concentrations
sort(unique(as.numeric(A)))
A_Results <- c(A_Results, paste0("Unique drug concentrations: ", 
                                  paste(sort(unique(as.numeric(A))), collapse = ", ")))

# define environments based on drug doses for evaluation
anchor_dose <- sapply(strsplit(perturbations, ' '), function(x) x[1])
Groups <- unique(anchor_dose)

envs <- rep(0, n)
for(group in Groups){
  test_pert <- perturbations[anchor_dose == group]
  envs[which(pertLabel %in% test_pert)] <- group
}

# load results from different gamma values
preds <- list.files("results/anchorG")
res <- matrix(nrow = n, ncol = length(preds))
gamma_vec <- rep(NA, length(preds))

for(i in 1:length(preds)){
  load(paste0("results/anchorG/", preds[i]))
  res[, i] <- re
  gamma_vec[i] <- gamma
}

# compute OOD performance
mse <- apply(res, 2, function(resj) tapply(resj, as.factor(envs), mean))
mse <- mse[-which(rownames(mse) == "0"), ]
qPerf <- mse

dfPerf <- data.frame(t(qPerf))
names(dfPerf) <- rownames(mse)
dfPerf$gamma <- gamma_vec ** 2
dfPerf <- dfPerf[dfPerf$gamma >= 1, ]

# plot OOD performance vs gamma
xseq <- seq(min(dfPerf$gamma), max(dfPerf$gamma), length.out = 10000)
dfPerf_g <- gather(dfPerf, 'environment', 'mse', -gamma)
dfPerf_g$environment <- as.factor(dfPerf_g$environment)
levels(dfPerf_g$environment) <- 1:8

# most difficult environment, the one with the largest OOD MSE of the standard Random Forest (gamma = 1)
worst_env <- rownames(mse)[which.max(mse[, gamma_vec == 1])]
print(worst_env)

A_Results <- c(A_Results, paste0("Most difficult environment: ", worst_env))

# optimal gamma for most difficult environment
gamma_opt <- dfPerf[which.min(dfPerf[, worst_env]), "gamma"]
gamma_opt

A_Results <- c(A_Results, paste0("Optimal gamma for most difficult environment: ", 
                                  round(gamma_opt, 6)))

# anchorG_opt.R reads gamma_opt from here
dir.create("results/anchor_opt", showWarnings = FALSE, recursive = TRUE)
save(gamma_opt, worst_env, file = "results/anchor_opt/gamma_opt.RData")

# relative improvement over random forest (gamma = 1)
# relative improvement over random forest (gamma = 1) in the most difficult environment
perfRF  <- dfPerf[dfPerf$gamma == 1, worst_env]
perfARF <- min(dfPerf[, worst_env])
(perfRF - perfARF) / perfRF
A_Results <- c(A_Results, paste0("Relative improvement over random forest (gamma = 1): ", 
                                  round((perfRF - perfARF) / perfRF, 6)))


ggAnchor <- ggplot(dfPerf_g, aes(x = gamma, y = mse, group = environment)) + 
  theme_bw() + 
  geom_smooth(xseq = xseq, se = FALSE, linewidth = 0.4, 
              aes(linetype = environment, col = environment)) + 
  geom_point(aes(col = environment, shape = environment), size = 0.6) + 
  scale_shape_manual(values=0:7) +
  geom_vline(xintercept = 1) +
  ylab("OOD MSE") + xlab(expression(gamma))

ggAnchor
ggsave("figures/AnchorCV.jpeg", width = 6, height = 3)

# compute performance of mean prediction per environment
mean_perf <- sapply(unique(envs), function(env) mean((Y[envs == env] - mean(Y[envs != env]))**2))
mean_perf

A_Results <- c(A_Results, paste0("Performance of mean prediction per environment: ", 
                                  paste(round(mean_perf, 4), collapse = ", ")))

# standard Random Forest (gamma = 1) against the constant prediction by the training mean, per environment
env_names  <- rownames(mse)
rf_perf    <- unlist(dfPerf[dfPerf$gamma == 1, env_names])
const_perf <- mean_perf[env_names]
gain_rf    <- (const_perf - rf_perf) / const_perf
A_Results <- c(A_Results, paste0("Relative improvement of the standard Random Forest over the training mean per environment: ",
                                 paste(env_names, round(gain_rf, 3), sep = ": ", collapse = ", ")))
A_Results <- c(A_Results, paste0("Range of this improvement: ", round(min(gain_rf), 3), " to ", round(max(gain_rf), 3)))

# spread of the out-of-environment MSE over gamma >= 1 in the most difficult environment
A_Results <- c(A_Results, paste0("OOE MSE of the most difficult environment over gamma >= 1: ",
                                 round(min(dfPerf[, worst_env]), 3), " to ", round(max(dfPerf[, worst_env]), 3)))

# first pass after anchorG_CV.R ends here, second pass after anchorG_opt_res.R continues
if(!file.exists("results/anchor_opt/var_importance.RData")){
  message("gamma_opt saved, run anchorG_opt.R and anchorG_opt_res.R, then this script again")
  quit(save = "no")
}
# analyze variable importance and selected proteins of anchor forest with optimal gamma
load("results/anchor_opt/var_importance.RData")
plot(var_importance)
plot(sort(var_importance, decreasing = T))

times <- sapply(strsplit(names(sort(var_importance, decreasing = T)), '_'), function(x) x[length(x)])
par(mfrow = c(1, 4))
boxplot(which(times == "0"), ylim = c(0, 20000))
boxplot(which(times == "6"), ylim = c(0, 20000))
boxplot(which(times == "24"), ylim = c(0, 20000))
boxplot(which(times == "48"), ylim = c(0, 20000))

# 40 most important proteins
imp_s <- which(var_importance >= sort(var_importance, decreasing = T)[40])
# most important proteins
sort(var_importance, decreasing = T)[1:10]
  
# regularization path
load("results/anchor_opt/regPath.RData")
gg_path <- plot(path, sqrt_scale = T)
gg_path

plot(path, plotly = TRUE, imp_s)


plotOOB(path)

ggsave("figures/regPath.jpeg", 
       gg_path + theme_bw() + theme(legend.position = "none"), 
       width = 6, height = 3)

# selected proteins at cp > 0.5
cp <- which(path$cp == min(path$cp[path$cp > 0.5]))
path_s <- which(path$varImp_path[cp, ] > 0)

# stability selection
load("results/anchor_opt/stability_selection.RData")
plot(stab, sqrt_scale = T)

# selected proteins at cp > 0.5
cp <- which(stab$cp == min(stab$cp[stab$cp > 0.5]))
stab_s <- which(stab$varImp_path[cp, ] > 0)

# check overlaps
all(stab_s == path_s)
# stability selection results in the same proteins as regularization path

# overlap with important proteins
sum(stab_s %in% imp_s)

sort(var_importance[path_s])

##check timepoints of regularized selection
#number of important proteins
length(path_s)
A_Results <- c(A_Results, paste0("Number of selected proteins (regularization path, cp > 0.5): ", 
                                  length(path_s)))

# get the last elemetns after splitting by '_'
times <- table(sapply(strsplit(names(path_s), '_'), function(x) x[length(x)]))
times
A_Results <- c(A_Results, paste0("Timepoints of selected proteins (regularization path, cp > 0.5): ", 
                                  paste(names(times), times, sep = ": ", collapse = ", ")))

# map to short names
load('data/order.RData')

imp_to_shortnames <- function(imp){
  imp <- sapply(names(imp), function(p) {
    sp <- strsplit(p, '_')[[1]]
    paste(sp[-length(sp)], collapse =  "_")
  })
  prot_names_short[imp]
}

imp_s <- imp_to_shortnames(imp_s)
path_s <- imp_to_shortnames(path_s)
stab_s <- imp_to_shortnames(stab_s)

# number of unique proteins selected
unique(path_s)
A_Results <- c(A_Results, paste0("Number of unique selected proteins (regularization path, cp > 0.5): ", 
                                  length(unique(path_s))))
A_Results <- c(A_Results, paste0("Selected proteins (regularization path, cp > 0.5): ", 
                                  paste(unique(path_s), collapse = ", ")))

# most important proteins
most_imp <- sort(var_importance[var_importance > 2], decreasing = T)
# timepoints of most important proteins
times_imp <- sapply(strsplit(names(most_imp), '_'), function(x) x[length(x)])
most_imp <- imp_to_shortnames(most_imp)
most_imp_h <- paste(most_imp, paste0(times_imp, 'h'))

times_imp <- table(times_imp)

A_Results <- c(A_Results, paste0("Timepoints of most important proteins: ", 
                                  paste(names(times_imp), times_imp, sep = ": ", collapse = ", ")))

A_Results <- c(A_Results, paste0("Most important proteins: ", 
                                  paste(most_imp_h, collapse = ", ")))

# save most important proteins to txt as example
fileConn2 <- file("results/most_important_proteins.txt")
writeLines(most_imp_h, fileConn2)
close(fileConn2)

# save selected proteins
save(most_imp, imp_s, path_s, stab_s, file = "results/anchor_opt/proteinSelection.RData")

# largest variable importances, to support the gap between the eighth and ninth predictor
imp_sorted <- sort(var_importance, decreasing = TRUE)[1:10]
imp_time <- sapply(strsplit(names(imp_sorted), '_'), function(x) x[length(x)])
A_Results <- c(A_Results, paste0("Ten largest variable importances: ",
                                 paste(paste0(imp_to_shortnames(imp_sorted), " ", imp_time, "h"),
                                       round(imp_sorted, 2), sep = ": ", collapse = ", ")))

# partial dependence plots for the 2 most important proteins
load("results/anchor_opt/partial_dependence.RData")
deps <- list(dep1, dep2, dep3)
dep_imp <- sapply(deps, function(dep) var_importance[names(dep$j)])
stopifnot(!anyNA(dep_imp))
deps <- deps[order(dep_imp, decreasing = TRUE)[1:2]]

# label from the protein stored in the partial dependence object
dep_lab <- function(dep) {
  sp <- strsplit(names(dep$j), '_')[[1]]
  paste0(imp_to_shortnames(dep$j), " ", sp[length(sp)], "h")
}

# mean prediction at the 5%, 50% and 95% quantiles of the observed expression, for the direction in the text
dep_direction <- function(dep) {
  q <- quantile(dep$xj, c(0.05, 0.5, 0.95))
  f <- approx(dep$x_seq, dep$preds_mean, xout = q)$y
  paste0(dep_lab(dep), ": mean predicted IC50 ", paste(round(f, 2), collapse = ", "),
         " at the 5%, 50%, 95% quantiles of expression (", paste(round(q, 2), collapse = ", "), ")")
}
A_Results <- c(A_Results, "Partial dependence of the two most important proteins:",
               sapply(deps, dep_direction))

# common y range for both panels, x range of the central 98% of the data
# coord_cartesian zooms without dropping points, unlike xlim and ylim
y_range <- range(sapply(deps, function(dep) range(dep$preds_mean)))
pd_plot <- function(dep) {
  plot(dep, n_examples = 0) + xlab(dep_lab(dep)) + theme_bw() + ggtitle("") + ylab("") +
    coord_cartesian(xlim = quantile(dep$xj, c(0.01, 0.99)), ylim = y_range)
}
gg1 <- pd_plot(deps[[1]]) + ylab(expression(widehat(IC50)))
gg2 <- pd_plot(deps[[2]])

ggdep <- grid.arrange(gg1, gg2, nrow = 1, bottom = "Partial dependence on differential log-protein expression")
ggsave("figures/AnchorDep.jpeg", ggdep, width = 7, height = 3)

# file to save results to
fileConn <- file("results/A_Results.txt")
writeLines(A_Results, fileConn)
close(fileConn)
