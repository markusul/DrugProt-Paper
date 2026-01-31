library(tidyr)
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(SDModels)
set.seed(42)

# load and preprocess data
source("R/utils.R")
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

# plot OOD performance vs gamma
xseq <- seq(min(dfPerf$gamma), max(dfPerf$gamma), length.out = 10000)
dfPerf_g <- gather(dfPerf, 'environment', 'mse', -gamma)
dfPerf_g$environment <- as.factor(dfPerf_g$environment)
levels(dfPerf_g$environment) <- 1:8

# smooth worst case environment curve
fit <- stats::loess(`#20` ~ gamma, dfPerf)
fitted <- predict(fit, newdata = data.frame(gamma = xseq))
# optimal gamma for most difficult environment
xseq[which.min(fitted)]
A_Results <- c(A_Results, paste0("Optimal gamma for most difficult environment: ", 
                                  round(xseq[which.min(fitted)], 6)))

# relative improvement over random forest (gamma = 1)
perfRF <- dfPerf[dfPerf$gamma == 1, "#20"]
perfARF <- min(fitted)
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

# analyze variable importance and selected proteins of anchor forest with optimal gamma
load("results/anchor_opt/var_importance.RData")
plot(var_importance)
plot(sort(var_importance, decreasing = T))

# 40 most important proteins
imp_s <- which(var_importance >= sort(var_importance, decreasing = T)[40])
# 3 most important proteins
sort(var_importance, decreasing = T)[1:3]
  
# regularization path
load("results/anchor_opt/regPath.RData")
gg_path <- plot(path, sqrt_scale = T)
gg_path

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

##check timeoints of regularized selection
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

# 3 most important proteins
most_imp <- which(var_importance >= sort(var_importance, decreasing = TRUE)[3])
# timepoints of most important proteins
times_imp <- sapply(strsplit(names(most_imp), '_'), function(x) x[length(x)])
most_imp <- imp_to_shortnames(most_imp)
most_imp <- paste(most_imp, paste0(times_imp, 'h'))

times_imp <- table(times_imp)

A_Results <- c(A_Results, paste0("Timepoints of 3 most important proteins: ", 
                                  paste(names(times_imp), times_imp, sep = ": ", collapse = ", ")))

A_Results <- c(A_Results, paste0("3 most important proteins: ", 
                                  paste(most_imp, collapse = ", ")))

# save most important proteins to txt as example
fileConn2 <- file("results/most_important_proteins.txt")
writeLines(most_imp, fileConn2)
close(fileConn2)

# save selected proteins
save(most_imp, imp_s, path_s, stab_s, file = "results/anchor_opt/proteinSelection.RData")

# partial dependence plots for the 3 most important proteins
load("results/anchor_opt/partial_dependence.RData")
gg3 <- plot(dep3) + xlab(most_imp[3]) + theme_bw() + ylim(5, 8) + xlim(-1, 1) + 
  ggtitle("") + ylab("") 
gg2 <- plot(dep2) + xlab(most_imp[2]) + theme_bw() + ylim(5, 8) + xlim(-1, 1) + 
  ggtitle("") + ylab("")
gg1 <- plot(dep1) + xlab(most_imp[1]) + theme_bw() + ylim(5, 8) + xlim(-1, 1) + 
  ylab(expression(widehat(IC50)))

library(gridExtra)
ggdep <- grid.arrange(gg1, gg2, gg3, nrow =1)
ggsave("figures/AnchorDep.jpeg", ggdep, width = 7, height = 3)

# comparison to plain model
load("results/plainFit_full.RData")
imp_plain_full <- sort(fit_plain_full$var_importance, decreasing = T)
imp_plain_full[1:10]

plot(var_importance, fit_plain_full$var_importance)

# file to save results to
fileConn <- file("results/A_Results.txt")
writeLines(A_Results, fileConn)
close(fileConn)
