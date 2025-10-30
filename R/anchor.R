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

n <- length(Y)
n
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
dfPerf_g <- gather(dfPerf, 'environment', 'mse', -gamma)
dfPerf_g$environment <- as.factor(dfPerf_g$environment)
levels(dfPerf_g$environment) <- 1:8



fit <- stats::loess(`#20` ~ gamma, dfPerf)
fitted <- predict(fit, newdata = data.frame(gamma = xseq))
# optimal gamma
xseq[which.min(fitted)]

perfRF <- dfPerf$`#20`[dfPerf$gamma == 1]
perfARF <- min(fitted)
(perfRF - perfARF) / perfRF

data2 <- data.frame(mse = fitted, gamma = xseq)

ggAnchor <- ggplot(dfPerf_g, aes(x = gamma, y = mse, group = environment)) + 
  theme_bw() + 
  geom_smooth(xseq = xseq, se = FALSE, linewidth = 0.4, 
              aes(linetype = environment, col = environment)) + 
  geom_point(aes(col = environment, shape = environment), size = 0.6) + 
  scale_shape_manual(values=0:7) +
  #geom_line(data = data2, colour = "blue", group = "3") + 
  geom_vline(xintercept = 1) +
  ylab("OOD MSE") + xlab(expression(gamma))
  # + coord_cartesian(ylim = c(6, 1), xlim = c(1, 4))

ggAnchor
ggsave("figures/AnchorCV.jpeg", width = 6, height = 3)

mean_perf <- sapply(unique(envs), function(env) mean((Y[envs == env] - mean(Y[envs != env]))**2))
mean_perf




load("results/anchor_opt/var_importance.RData")
plot(var_importance)
plot(sort(var_importance, decreasing = T))
imp_s <- which(var_importance >= sort(var_importance, decreasing = T)[40])

load("results/anchor_opt/regPath.RData")
gg_path <- plot(path, sqrt_scale = T)
gg_path + theme_bw() + theme(legend.position = "none")

cp <- which(path$cp == min(path$cp[path$cp > 0.5]))
path_s <- which(path$varImp_path[cp, ] > 0)

load("results/anchor_opt/stability_selection.RData")
plot(stab, sqrt_scale = T)
cp <- which(stab$cp == min(stab$cp[stab$cp > 0.5]))
stab_s <- which(stab$varImp_path[cp, ] > 0)

all(stab_s == path_s)
sum(stab_s %in% imp_s)

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

save(imp_s, path_s, stab_s, file = "results/anchor_opt/proteinSelection.RData")

load("results/anchor_opt/partial_dependence.RData")


plot(dep2)

