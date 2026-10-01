library(SDModels)
options(future.globals.maxSize = 30.0 * 1e9) # increase max global size for parallel processing
library(ggplot2)

print("load data")
load("results/anchorG_opt.RData")
print(fit_anchor)
var_importance <- fit_anchor$var_importance

print("Variable importance:")
save(var_importance, file = "results/anchor_opt/var_importance.RData")

# check whether enough trees were grown
print("plot anchor fit")
png("results/anchor_opt/anchor_plot.png", width = 800, height = 600)
plot(fit_anchor)
dev.off()

print("regpaths")
path <- regPath(fit_anchor)
save(path, file = "results/anchor_opt/regPath.RData")

print("stability selection")
stab <- stabilitySelection(fit_anchor)
save(stab, file = "results/anchor_opt/stability_selection.RData")

# three most important proteins
most_imp <- which(fit_anchor$var_importance >= sort(fit_anchor$var_importance, decreasing = TRUE)[3])

print("partial dependence")
dep1 <- partDependence(fit_anchor, most_imp[1])
dep2 <- partDependence(fit_anchor, most_imp[2])
dep3 <- partDependence(fit_anchor, most_imp[3])

save(dep1, dep2, dep3, file = "results/anchor_opt/partial_dependence.RData")