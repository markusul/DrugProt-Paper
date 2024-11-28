library(SDForest)

load('results/fit_plain.Rdata')
load('results/fit_sdf.Rdata')

fit_sdf <- fromList(fit_sdf)

path_sdf <- regPath(fit_sdf)
plotOOB(path_sdf)

fit_sdf <- prune(fit_sdf, path_sdf$cp_min)


imp_plain <- fit_plain$var_importance
imp_sdf <- fit_sdf$var_importance

plot(log(imp_sdf), log(imp_plain))

length(imp_sdf)
var_names <- fit_sdf$var_names
var_names <- sub("_[0-9]+$", "", var_names)

imp <- data.frame(sdf = imp_sdf / max(imp_sdf), 
                  plain = imp_plain / max(imp_plain))

agg_imp <- aggregate(imp, by = list(var_names), FUN = sum)
agg_imp

library(tidyr)
short_agg_imp <- gather(agg_imp, key = 'type', value = 'importance', -Group.1)

library(ggplot2)
ggplot(agg_imp, aes(x = plain, y = sdf)) +
  geom_point() + 
  theme_bw()

ggplot(agg_imp, aes(x = log(plain), y = log(sdf))) +
  geom_point() + 
  theme_bw()

mean(agg_imp$sdf == 0)
