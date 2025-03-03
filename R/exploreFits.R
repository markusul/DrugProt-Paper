library(SDModels)

load('results/fit_plain.Rdata')
load('results/fit_sdf.Rdata')
load('results/fit_anchor.Rdata')

imp_plain <- fit_plain$var_importance
imp_sdf <- fit_sdf$var_importance
imp_anchor <- fit_anchor$var_importance

important_prot <- sort(imp_anchor, decreasing = T)[1:(0.01*length(imp_anchor)/3)]
save(important_prot, file = 'results/impProt.RData')

plot(log(imp_plain), log(imp_sdf))
plot(log(imp_plain), log(imp_anchor))
plot(log(imp_anchor), log(imp_sdf))

length(imp_sdf)
var_names <- fit_sdf$var_names
var_names <- sub("_[0-9]+$", "", var_names)

imp <- data.frame(sdf = imp_sdf / max(imp_sdf), 
                  plain = imp_plain / max(imp_plain), 
                  anchor = imp_anchor / max(imp_anchor))

agg_imp <- aggregate(imp, by = list(var_names), FUN = sum)
agg_imp

library(tidyr)
short_agg_imp <- gather(agg_imp, key = 'type', value = 'importance', -Group.1)

library(ggplot2)
gg_imp1 <- ggplot(agg_imp, aes(x = log(plain), y = log(sdf))) +
  geom_point(size =  0.1) + 
  theme_bw()
gg_imp1
ggsave('figures/imp1.png', gg_imp1, width = 4, height = 3)

q <- 0.99
gg_imp2 <- ggplot(agg_imp, aes(x = log(plain), y = log(anchor))) +
  geom_point(size =  0.1) + 
  theme_bw() + annotate("rect", 
                        xmin = min(log(agg_imp$plain)), 
                        xmax = max(log(agg_imp$plain)), 
                        ymin = quantile(log(agg_imp$anchor), q), 
                        ymax = max(log(agg_imp$anchor)),
                     alpha = .2, fill = 'green') + 
  annotate("rect", 
           xmin = quantile(log(agg_imp$plain), q), 
           xmax = max(log(agg_imp$plain)), 
           ymin = min(log(agg_imp$anchor)), 
           ymax = quantile(log(agg_imp$anchor), q),
           alpha = .2, fill = 'orange')
gg_imp2
ggsave('figures/imp2.png', gg_imp2, width = 4, height = 3)

gg_imp3 <- ggplot(agg_imp, aes(x = log(sdf), y = log(anchor))) +
  geom_point(size = 0.1) + 
  theme_bw()
gg_imp3
ggsave('figures/imp3.png', gg_imp3, width = 4, height = 3)




plot(fit_anchor$ooEnv_predictions, Y)

err_outEnv <- data.frame(anchor = (fit_anchor$ooEnv_predictions - fit_anchor$Y)**2, 
                         sdf = (fit_sdf$ooEnv_predictions - fit_sdf$Y)**2, 
                         plain = (fit_plain$ooEnv_predictions - fit_plain$Y)**2, 
                         env = fit_plain$envs)
err_outEnv <- gather(err_outEnv, value = 'error', key = 'method', -env)

ggplot(err_outEnv, aes(y = error, x = env, col = method)) + 
  geom_boxplot()


err_OOB <- data.frame(anchor = (fit_anchor$oob_predictions - fit_anchor$Y)**2, 
                         sdf = (fit_sdf$oob_predictions - fit_sdf$Y)**2, 
                         plain = (fit_plain$oob_predictions - fit_plain$Y)**2, 
                         env = fit_plain$envs)
err_OOB <- gather(err_OOB, value = 'error', key = 'method', -env)

ggplot(err_OOB, aes(y = error, x = env, col = method)) + 
  geom_boxplot()


err_outEnv$type <- 'oEnv'
err_OOB$type <- 'oBag'

err <- rbind(err_outEnv, err_OOB)

ggplot(err, aes(y = error, x = env, col = type)) + 
  geom_boxplot() + facet_grid(~method) + 
  xlab('anchor drug of combination') + 
  teh
