library(plotly)

load('data/order.RData')

# protein of interest
P <- 1
t <- 6
alpha <- 0.5


load(file = paste0('pvals/', P , '_', t, '.RData'))
pMat <- as.matrix(pMat)
rownames(pMat) <- colnames(pMat)

pMat[is.na(pMat)] <- 2

pMat <- pMat[drugOrder, drugOrder]
#pMat <- -log(pMat)
pMat <- round(pMat, 3)
min(pMat)


ht <- plot_ly(z = pMat, x = colnames(pMat), y = colnames(pMat), 
              type = "heatmap", colors = "Greys") %>%
              layout(title = prot_names_short[which(prot_names == P)])
ht

htmlwidgets::saveWidget(as_widget(ht), paste0(prot_names_short[which(prot_names == P)], ".html"))


plot(pval.corr, main = prot_names_short[which(prot_names == P)])

prot_names_short[which(prot_names == P)]
unname(prot_names_short[names(pval.corr[pval.corr < alpha])])

t <- 6
load(file = paste0('models/', P , '_', t, '.RData'))

plot(Y - design %*% fit$bhat, x = design %*% fit$bhat)
car::qqPlot(design %*% fit$bhat)
plot(sqrt(abs(Y - design %*% fit$bhat)), x = design %*% fit$bhat)



library(hdi)
?lasso.proj

res <- fit$clusterGroupTest()
names(res)






res$pval

unname(fit$pval)

plot(res)
res$pval
res$clusters

resfit$pval[1:2]
fit$pval.corr[1]
fit$groupTest(1:2, F)
?fit$clusterGroupTest()
res$rightCh
plot(res)


Net <- lapply(c(24, 48), function(t){
  net <- NULL
  targets <- c()
  
  for (P in 1:length(prot_names)) {
    path <- paste0('pvals/', P , '_', t, '.RData')
    if(file.exists(path)){
      load(file = path)
      net <- cbind(net, pval.corr[(length(pval.corr)-length(prot_names_short) + 1):length(pval.corr)])
      targets <- c(targets, prot_names_short[P])
    }
  }
  list(net, targets)
})

Links_all <- lapply(Net, function(net){
  res <- apply(net[[1]], 2, function(pval) which(pval < alpha))
  links <- NULL
  for(i in 1:ncol(net[[1]])){
    if(length(res[[i]]) != 0)
      links <- rbind(links, data.frame(res[[i]], which(net[[2]][i] == prot_names_short)))
  }
  rownames(links) <- 1:nrow(links)
  colnames(links) <- c('source', 'target')
  links
})


library(networkD3)

alpha <- 0.8

# summary graph
Links_sum <- do.call(rbind, Links_all)
Links_sum$source <- Links_sum$source - 1
Links_sum$target <- Links_sum$target - 1
Links_sum$value <- 1
Nodes_sum <- data.frame(name = prot_names_short, group = 1, size = 1)

forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99, 
             arrows = T, zoom = T, charge = -100,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))


# temporal graph
expTimes <- c(6, 24, 48)
nodenames <- c(paste(prot_names_short, expTimes[1], sep = '_'), 
               paste(prot_names_short, expTimes[2], sep = '_'), 
               paste(prot_names_short, expTimes[3], sep = '_'))
nodegroups <- rep(paste0(expTimes, 'h'), each = length(prot_names_short))

Links_temp <- lapply(1:2, function(t){
  links <- Links_all[[t]]
  links[, 1] <- links[, 1] + (t-1) * length(prot_names_short)
  links[, 2] <- links[, 2] + (t) * length(prot_names_short)
  links
})
Links_temp <- do.call(rbind, Links_temp)
used_nodes <- sort(unique(unname(unlist(Links_temp))))

Links_temp$source <- Links_temp$source - 1
Links_temp$target <- Links_temp$target - 1
Links_temp$value <- 1
Nodes_temp <- data.frame(name = nodenames, group = nodegroups, size = 1)
Nodes_temp$size[used_nodes] <- 100

forceNetwork(Links = Links_temp, Nodes = Nodes_temp,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99, Nodesize = 3,
             arrows = T, zoom = T, legend=T, charge = -100,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))

#ht
#htmlwidgets::saveWidget(as_widget(ht), "temporalGraph.html")
