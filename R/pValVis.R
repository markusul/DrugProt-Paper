library(plotly)

load('data/order.RData')

# protein of interest
P <- 3753
t <- 24
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


load(file = paste0('models/', P , '_', t, '.RData'))
fit$
res <- fit$clusterGroupTest()
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


library(networkD3)

alpha <- 0.8
P <- 2

which(targets == prot_names_short[P])

ancPval <- Net[[1]][[1]][, targets == prot_names_short[P]]
decPval <- Net[[2]][[1]][prot_names[P], ]

min(ancPval)

anc <- which(ancPval < alpha)
dec <- which(decPval < alpha) + length(prot_names_short)

source <- unname(c(anc, rep(P, length(dec)))) - 1
target <- c(rep(P, length(anc)), dec) - 1

nodenames <- c(prot_names_short, paste0(Net[[2]][[2]], '_48'))
groups <- c(rep(0, length(prot_names_short)), rep(10, length(Net[[2]][[2]])))

Links <- data.frame(source = source, target = target, value = 1)
Nodes <- data.frame(name = nodenames, group = groups, size = 1)

tail(Nodes)

forceNetwork(Links = Links, Nodes = Nodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99, 
             arrows = T, zoom = T, 
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))


# summary graph

net <- Net[[1]]

Links <- lapply(Net, function(net){
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
Links <- do.call(rbind, Links)
Links$source <- Links$source - 1
Links$target <- Links$target - 1
Links$value <- 1
Nodes <- data.frame(name = prot_names_short, group = 1, size = 1)


Links$source == 2

forceNetwork(Links = Links, Nodes = Nodes,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99, 
             arrows = T, zoom = T, 
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
