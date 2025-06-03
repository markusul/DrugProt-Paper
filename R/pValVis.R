library(plotly)

load('data/order.RData')

# protein of interest
P <- 3753
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


load(file = paste0('models/', P , '_', t, '.RData'))
fit$
res <- fit$clusterGroupTest()
res$rightCh
plot(res)

P <- 1


t <- 24

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
dim(net)
net <- net
rownames(net) <- prot_names_short
colnames(net) <- targets

net <- t(net)

net <- net < alpha
net <- net[, colSums(net) > 0]
missingRows <- colnames(net)[!colnames(net) %in% rownames(net)]
netRows <- matrix(2, ncol = ncol(net), nrow = length(missingRows))
rownames(netRows) <- missingRows
net <- rbind(net, netRows)

missingCols <- rownames(net)[!rownames(net) %in% colnames(net)]
netCols <- matrix(2, ncol = length(missingCols), nrow = nrow(net))
colnames(netCols) <- missingCols
net <- cbind(net, netCols)

withEdge <- colSums(net < alpha) + rowSums(net < alpha) > 0
net <- net[withEdge, withEdge]


Graph <- igraph::graph_from_adjacency_matrix(t(net<alpha))
plot(Graph, layout = layout_in_circle(net.bg))
