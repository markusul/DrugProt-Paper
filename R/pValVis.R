library(plotly)

load('data/order.RData')
nDrugs <- length(drugOrder)
# protein of interest
P <- 440

allPvecs <- lapply(1:length(prot_names_short), function(P){
pvec <- sapply(c(6, 24, 48), function(t){
  if(file.exists(paste0('results/DrugEffects/', P , '_', t, '.RData'))){
    load(file = paste0('results/DrugEffects/', P , '_', t, '.RData'))
    pval.drugs <- p.adjust(pval.drugs, method = 'holm')
    return(pval.drugs)
  }else{
    print(paste0(P, '_', t, " not found!"))
    return(rep(1, 122))
  }
})
if(!is.null(dim(pvec))) pvec <- apply(pvec, 1, function(p) min(min(p) * 3, 1))
pvec
})
allPvecs


allPvecs <- allPvecs[unlist(lapply(allPvecs, length)) == 122]
allPvecs <- do.call(rbind, allPvecs)
dim(allPvecs)
length(prot_names_short)

d <- dist(allPvecs)
fit <- hclust(d)
plot(fit)
clusters <- cutree(fit, 7)

pMat <- matrix(NA, nrow = nDrugs, ncol = nDrugs)
rownames(pMat) <- colnames(pMat) <- names(pvec)[1:nDrugs]
for(l in names(pvec)){
  drugs <- strsplit(l, ":")[[1]]
  if(length(drugs) == 1) drugs <- c(drugs, drugs)
  pMat[drugs[1], drugs[2]] <- pvec[l]
  pMat[drugs[2], drugs[1]] <- pvec[l]
}

pMat <- as.matrix(pMat)
#rownames(pMat) <- colnames(pMat)

pMat[is.na(pMat)] <- 2

pMat <- pMat[drugOrder, drugOrder]
#pMat <- -log(pMat)
#pMat <- round(pMat, 3)
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
    path <- paste0('results/ProteinEffects/', P , '_', t, '.RData')
    if(file.exists(path)){
      load(file = path)
      pval.corr <- p.adjust(pval[(length(pval)-length(prot_names_short) + 1):length(pval)], 
                            method = "BH")
      net <- cbind(net, pval.corr)
      targets <- c(targets, prot_names_short[P])
    }
  }
  list(net, targets)
})

alpha <- 0.05 / length(prot_names_short)

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

save(Net, Links_all, file = "results/proteinNetwork.RData")
load("results/proteinNetwork.RData")


library(networkD3)
P <- c(98, 33, 134)


# summary graph
Links_sum <- do.call(rbind, Links_all)
Links_sum$source <- Links_sum$source - 1
Links_sum$target <- Links_sum$target - 1
Links_sum$value <- 1
Nodes_sum <- data.frame(name = prot_names_short, group = clusters, size = 1)

P.ind <- P - 1
rel.Links <- Links_sum$source %in% P.ind | Links_sum$target %in% P.ind
Links_sum <- Links_sum[rel.Links, ]
rel.Nodes <- sort(unique(unlist(Links_sum[, c('source', 'target')])))
Nodes_sum <- Nodes_sum[rel.Nodes+1, ]

#reorganize link index
for(i in 1:length(rel.Nodes)){
  Links_sum[Links_sum == rel.Nodes[i]] <- i - 1
}
Links_sum$value <- 1


fN <- forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99, 
             arrows = T, zoom = T, charge = -20,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
fN
htmlwidgets::saveWidget(as_widget(fN), "SummaryGraph.html")


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
Nodes_temp$radius <- as.numeric(rep(1:length(prot_names_short), 3))

P.ind <- unlist(lapply(P, function(p) p + 0:2 * length(prot_names_short))) - 1
rel.Links <- Links_temp$source %in% P.ind | Links_temp$target %in% P.ind
Links_temp <- Links_temp[rel.Links, ]
rel.Nodes <- sort(unique(unlist(Links_temp[, c('source', 'target')])))
Nodes_temp <- Nodes_temp[rel.Nodes+1, ]

#reorganize link index
for(i in 1:length(rel.Nodes)){
  Links_temp[Links_temp == rel.Nodes[i]] <- i - 1
}
Links_temp$value <- 1


library(HiveR)
require("grid")
edges <- Links_temp
edges[, 1:2] <- edges[, 1:2] + 1
names(edges) <- c("id1", "id2", "weight")
# color
row.names(edges) <- NULL
edges$id1 <- as.integer(edges$id1)
edges$id2 <- as.integer(edges$id2)
edges$color <- "black"
edges$weight <- 0.1

str(edges)

nodes <- Nodes_temp
nodes$size <- 0.01
names(nodes) <- c("lab", "axis", "size", "radius")
#nodes$axis <- as.integer(as.numeric(nodes$axis == "24h") + 1)
nodes$axis[nodes$axis == "6h"] <- 2
nodes$axis[nodes$axis == "24h"] <- 1
nodes$axis[nodes$axis == "48h"] <- 3
nodes$axis <- as.integer(nodes$axis)
nodes$id <- 1:nrow(nodes)
nodes$radius <- nodes$radius * 3
nodes$color <- c("black", "red", "blue", "green", "violet", "yellow", "darkgreen")[clusters]
rep(c("black", "red", "blue", "green", "violet", "yellow", "darkgreen"), times = as.numeric(table(clusters)))



#sort by cliusters
nodes$radius <- sort_by(nodes$radius * 3, rep(clusters, 3))
length(nodes$radius)
length(clusters) * 3

str(nodes)

HEC <- list()
HEC$nodes <- nodes
HEC$edges <- edges
HEC$type <- "2D"
HEC$desc <- "HairEyeColor data set"
HEC$axis.cols <- c("grey", "grey")
class(HEC) <- "HivePlotData"

chkHPD(HEC) # answer of FALSE means there are no problems
#sumHPD(HEC)

plotHive(HEC, ch = 0.001, bkgnd = "white", 
         axLabs = c("24h", "6h", "48h"), 
         axLab.gpar = gpar(col = "black", fontsize = 24))

forceNetwork(Links = Links_temp, Nodes = Nodes_temp,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99,# Nodesize = 3,
             arrows = T, zoom = T, legend=T, charge = -20,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))

#ht
#htmlwidgets::saveWidget(as_widget(ht), "temporalGraph.html")
