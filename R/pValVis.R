library(plotly)
library(networkD3)

# load ordering of drugs (sort drugs with experiments together)
load('data/order.RData')
nDrugs <- length(drugOrder)

load("results/impProt.RData")
P <- names(important_prot[1:20])
P <- sapply(P, function(p) which(prot_names_short == strsplit(p, '[.]')[[1]][1]))


##### Drug Effects ####
load("results/DrugEffects.RData")

P <- c(1:nrow(allPvecs))
#P <- 1387

# collect min p value of drug effect over proteins
pvec <- apply(matrix(allPvecs[P, ], nrow = length(P)), 2, min)
pvec <- pmin(pvec * length(P), 1)
names(pvec) <- colnames(allPvecs)

pMat <- matrix(NA, nrow = nDrugs, ncol = nDrugs)
rownames(pMat) <- colnames(pMat) <- names(pvec)[1:nDrugs]
for(l in names(pvec)){
  drugs <- strsplit(l, ":")[[1]]
  if(length(drugs) == 1) drugs <- c(drugs, drugs)
  pMat[drugs[1], drugs[2]] <- pvec[l]
  pMat[drugs[2], drugs[1]] <- pvec[l]
}

pMat <- as.matrix(pMat)
pMat[is.na(pMat)] <- 2
pMat <- pMat[drugOrder, drugOrder]

ht <- plot_ly(z = pMat, x = colnames(pMat), y = colnames(pMat), 
              type = "heatmap", colors = "Greys") %>%
              layout(title = prot_names_short[P])
ht

#htmlwidgets::saveWidget(as_widget(ht), paste0(prot_names_short[which(prot_names == P)], ".html"))


##### Protein Network ####

load("results/proteinNetwork.RData")

# proteins of interest
#P <- 200:260

# significance level for full protein network
alpha <- 0.05 / length(P) / 3

# transform p value to links using alpha
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

# summary graph
Links_sum <- do.call(rbind, Links_all)
Links_sum$source <- Links_sum$source - 1
Links_sum$target <- Links_sum$target - 1
Links_sum$value <- 1
Nodes_sum <- data.frame(name = prot_names_short, group = 1, size = 1)

P.ind <- P - 1
rel.Links <- Links_sum$source %in% P.ind | Links_sum$target %in% P.ind
Links_sum <- Links_sum[rel.Links, ]
rel.Nodes <- sort(unique(unlist(Links_sum[, c('source', 'target')])))
Nodes_sum$group[P] <- 2
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
             arrows = T, zoom = T, charge = -5,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
fN
#htmlwidgets::saveWidget(as_widget(fN), "SummaryGraph.html")


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
nodes$color <- "black"
#nodes$color <- c("black", "red", "blue", "green", "violet", "yellow", "darkgreen")[clusters]
#rep(c("black", "red", "blue", "green", "violet", "yellow", "darkgreen"), times = as.numeric(table(clusters)))


HEC <- list()
HEC$nodes <- nodes
HEC$edges <- edges
HEC$type <- "2D"
HEC$desc <- "HairEyeColor data set"
HEC$axis.cols <- c("grey", "grey")
class(HEC) <- "HivePlotData"

chkHPD(HEC) # answer of FALSE means there are no problems

plotHive(HEC, ch = 0.001, bkgnd = "white", 
         axLabs = c("24h", "6h", "48h"), 
         axLab.gpar = gpar(col = "black", fontsize = 24))

forceNetwork(Links = Links_temp, Nodes = Nodes_temp,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99,# Nodesize = 3,
             arrows = T, zoom = T, legend=T, charge = -5,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))

#ht
#htmlwidgets::saveWidget(as_widget(ht), "temporalGraph.html")
