# load libraries
library(plotly)
library(networkD3)
library(webshot2)
library(htmlwidgets)
library(HiveR)
library(grid)

# significance level for tests
alpha <- 0.05

# replace drug ids by names
load("data/drugLookup.RData")
replace_drug_ids <- function(x) {
  ids <- strsplit(x, "_|`")[[1]]
  ids <- ids[grepl("#", ids)]
  names <- sapply(ids, function(id) drug_lookup[[id]])
  paste(names, collapse = ":")
}

# load ordering of drugs (sort drugs with experiments together)
load('data/order.RData')
nDrugs <- length(drugOrder)
drugOrder <- sapply(drugOrder, replace_drug_ids)
nProt <- length(prot_names_short)

# load data for Drug Effects
load("results/DrugEffects.RData")

# load p-value network
load("results/proteinNetworkPval.RData")

# collect results
P_Results <- c("Results from the P-Value analysis of DrugProt", "-----------------------------------")
P_Results <- c(P_Results, paste0("Number of proteins: ", nProt))
P_Results <- c(P_Results, paste0("Number of treatments: ", length(treatment)))

drugComb <- sapply(treatment, function(x) grepl(":", x))
# number of single drugs
nSingleDrugs <- sum(!drugComb)
# number of drug combinations
nComboDrugs <- sum(drugComb)

P_Results <- c(P_Results, paste0("Number of single drugs: ", nSingleDrugs))
P_Results <- c(P_Results, paste0("Number of drug combinations: ", nComboDrugs))

# number of datapoints in the models
P_Results <- c(P_Results, "Number of data points in the models:")
for(i in c("6", "24", "48")){
  load(paste0("Z/", i, ".RData"))
  P_Results <- c(P_Results, paste0("Time point ", i, "h:", nrow(Z)))
  print(dim(Z))
}

# total number of parameters in the model
alpha_n <- nSingleDrugs * 2 * 3
beta_n <- nComboDrugs * 2 * 3
gamma_n <- nProt * 2
(alpha_n + beta_n + gamma_n) * nProt #times the number of models
P_Results <- c(P_Results, paste0("Total number of parameters in the overall model: ", 
                                   (alpha_n + beta_n + gamma_n) * nProt))

# number of parameters in equation 1
nSingleDrugs * 2 + nComboDrugs * 2
P_Results <- c(P_Results, paste0("Number of parameters in equation 1: ", 
                                   nSingleDrugs * 2 + nComboDrugs * 2))

# number of parameters in equation 2 and 3
nSingleDrugs * 2 + nComboDrugs * 2 + nProt
P_Results <- c(P_Results, paste0("Number of parameters in equation 2 and 3: ", 
                                  nSingleDrugs * 2 + nComboDrugs * 2 + nProt))

# proteins of interest
load("results/anchor_opt/proteinSelection.RData")
P_selection <- which(prot_names_short %in% path_s)
nSelection <- length(P_selection)
P_Results <- c(P_Results, paste0("Number of selected proteins by AnchorForest: ", nSelection))

print(prot_names_short[P_selection])
print(nSelection)

##### P-value analysis of drug effects #####
P_Results <- c(P_Results, "", "P-Value analysis of drug effects", "-----------------------------------")

# total number of potential drug effects
totalEffects <- length(treatment) * 3 * nProt
totalEffects
prod(dim(allPvecs))
P_Results <- c(P_Results, paste0("Total number of potential drug effects in the models: ", 
                                  totalEffects))

# total significant drug effects
allPvecsAdj <- array(p.adjust(allPvecs, method = "holm"), dim = dim(allPvecs))
nSigDrug <- apply(allPvecsAdj, 2, function(p) sum(p < alpha))
nSigDrug
P_Results <- c(P_Results, "Number of significant drug effects in the models over time")
P_Results <- c(P_Results, paste0("Time point 6h: ", nSigDrug[1]))
P_Results <- c(P_Results, paste0("Time point 24h: ", nSigDrug[2]))
P_Results <- c(P_Results, paste0("Time point 48h: ", nSigDrug[3]))

#adjust p values of drug effects on selected proteins
selPvecs <- allPvecs[, , P_selection]
selPvecs <- array(p.adjust(selPvecs, method = "holm"), dim = dim(selPvecs))

# collect min p value of drug effect over proteins and time points
pvec <- apply(selPvecs, 1, function(p) min(p))
names(pvec) <- sapply(treatment, replace_drug_ids)

# number of significant drug effects on selected proteins
sum(pvec < alpha)
P_Results <- c(P_Results, paste0("Number of significant drug effects on selected proteins: ", sum(pvec < alpha)))

# significant single drugs on selected proteins
sum(pvec[!drugComb] < alpha)
P_Results <- c(P_Results, paste0("Number of significant single drugs on selected proteins: ", 
                                  sum(pvec[!drugComb] < alpha)))

# significant drug combinations on selected proteins
sum(pvec[drugComb] < alpha)
P_Results <- c(P_Results, paste0("Number of significant drug combinations on selected proteins: ", 
                                  sum(pvec[drugComb] < alpha)))

##### P-value network analysis #####
P_Results <- c(P_Results, "", "P-Value network analysis", "-----------------------------------")

# analyze edges in the full graph
# p-value correction
Pval_full <- Pval_all
Tpval <- sapply(Pval_all, function(links) links$pvalue)
Tpval.corr <- matrix(p.adjust(Tpval, method = "BH"), ncol = 2)
Pval_full[[1]][, "pvalue"] <- Tpval.corr[, 1]
Pval_full[[2]][, "pvalue"] <- Tpval.corr[, 2]

# convert p-values to links
Links_all <- lapply(Pval_full, function(links) links[links$pvalue < alpha, ])

# Prepare Summary Graph
Links_sum <- do.call(rbind, Links_all)
Links_sum$source <- Links_sum$source - 1
Links_sum$target <- Links_sum$target - 1
Links_sum$value <- 1

# number of links in the summary graph
nLinksFull <- nrow(unique(Links_sum[, -3])) - nProt
nLinksFull
P_Results <- c(P_Results, paste0("Number of links in the FULL summary graph: ", 
                                  nLinksFull))

# number of links in the temporal graph
nrow(Links_sum)
P_Results <- c(P_Results, paste0("Number of links in the FULL temporal graph: ", 
                                  nrow(Links_sum)))

# select relevant p-values
Pval_sel <- lapply(Pval_all, function(links){
  rel.Links <- links$source %in% P_selection | links$target %in% P_selection
  links[rel.Links, ]
})

# p-value correction
Tpval <- sapply(Pval_sel, function(links) links$pvalue)
Tpval.corr <- matrix(p.adjust(Tpval, method = "BH"), ncol = 2)
Pval_sel[[1]][, "pvalue"] <- Tpval.corr[, 1]
Pval_sel[[2]][, "pvalue"] <- Tpval.corr[, 2]

# convert p-values to links
Links_all <- lapply(Pval_sel, function(links) links[links$pvalue < alpha, ])

# Prepare Summary Graph
Links_sum <- do.call(rbind, Links_all)
Links_sum$source <- Links_sum$source - 1
Links_sum$target <- Links_sum$target - 1
Links_sum$value <- 1

# number of links in the summary graph
nLinksSel <- nrow(unique(Links_sum[, -3])) - length(P_selection)
nLinksSel
P_Results <- c(P_Results, paste0("Number of links in the SELECTED summary graph: ", 
                                  nLinksSel))

# number of links in the temporal graph
nrow(Links_sum)
P_Results <- c(P_Results, paste0("Number of links in the SELECTED temporal graph: ", 
                                  nrow(Links_sum)))

P_Results <- c(P_Results, "")

# number of possible edges in the full summary graph
nProt*(nProt-1)
P_Results <- c(P_Results, paste0("Number of possible edges in the FULL summary graph: ", 
                                  nProt*(nProt-1)))

# number of possible edges in the full temporal graph
nProt * nProt * 2
P_Results <- c(P_Results, paste0("Number of possible edges in the FULL temporal graph: ", 
                                  nProt * nProt * 2))

# sanity check
nProt * nProt * 2 + prod(dim(allPvecs)) * 2 == (alpha_n + beta_n + gamma_n) * nProt

# number of possible edges in the selected summary graph
((nProt-nSelection) * nSelection + nSelection * (nSelection-1)/2) * 2
P_Results <- c(P_Results, paste0("Number of possible edges in the SELECTED summary graph: ", 
                                  ((nProt-nSelection) * nSelection + nSelection * (nSelection-1)/2) * 2))

# number of possible edges in the selected temporal graph
nSelection * nProt * 2 + nSelection * (nProt - nSelection) * 2
P_Results <- c(P_Results, paste0("Number of possible edges in the SELECTED temporal graph: ", 
                                  nSelection * nProt * 2 + nSelection * (nProt - nSelection) * 2))

# degree of summary graphs
P_Results <- c(P_Results, "")

# full graph
nLinksFull / nProt * 2 # every edge affects two nodes
P_Results <- c(P_Results, paste0("Average degree of the FULL summary graph: ", 
                                  nLinksFull / nProt * 2))

# selected graph
SelToSel <- sum(apply(Links_sum,1, function(links){
  links["source"] %in% P_selection & links["target"] %in% P_selection
}))

(nLinksSel - SelToSel) / nSelection + SelToSel / nSelection * 2 # only some edges affect two nodes
P_Results <- c(P_Results, paste0("Average degree of the SELECTED summary graph: ", 
                                  (nLinksSel - SelToSel) / nSelection + SelToSel / nSelection * 2))

#selected nodes + connected nodes
rel.Nodes <- sort(unique(c(P_selection-1, unlist(Links_sum[, c('source', 'target')]))))
Nodes_sum <- data.frame(name = prot_names_short[rel.Nodes+1], group = "Connected", size = 1)
Nodes_sum$group[rel.Nodes %in% (P_selection-1)] <- "Selected"

#reorganize link index
for(i in 1:length(rel.Nodes)){
  Links_sum[, c('source', 'target')][Links_sum[, c('source', 'target')] == rel.Nodes[i]] <- i - 1
}

# Prepare Temporal Graph
expTimes <- c(6, 24, 48)
rel6 <- sort(unique(c(P_selection, Links_all[[1]][, "source"])))
rel24 <- sort(unique(c(P_selection, Links_all[[1]][, "target"], Links_all[[2]][, "source"])))
rel48 <- sort(unique(c(P_selection, Links_all[[2]][, "target"])))
rel <- list(rel6, rel24, rel48)
lenRel <- c(0, length(rel6), length(rel24), length(rel48))

nodenames <- c(paste(prot_names_short[rel6], expTimes[1], sep = '_'), 
               paste(prot_names_short[rel24], expTimes[2], sep = '_'), 
               paste(prot_names_short[rel48], expTimes[3], sep = '_'))
nodegroups <- rep(paste0(expTimes, "h"), times = c(length(rel6), length(rel24), length(rel48)))

# reorganize link index to match new nodenames
Links_temp <- lapply(1:2, function(t){
  links <- Links_all[[t]]
  for(i in 1:length(rel[[t]])){
    links[links[, 1] == rel[[t]][i], 1] <- i - 1 + sum(lenRel[1:t])
  }
  for(i in 1:length(rel[[t+1]])){
    links[links[, 2] == rel[[t+1]][i], 2] <- i - 1 + sum(lenRel[1:(t+1)])
  }
  links
})
Links_temp <- do.call(rbind, Links_temp)
Links_temp$value <- 1
Nodes_temp <- data.frame(name = nodenames, group = nodegroups, size = 0.1)
Nodes_temp$radius <- as.numeric(c(rel6, rel24, rel48))

# Prepare for HivePlot
edges <- Links_temp
edges[, 1:2] <- edges[, 1:2] + 1

names(edges) <- c("id1", "id2", "weight")
# color
row.names(edges) <- NULL
edges$id1 <- as.integer(edges$id1)
edges$id2 <- as.integer(edges$id2)
edges$color <- "black"
edges$weight <- 0.1

nodes <- Nodes_temp
names(nodes) <- c("lab", "axis", "size", "radius")

nodes$axis[nodes$axis == "6h"] <- 2
nodes$axis[nodes$axis == "24h"] <- 1
nodes$axis[nodes$axis == "48h"] <- 3
nodes$axis <- as.integer(nodes$axis)
nodes$id <- 1:nrow(nodes)
nodes$radius <- nodes$radius * 3
nodes$color <- "grey"
nodes$color[unlist(rel) %in% P_selection] <- "black"

HEC <- list()
HEC$nodes <- nodes
HEC$edges <- edges
HEC$type <- "2D"
HEC$desc <- "HairEyeColor data set"
HEC$axis.cols <- c("grey", "grey")
class(HEC) <- "HivePlotData" 

# prepare p value of drug effects matrix for heatmap
pMat <- matrix(NA, nrow = nDrugs, ncol = nDrugs)
rownames(pMat) <- colnames(pMat) <- names(pvec)[1:nDrugs]
for(l in names(pvec)){
  drugs <- strsplit(l, ":")[[1]]
  if(length(drugs) == 1) drugs <- c(drugs, drugs)
  pMat[drugs[1], drugs[2]] <- pvec[l]
  pMat[drugs[2], drugs[1]] <- pvec[l]
}

pMat <- as.matrix(pMat)
pMat[is.na(pMat)] <- 2 # set no data to 2
pMat <- pMat[drugOrder, drugOrder]


##### Visualizations #####
# plot heatmap of drug effects
ht <- plot_ly(z = pMat, x = colnames(pMat), y = colnames(pMat), 
              type = "heatmap", colors = "Greys",
              colorbar = list(title = "<b> P-Values </b>"))
ht

#save heatmap
saveWidget(as_widget(ht), "figures/ht_sel.html")
webshot("figures/ht_sel.html", file = "figures/ht_sel.png")

# zoom in on significant part of heatmap
ht_zoom <- ht
ht_zoom <- ht_zoom %>% layout(xaxis = list(range = c(0, 11.5)), 
                              yaxis = list(range = c(0, 11.5))) %>%
  layout(xaxis = list(tickfont = list(size = 20)), 
         yaxis = list(tickfont = list(size = 20)))
ht_zoom
saveWidget(as_widget(ht_zoom), "figures/ht_zoom_sel.html")
webshot("figures/ht_zoom_sel.html", file = "figures/ht_zoom_sel.png")

Nodes_sum$size <- 10
fN <- forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
                   Nodesize = "size", radiusCalculation = JS("Math.sqrt(d.nodesize)"),
                   fontSize = 0.1, bounded = T,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   Group = "group", opacity = 0.99, 
                   arrows = T, zoom = T, charge = -20,
                   opacityNoHover = TRUE, legend = T,
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
fN

saveWidget(fN, "figures/summary_sel.html")
webshot("figures/summary_sel.html", file = "figures/summary_sel.png", zoom = 1)

Nodes_temp$size <- 10
fN <- forceNetwork(Links = Links_temp, Nodes = Nodes_temp,
                   Nodesize = "size", radiusCalculation = JS("Math.sqrt(d.nodesize)"),
                   fontSize = 0.1, bounded = T,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   Group = "group", opacity = 0.99,# Nodesize = 3,
                   arrows = T, zoom = T, legend=T, charge = -10,
                   opacityNoHover = TRUE,
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
fN

saveWidget(fN, "figures/temp_sel.html")
webshot("figures/temp_sel.html", file = "figures/temp_sel.png")


plotHive(HEC, ch = 0.001, bkgnd = "white", 
         axLabs = c("24h", "6h", "48h"), 
         axLab.gpar = gpar(col = "black", fontsize = 24))

png(filename="figures/hive_sel.png", width = 1000, height = 700)
plotHive(HEC, ch = 0.001, bkgnd = "white", 
         axLabs = c("24h", "6h", "48h"), 
         axLab.gpar = gpar(col = "black", fontsize = 24))
dev.off()


# file to save results to
fileConn <- file("results/P_Results.txt")
writeLines(P_Results, fileConn)
close(fileConn)
