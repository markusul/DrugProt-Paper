# load libraries
library(plotly)
library(networkD3)
library(webshot2)
library(htmlwidgets)
library(HiveR)
library(grid)
source("R/figureUtils.R")

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

# number of parameters per model, read from the projections (drug and lagged protein columns)
nPar <- sapply(c("6", "24", "48"), function(i){ load(paste0("Z/", i, ".RData")); ncol(Z) })

# total number of parameters in the model, over all three time points and all response proteins
sum(nPar) * nProt
P_Results <- c(P_Results, paste0("Total number of parameters in the overall model: ", 
                                 sum(nPar) * nProt))

# number of parameters in equation 1
P_Results <- c(P_Results, paste0("Number of parameters in equation 1: ", nPar["6"]))

# number of parameters in equation 2 and 3
P_Results <- c(P_Results, paste0("Number of parameters in equation 2 and 3: ", nPar["24"]))

# proteins of interest
load("results/anchor_opt/proteinSelection.RData")
P_selection <- which(prot_names_short %in% most_imp)
stopifnot(all(unique(path_s) %in% prot_names_short))
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

nSigDrug_oComb <- apply(allPvecsAdj[(nSingleDrugs+1):(nSingleDrugs + nComboDrugs), , ], 2, function(p) sum(p < alpha))
nDrug_oComb <- nComboDrugs * nProt
P_Results <- c(P_Results, "Number of significant drug effects (only combination) in the models over time")
P_Results <- c(P_Results, paste0("Time point 6h: ", nSigDrug_oComb[1]))
P_Results <- c(P_Results, paste0("Time point 24h: ", nSigDrug_oComb[2]))
P_Results <- c(P_Results, paste0("Time point 48h: ", nSigDrug_oComb[3]))
P_Results <- c(P_Results, "Percentage of significant drug combinations out of all interactions in the models over time")
P_Results <- c(P_Results, paste0("Time point 6h: ", nSigDrug_oComb[1]/nDrug_oComb))
P_Results <- c(P_Results, paste0("Time point 24h: ", nSigDrug_oComb[2]/nDrug_oComb))
P_Results <- c(P_Results, paste0("Time point 48h: ", nSigDrug_oComb[3]/nDrug_oComb))

# total effects of the drugs, Holm correction over all total effects
load("results/DrugEffects_total.RData")
allPvecsTotalAdj <- array(p.adjust(allPvecsTotal, method = "holm"), dim = dim(allPvecsTotal))
nSigTotal <- apply(allPvecsTotalAdj, 2, function(p) sum(p < alpha))
nSigTotal
P_Results <- c(P_Results, "Number of significant total drug effects over time (separate Holm correction)")
P_Results <- c(P_Results, paste0("Time point ", c(6, 24, 48), "h: ", nSigTotal))

# concentration of the significant single-drug effects at 6 h on few proteins
sig6 <- t(allPvecsAdj[!drugComb, 1, ] < alpha)       # proteins x single drugs
pos6 <- sapply(which(!drugComb), function(d){
  load(paste0("results/Coef/drugs/", d, ".RData"))
  dEff[, 1] > 0
})
nSigPerProt <- rowSums(sig6)
P_Results <- c(P_Results, paste0("Number of significant single-drug effects at 6h: ", sum(sig6)))
P_Results <- c(P_Results, paste0("Share of positive significant single-drug effects at 6h: ",
                                 round(mean(pos6[sig6]), 3)))
P_Results <- c(P_Results, paste0("Proteins with 0, 1, 2-5, 6-20, more than 20 significant single drugs at 6h: ",
                                 paste(table(cut(nSigPerProt, c(-1, 0, 1, 5, 20, Inf))), collapse = ", ")))
P_Results <- c(P_Results, paste0("Share of significant single-drug effects at 6h on proteins with more than 5 significant single drugs: ",
                                 round(sum(nSigPerProt[nSigPerProt > 5]) / sum(sig6), 3)))


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
P_Results <- c(P_Results, paste0("Significant single drugs on selected proteins: ", 
                                paste(names(pvec[!drugComb][pvec[!drugComb] < alpha]), 
                                      collapse = ", ")))


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
nonself <- Links_sum$source != Links_sum$target
nLinksFull <- nrow(unique(Links_sum[nonself, c("source", "target")]))
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
nonself <- Links_sum$source != Links_sum$target
nLinksSel <- nrow(unique(Links_sum[nonself, c("source", "target")]))
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

# number of p-values in Drug-Prot: protein edges at two transitions, group p-values of the direct
# drug effects at three time points, and of the total effects at 24 and 48 hours
nPval <- nProt * nProt * 2 + prod(dim(allPvecs)) + nTreatment * 2 * nProt
nPval
P_Results <- c(P_Results, paste0("Number of p-values: ", nPval))

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
##    edges between two selected proteins count for both of them
sel0 <- P_selection - 1   # Links_sum uses zero-based indices
SelToSel <- nrow(unique(Links_sum[nonself & Links_sum$source %in% sel0 & Links_sum$target %in% sel0,
                                  c("source", "target")]))
(nLinksSel + SelToSel) / nSelection
P_Results <- c(P_Results, paste0("Average degree of the SELECTED summary graph: ",
                                 (nLinksSel + SelToSel) / nSelection))

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
pMat[pMat == 2] <- NA

##### Figure 2: heatmap of the drug effects on the selected proteins #####
# pMat holds the corrected p-values, with NA for drug pairs without data
ht_df <- data.frame(x = factor(rep(colnames(pMat), each = nrow(pMat)), levels = colnames(pMat)),
                    y = factor(rep(rownames(pMat), times = ncol(pMat)), levels = rownames(pMat)),
                    p = as.vector(pMat))

ht_plot <- function(df) {
  ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = p)) +
    geom_tile(data = df[is.na(df$p), ], aes(colour = "No data"), fill = "#10bee0", linewidth = 0) +
    scale_fill_gradientn("P-value", limits = c(0, 1), breaks = c(0.05, 0.25, 0.5, 1), labels = c("0.05", "0.25", "0.5", "1"),
                         colours = c("#ffffff", "#d9d9d9", "#969696", "#737373", "#404040", "#1a1a1a", "#000000"),
                         values = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 1), na.value = "#10bee0",
                         guide = guide_colourbar(order = 1, frame.colour = "grey30", ticks.colour = "grey30",
                                                 barheight = unit(1.1, "in"), barwidth = unit(0.15, "in"))) +
    scale_colour_manual(NULL, values = c("No data" = "#10bee0"),
                        guide = guide_legend(order = 2, override.aes = list(fill = "#10bee0"))) +
    coord_fixed(expand = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_fig +
    theme(panel.grid = element_blank())
}

# (a) all 63 drugs, without drug names
ht_a <- ht_plot(ht_df) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")

# (b) the drugs that are part of a measured drug pair
in_pair <- rownames(pMat)[rowSums(!is.na(pMat)) > 1]
ht_b <- ht_plot(droplevels(ht_df[ht_df$x %in% in_pair & ht_df$y %in% in_pair, ])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

save_fig(list(ht_a, ht_b), 2, width = 7.5, height = 3.8, widths = c(0.42, 0.58), tags = c("(a)", "(b)"))

Nodes_sum$size <- 10
fN <- forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
                   Nodesize = "size", radiusCalculation = JS("Math.sqrt(d.nodesize)"),
                   fontSize = 0.1, bounded = T,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   Group = "group", opacity = 0.99, 
                   arrows = T, zoom = T, charge = -20,
                   opacityNoHover = TRUE, legend = F,
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
fN

saveWidget(fN, "figures/summary_sel.html")
webshot("figures/summary_sel.html", file = "figures/summary_sel.png", zoom = 3, 
        vwidth = 700,
        vheight = 700)

Nodes_temp$size <- 10
fN <- forceNetwork(Links = Links_temp, Nodes = Nodes_temp,
                   Nodesize = "size", radiusCalculation = JS("Math.sqrt(d.nodesize)"),
                   fontSize = 0.1, bounded = T,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   Group = "group", opacity = 0.99,# Nodesize = 3,
                   arrows = T, zoom = T, legend=F, charge = -10,
                   opacityNoHover = TRUE,
                   colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
fN

saveWidget(fN, "figures/temp_sel.html")
webshot("figures/temp_sel.html", file = "figures/temp_sel.png", zoom  = 3,
        vwidth = 700,
        vheight = 700)


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


##### Protein network of the selected proteins, Figure 3 #####
P_selection <- which(prot_names_short %in% most_imp)
stopifnot(all(unique(most_imp) %in% prot_names_short))

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


# label the selected proteins and the proteins connected to more than one of them
selected <- which(Nodes_sum$group == "Selected") - 1
neighbors <- table(unlist(lapply(selected, function(sel){
  unique(c(Links_sum[Links_sum$source == sel, "target"], 
           Links_sum[Links_sum$target == sel, "source"]))
})))
connectingNodes <- as.integer(names(neighbors)[neighbors > 1]) + 1
Nodes_sum[-c(connectingNodes, selected + 1), "name"] <- ""

connectedTemp <- which(unlist(lapply(Nodes_temp$name, 
                                      function(x)strsplit(x, "_")[[1]][1])) %in% 
                          Nodes_sum$name)
Nodes_temp[-connectedTemp, "name"] <- ""

# positions of the node labels
# each label is put at the candidate position around its node that covers no node and no other label
# and that crosses as little of the edges as possible, a label that is not next to its node gets a line to it
# nd and ed hold the nodes and edges in inches, the panel has the given width and height in inches
place_labels <- function(nd, ed, width, height) {
  lab <- which(nd$name != "")
  # half of the width and height of the text, with some space around it
  hw <- systemfonts::string_width(nd$name[lab], family = fig_family, size = 8) / 72 / 2 + 0.04
  hh <- 0.045 + 0.025
  # candidate positions, given by the direction from the node to the label in degrees and their distance
  cand <- expand.grid(angle = seq(0, 345, by = 15), d = c(0.06, 0.1, 0.16, 0.24, 0.34, 0.46, 0.6, 0.8))
  ux <- cos(cand$angle * pi / 180)
  uy <- sin(cand$angle * pi / 180)
  # only a label at the shortest distance is drawn without a line to its node
  has_line <- cand$d > 0.06

  # length of the part of a segment (x, y, xend, yend) that lies inside a box (xmin, xmax, ymin, ymax)
  # for all pairs of a box and a segment
  inside <- function(seg, box) {
    n <- nrow(box)
    m <- nrow(seg)
    if (n == 0 || m == 0) return(matrix(0, n, m))
    b <- box[rep(seq_len(n), times = m), , drop = FALSE]
    s <- seg[rep(seq_len(m), each = n), , drop = FALSE]
    dx <- s[, 3] - s[, 1]
    dy <- s[, 4] - s[, 2]
    p <- cbind(-dx, dx, -dy, dy)
    q <- cbind(s[, 1] - b[, 1], b[, 2] - s[, 1], s[, 2] - b[, 3], b[, 4] - s[, 2])
    t0 <- rep(0, n * m)
    t1 <- rep(1, n * m)
    for (k in 1:4) {
      r <- q[, k] / p[, k]
      t0 <- ifelse(p[, k] < 0, pmax(t0, r), t0)
      t1 <- ifelse(p[, k] > 0, pmin(t1, r), t1)
      t1 <- ifelse(p[, k] == 0 & q[, k] < 0, -1, t1)
    }
    matrix(pmax(t1 - t0, 0) * sqrt(dx^2 + dy^2), n, m)
  }

  # do two segments cross, for all pairs of a segment in a and a segment in b
  # the segments in a all start at the same point
  crosses <- function(a, b) {
    left <- function(px, py, qx, qy, rx, ry) (ry - py) * (qx - px) > (qy - py) * (rx - px)
    one <- function(v) matrix(v, nrow(a), nrow(b), byrow = TRUE)
    start <- left(a[1, 1], a[1, 2], b[, 1], b[, 2], b[, 3], b[, 4])
    end <- left(a[, 3], a[, 4], one(b[, 1]), one(b[, 2]), one(b[, 3]), one(b[, 4]))
    one(start) != end &
      left(a[, 1], a[, 2], a[, 3], a[, 4], one(b[, 1]), one(b[, 2])) !=
      left(a[, 1], a[, 2], a[, 3], a[, 4], one(b[, 3]), one(b[, 4]))
  }

  edges <- as.matrix(ed)
  # places the labels one after the other in the given order, and then twice again given all other labels
  # returns the boxes of the labels, their lines to the nodes, and the total cost
  place <- function(first) {
    box <- matrix(NA_real_, length(lab), 4)
    line <- matrix(NA_real_, length(lab), 4)
    total <- numeric(length(lab))
    for (i in rep(first, 3)) {
      node <- c(nd$x[lab[i]], nd$y[lab[i]])
      # boxes of all candidates, the line of a candidate ends at the point of its box closest to the node
      r <- cand$d + pmin(hw[i] / abs(ux), hh / abs(uy))
      b <- cbind(node[1] + r * ux - hw[i], node[1] + r * ux + hw[i], node[2] + r * uy - hh, node[2] + r * uy + hh)
      l <- cbind(node[1], node[2], pmin(pmax(node[1], b[, 1]), b[, 2]), pmin(pmax(node[2], b[, 3]), b[, 4]))
      # short distances and labels without a line are preferred, and the label stays inside the panel
      cost <- 12 * cand$d + 30 * pmax(cand$d - 0.3, 0) + 6 * has_line +
        1e4 * (b[, 1] < 0.05 | b[, 2] > width - 0.05 | b[, 3] < 0.05 | b[, 4] > height - 0.05)
      # the label of a node without edges goes above the node if there is space
      far <- sqrt(pmin((edges[, 1] - node[1])^2 + (edges[, 2] - node[2])^2,
                       (edges[, 3] - node[1])^2 + (edges[, 4] - node[2])^2)) > 0.06
      if (all(far)) cost <- cost + 5 * (cand$angle != 90)
      # the label covers no node
      dist <- sqrt(pmax(outer(b[, 1], nd$x, "-"), 0, -outer(b[, 2], nd$x, "-"))^2 +
                     pmax(outer(b[, 3], nd$y, "-"), 0, -outer(b[, 4], nd$y, "-"))^2)
      cost <- cost + 300 * rowSums(dist < 0.035)
      # no other node is close to the label, so that it is clear which node it belongs to
      # this matters most for a label without a line
      cost <- cost + rowSums(dist[, -lab[i], drop = FALSE] < ifelse(has_line, 0.08, 0.2)) * ifelse(has_line, 40, 100)
      # the label crosses as little of the edges as possible
      cost <- cost + 60 * rowSums(inside(edges, b))
      # the label does not overlap the other labels and their lines
      others <- setdiff(which(!is.na(box[, 1])), i)
      ob <- box[others, , drop = FALSE]
      ol <- line[others, , drop = FALSE]
      ol <- ol[!is.na(ol[, 1]), , drop = FALSE]
      cost <- cost + 1e3 * rowSums(outer(b[, 1], ob[, 2], "<") & outer(b[, 2], ob[, 1], ">") &
                                     outer(b[, 3], ob[, 4], "<") & outer(b[, 4], ob[, 3], ">"))
      cost <- cost + 300 * rowSums(inside(ol, b) > 0)
      # the line of the label passes no other node, no other label, and no edge of another node
      vx <- l[, 3] - node[1]
      vy <- l[, 4] - node[2]
      t <- pmin(pmax((outer(vx, nd$x - node[1]) + outer(vy, nd$y - node[2])) / (vx^2 + vy^2), 0), 1)
      dist <- sqrt((outer(rep(1, nrow(cand)), nd$x - node[1]) - t * vx)^2 +
                     (outer(rep(1, nrow(cand)), nd$y - node[2]) - t * vy)^2)
      cost <- cost + 150 * has_line * rowSums(dist[, -lab[i], drop = FALSE] < 0.06)
      cost <- cost + 300 * has_line * colSums(inside(l, ob) > 0)
      cost <- cost + 40 * has_line * rowSums(crosses(l, edges[far, , drop = FALSE]))

      best <- which.min(cost)
      box[i, ] <- b[best, ]
      line[i, ] <- if (has_line[best]) l[best, ] else NA
      total[i] <- cost[best]
    }
    list(box = box, line = line, total = sum(total))
  }

  # the result depends on the order of the labels, the best of six orders is used
  # the orders start in a crowded region, in an empty region, or on one of the four sides
  crowded <- sapply(lab, function(v) sum(sqrt((nd$x - nd$x[v])^2 + (nd$y - nd$y[v])^2) < 0.6))
  res <- lapply(list(-crowded, crowded, nd$x[lab], -nd$x[lab], nd$y[lab], -nd$y[lab]), function(v) place(order(v)))
  res <- res[[which.min(sapply(res, function(r) r$total))]]
  data.frame(name = nd$name[lab], lx = (res$box[, 1] + res$box[, 2]) / 2, ly = (res$box[, 3] + res$box[, 4]) / 2,
             x = res$line[, 1], y = res$line[, 2], xend = res$line[, 3], yend = res$line[, 4])
}

# static version of the two graphs, with a fixed layout so that the figure is reproducible
# width and height are the size of the panel in inches, they are needed to place the labels
net_plot <- function(Links, Nodes, colors, seed, width, height) {
  g <- igraph::graph_from_data_frame(data.frame(from = Links$source + 1, to = Links$target + 1),
                                     directed = TRUE, vertices = data.frame(id = seq_len(nrow(Nodes))))
  g <- igraph::simplify(g)
  isolated <- igraph::degree(g) == 0
  set.seed(seed)
  xy <- matrix(NA_real_, nrow(Nodes), 2)
  xy[!isolated, ] <- igraph::layout_components(igraph::induced_subgraph(g, which(!isolated)),
                                               layout = igraph::layout_with_fr)
  xy[!isolated, ] <- apply(xy[!isolated, , drop = FALSE], 2, function(v) (v - min(v)) / diff(range(v)))
  # nodes without an edge are placed in a row below the graph
  if (any(isolated)) xy[isolated, ] <- cbind(seq(0.08, 0.92, length.out = sum(isolated)), -0.14)
  # the coordinates are in inches from here on, with a margin around the graph
  margin <- 0.25
  xy[, 1] <- margin + (xy[, 1] - min(xy[, 1])) / diff(range(xy[, 1])) * (width - 2 * margin)
  xy[, 2] <- margin + (xy[, 2] - min(xy[, 2])) / diff(range(xy[, 2])) * (height - 2 * margin)
  nd <- data.frame(x = xy[, 1], y = xy[, 2], name = Nodes$name, group = factor(Nodes$group, levels = names(colors)))
  # edges end shortly before the node they point to, self-dependencies are not drawn in the summary graph
  ed <- data.frame(x = nd$x[Links$source + 1], y = nd$y[Links$source + 1],
                   xend = nd$x[Links$target + 1], yend = nd$y[Links$target + 1])
  ed <- ed[ed$x != ed$xend | ed$y != ed$yend, ]
  len <- sqrt((ed$xend - ed$x)^2 + (ed$yend - ed$y)^2)
  gap <- pmin(0.04 / len, 0.4)
  ed$xend <- ed$xend - gap * (ed$xend - ed$x)
  ed$yend <- ed$yend - gap * (ed$yend - ed$y)
  labels <- place_labels(nd, ed, width, height)
  # white copies of the labels, shifted in all directions, separate the text from the edges
  halo <- do.call(rbind, lapply(seq(0, 2 * pi, length.out = 17)[-1], function(a) {
    transform(labels, lx = lx + 0.012 * cos(a), ly = ly + 0.012 * sin(a))
  }))

  ggplot(nd, aes(x, y)) +
    geom_segment(data = ed, aes(xend = xend, yend = yend), colour = "grey70", linewidth = 0.3,
                 arrow = arrow(length = unit(0.04, "in"), type = "closed")) +
    geom_segment(data = labels[!is.na(labels$x), ], aes(xend = xend, yend = yend), linewidth = 0.25) +
    geom_text(data = halo, aes(lx, ly, label = name), family = fig_family, size = 8 / .pt, colour = "white") +
    geom_point(aes(colour = group), size = 1.3) +
    geom_text(data = labels, aes(lx, ly, label = name), family = fig_family, size = 8 / .pt) +
    scale_colour_manual(NULL, values = colors, drop = FALSE) +
    coord_cartesian(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
    theme_void(base_size = 10, base_family = fig_family) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

# size of the figure in inches and the share of the two panels in its width
# the tags and the legend take about 0.55 in of the height
fig3_size <- c(7.5, 4.8)
fig3_share <- c(0.42, 0.58)
net_a <- net_plot(Links_sum, Nodes_sum, c(Selected = fig_colors[["orange"]], Connected = fig_colors[["blue"]]),
                  seed = 1, width = fig3_share[1] * fig3_size[1], height = fig3_size[2] - 0.55)
net_b <- net_plot(Links_temp, Nodes_temp,
                  c("6h" = fig_colors[["blue"]], "24h" = fig_colors[["orange"]], "48h" = fig_colors[["green"]]),
                  seed = 1, width = fig3_share[2] * fig3_size[1], height = fig3_size[2] - 0.55)

save_fig(list(net_a, net_b), 3, width = fig3_size[1], height = fig3_size[2], widths = fig3_share, tags = c("(a)", "(b)"))



##### visualization of empirical validation of p-values #####

load("results/PvalAnalysis.RData")

library(ggplot2)

expTimes <- c(6, 24, 48)
alpha <- 0.05

df <- do.call(rbind, lapply(seq_along(res), function(t){
  m <- res[[t]]
  do.call(rbind, lapply(seq_len(nrow(m)), function(k){
    d <- m[[k, "coef"]]
    d$tp <- expTimes[t]
    d$protein <- k
    d
  }))
}))
df <- df[is.finite(df$rate), ]

df$metric <- factor(ifelse(df$true_effect, "Power", "Type I error"),
                    levels = c("Type I error", "Power"))
df$tp <- factor(df$tp, levels = expTimes, labels = paste(expTimes, "h"))
df$group <- factor(df$group, levels = c("single", "double", "protein"),
                   labels = c("Single drug", "Drug pair", "Lagged protein"))

alpha_line <- data.frame(metric = factor("Type I error", levels = levels(df$metric)),
                         y = alpha)

p <- ggplot(df, aes(x = group, y = rate, fill = group)) +
  geom_boxplot(linewidth = 0.1, outlier.size = 0.02) +
  #geom_violin(linewidth = 0.6) + 
  #geom_point() +
  geom_hline(data = alpha_line, aes(yintercept = y),
             linetype = "dashed", linewidth = 0.8, colour = "grey30") +
  facet_grid(metric ~ tp, scales = "free", space = "free_x") +
  scale_fill_manual(values = c("Single drug" = "#0072B2",
                               "Drug pair" = "#E69F00",
                               "Lagged protein" = "#009E73")) +
  labs(x = "Effect type", y = "Rejection rate") +
  theme_fig +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

p
save_fig(p, 4, width = 7.5, height = 5.6)