library(shinydashboard)
library(plotly)
library(networkD3)
library(HiveR)
library(grid)
library(readxl)
library(shinyWidgets)

#choises of timepoints for drug effects
t_choice <- c("6h", "24h", "48h")

# replace drug ids by names
load('data/drugLookup.RData')
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

alpha <- 0.05

load("results/proteinNetworkPval.RData")
Pval_full <- do.call(rbind, Pval_all)
Pval_full$pvalue <- p.adjust(Pval_full$pvalue, "BH")
Links_full <- Pval_full[Pval_full$pvalue < alpha, ]
dim(Links_full)

library(igraph)

Links_g <- Links_full[, c("source", "target")]
Links_g


# create I graph from Links
g <- graph_from_data_frame(Links_full[, c("source", "target")], 
                           directed = TRUE)
summary(g)



plot(g)

a <- which(prot_names_short %in% "LMNA")
b <- which(prot_names_short %in% "EMD")


k_shortest_paths(g, from = a, to = b, k = 3)

?k_shortest_paths

P_selection <- which(prot_names_short %in% c("LMNA", "CTNA1"))

# select relevant p-values
Pval_sel <- lapply(t_choice[-1], function(tp){
  load(paste0("results/proteinNetworkPval_", tp, ".RData"))
  links <- Pval_all
  rel.Links <- links$source %in% P_selection | links$target %in% P_selection
  print(mem_used())
  links[rel.Links, ]
})
print(mem_used())
# p-value correction
Tpval <- sapply(Pval_sel, function(links) links$pvalue)
Tpval.corr <- matrix(p.adjust(Tpval, method = "BH"), ncol = 2)
Pval_sel[[1]][, "pvalue"] <- Tpval.corr[, 1]
Pval_sel[[2]][, "pvalue"] <- Tpval.corr[, 2]

# convert p-values to links
Links_all <- lapply(Pval_sel, function(links) links[links$pvalue < alpha, ])


Links_full <- do.call(rbind, Links_all)

Links_g <- Links_full[, c("source", "target")]
Links_g$source <- prot_names_short[Links_g$source]
Links_g$target <- prot_names_short[Links_g$target]

g <- graph_from_data_frame(Links_g, 
                           directed = TRUE)
summary(g)

plot(g)
which(V(g) == a)

a <- which(prot_names_short %in% "LMNA")
b <- which(prot_names_short %in% "CTNA1")

from <- which(names(V(g)) == a)
to <-  which(names(V(g)) == b)

sPaths <- k_shortest_paths(g, from = "LMNA", to = "CTNA1", k = 3, mode = "all")
str(sPaths$epaths[[1]])
