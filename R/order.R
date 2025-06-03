library(hdi)
library(plotly)
library(parallel)
library(data.table)

load("data/prepData.RData")
load("data/protNames.RData")

data$pertLabel[data$pertLabel == 'no no'] <- 'no'

# remove proteins with less than 11 unique values
abundand <- apply(data[, prot_names], 2, function(x) length(unique(x)) > 10)
prot_names <- prot_names[abundand]

#mean per protein plate at time 0
dat0 <- data[data$pert_time == 0, ]
dat0 <- aggregate(dat0[, prot_names], by = list(dat0$protein_plate), FUN = mean)

# sort drugs by information
drugOrder <- names(sort(colSums(data[, pert_names] != 0), decreasing = T))
drugOrder <- paste0('`', drugOrder, '`')

# collect protein names in short
prot_names_short <- sapply(prot_names, function(p) strsplit(p, '[.]')[[1]][1])


save(file = 'data/order.RData', prot_names_short, drugOrder, prot_names)
