library(hdi)
library(plotly)

load("data/prepData.RData")
load("data/protNames.RData")
dim(data)

#mean per protein plate at time 0
dat0 <- data[data$pert_time == 0, ]
dat0 <- aggregate(dat0[, prot_names], by = list(dat0$protein_plate), FUN = mean)

# sort drugs by information
drugOrder <- names(sort(colSums(data[, pert_names] != 0), decreasing = T))
drugOrder <- paste0('`', drugOrder, '`')

# add mean at 0 to data
names(dat0)[-1] <- paste0(names(dat0)[-1], "_0")
row.names(dat0) <- dat0[, 1]
dat0 <- dat0[data$protein_plate, -1]

# add mean at 0 to data
dat <- data
dat <- cbind(dat, dat0)

# remove M453(ATCC) (no baseline protein expression)
dat <- dat[dat$protein_plate != "M453(ATCC)", ]

# add nois to imputation for computational und statistical stability
add_noise_to_imputation <- function(x){
  border <- sort(unique(x))[2]
  imp <- which(x < border)
  x_ <- x[-imp]
  y <- x
  y[imp] <- truncnorm::rtruncnorm(length(imp), a = -Inf, b = border, mean = mean(x_), sd = sd(x_))
  y
}

datI <- dat
plate_idx <- lapply(unique(dat$protein_plate), function(plate) which(dat$protein_plate == plate))

for(p in prot_names){
  for(plate in plate_idx)
    datI[plate, p] <- add_noise_to_imputation(datI[plate, p])
}

### collect aggregated lagged times
# unique experiments
uniqueExp <- unique(datI[!(datI$pertLabel %in% c('no', 'no no')), 
                         c("Anchor_dose", "Library_dose", 
                           "pertLabel", 'protein_plate')])
str(uniqueExp)
expTimes <- c(6, 24, 48)

aggData <- lapply(expTimes, function(prev) {
  sapply(1:nrow(uniqueExp), function(i){
    experiment <- uniqueExp[i, ]
  
    anchor <- datI$Anchor_dose == experiment$Anchor_dose
    libr <- datI$Library_dose == experiment$Library_dose
    pert <- datI$pertLabel == experiment$pertLabel
    plate <- datI$protein_plate == experiment$protein_plate
    tim <- datI$pert_time == prev
    sel <- anchor & libr & pert & plate & tim
  
    if(sum(sel) == 0) print(experiment)
    t(colMeans(datI[sel, prot_names]))
  })
})

dim(aggData[[1]])

aggData[1:3, 1680]

# protein of interest
P <- "P08621.P08621.RU17_HUMAN.SNRNP70.U1.small.nuclear.ribonucleoprotein.70.kDa"
t <- 48

# Data for model
Y <- datI[datI$pert_time == t, P] - datI[datI$pert_time == t, paste0(P, "_0")]
D <- datI[datI$pert_time == t, pert_names]

## prepare design matrix with interactions
drug_design <- model.matrix(~ -1 + .^2, data = D)
dLabels <- c(colnames(drug_design), colnames(drug_design))

# remove tratements wihtout data
with_data <- apply(drug_design, 2, function(x) length(unique(x))) != 1
drug_design <- drug_design[, with_data]
dLabels_measured <- colnames(drug_design)
dlabels_model <- c(dLabels_measured, dLabels_measured)

# add intercept to each treatment
drug_intercept <- drug_design
drug_intercept[drug_design != 0] <- 1
colnames(drug_intercept) <- paste0(colnames(drug_design), "_intercept")

# combine intercept and drug design
drug_design <- cbind(drug_design, drug_intercept)
design <- drug_design

# for time 24 and 48 use lagged protein expression

datI[datI$pert_time == t, ]

design <- cbind(design, protein_design)


#labels
single_effects <- rep(c(colnames(D), rep(NA, choose(ncol(D), 2))), 2)
colnames(drug_design)

#hdi fit
fit <- lasso.proj(x = design, y = Y, ncores = 11, parallel = TRUE)

pMat <- matrix(NA, nrow = ncol(D), ncol = ncol(D))
rownames(pMat) <- colnames(pMat) <- dLabels[1:ncol(D)]

for(l in dLabels_measured){
  pval <- fit$groupTest(which(dlabels_model == l))
  drugs <- strsplit(l, ":")[[1]]
  if(length(drugs) == 1) drugs <- c(drugs, drugs)
  pMat[drugs[1], drugs[2]] <- pval
  pMat[drugs[2], drugs[1]] <- pval
}

pMat
pMat[is.na(pMat)] <- 2

pMat <- pMat[drugOrder, drugOrder]
#pMat <- -log(pMat)
pMat <- round(pMat, 3)


ht <- plot_ly(z = pMat, x = colnames(pMat), y = colnames(pMat), 
                type = "heatmap", colors = "Greys")
ht
