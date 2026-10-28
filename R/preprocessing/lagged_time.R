print("start lagged time data preparation")
set.seed(22)

load("data/prepData.RData")
load("data/protNames.RData")

# remove proteins with less than 11 unique values
abundand <- apply(data[, prot_names], 2, function(x) length(unique(x)) > 10)
prot_names <- prot_names[abundand]

# sort drugs by information
drugOrder <- names(sort(colSums(data[, pert_names] != 0), decreasing = T))
drugOrder <- paste0('`', drugOrder, '`')

# collect protein names in short
prot_names_short <- sapply(prot_names, function(p) {
  parts     <- strsplit(p, '[.]')[[1]]
  human_idx <- grep("_HUMAN", parts)
  n         <- length(human_idx)
  # gene symbols are the n fields immediately following the LAST _HUMAN field
  last_human <- max(human_idx)
  genes      <- parts[(last_human + 1):(last_human + n)]
  paste0(genes, collapse = "/")
})

expTimes <- sort(unique(data$pert_time))[-1]
nTreatment <- length(unique(data[, "pertLabel"]))-2

dat <- data

# remove M453(ATCC) (no baseline protein expression)
dat <- dat[dat$protein_plate != "M453(ATCC)", ]

load("data/na_count.RData")   # share of missing values per protein, before the floor imputation

# impute values at the detection floor by the conditional mean of a normal fitted to the observed values, truncated at the smallest observed value
add_imputation <- function(x, floor){
  imp <- which(x <= floor + 1e-8)
  if(length(imp) == 0 || length(unique(x[-imp])) < 2) return(x)
  x_ <- x[-imp]
  x[imp] <- truncnorm::etruncnorm(a = -Inf, b = min(x_), mean = mean(x_), sd = sd(x_))
  x
}

datI <- dat
plate_idx <- lapply(unique(dat$protein_plate), function(plate) which(dat$protein_plate == plate))

for(p in prot_names){
  if(na_count[p] == 0) next          # protein without values at the floor
  floor_p <- min(data[, p])          # 80% of the minimum, global per protein
  for(plate in plate_idx){
    datI[plate, p] <- add_imputation(datI[plate, p], floor_p)
  }
  if(anyNA(datI[, p])) stop('nas produced')
}

# E[P0|C] from the imputed controls, so baseline and treated samples share the same imputation
dat0 <- datI[datI$pert_time == 0, ]
dat0 <- aggregate(dat0[, prot_names], by = list(dat0$protein_plate), FUN = median)
names(dat0)[-1] <- paste0(names(dat0)[-1], "_0")
row.names(dat0) <- dat0[, 1]
datI <- cbind(datI, dat0[datI$protein_plate, -1])

### collect aggregated lagged times
# unique experiments E[Pt|C,D]
uniqueExp <- unique(datI[, c("Anchor_dose", "Library_dose", "pertLabel", 'protein_plate')])

expTimes <- c(6, 24, 48)

aggData <- lapply(expTimes, function(prev) {
  t(sapply(1:nrow(uniqueExp), function(i){
    experiment <- uniqueExp[i, ]
    
    anchor <- datI$Anchor_dose == experiment$Anchor_dose
    libr <- datI$Library_dose == experiment$Library_dose
    pert <- datI$pertLabel == experiment$pertLabel
    plate <- datI$protein_plate == experiment$protein_plate
    tim <- datI$pert_time == prev
    sel <- anchor & libr & pert & plate & tim
    
    if(sum(sel) == 0) {
      # experiments that do not have this time point
      print(prev)
      print(experiment)
    }

    # median over replicates
    if(sum(sel) == 0) return(rep(NA_real_, length(prot_names)))
    apply(as.matrix(datI[sel, prot_names, drop = FALSE]), 2, median)
  }))
})

# add experiment label
datInfo <- datI[, c("Anchor_dose", "Library_dose", "pertLabel", 'protein_plate')]
label <- rep(NA, nrow(datInfo))

for(i in 1:nrow(uniqueExp)){
  experiment <- uniqueExp[i, ]
  anchor <- datI$Anchor_dose == experiment$Anchor_dose
  libr <- datI$Library_dose == experiment$Library_dose
  pert <- datI$pertLabel == experiment$pertLabel
  plate <- datI$protein_plate == experiment$protein_plate
  sel <- anchor & libr & pert & plate
  
  label[sel] <- i
}
datI$label <- label



# keep proteins whose differential expression varies at every time point, in the samples and in the lagged design
# proteins below the detection limit almost everywhere have a differential expression of exactly 0 and carry no information
# a column needs at least minInfo samples away from its most frequent value, so that it
# still varies in every training fold of the cross-validation in hdi
minInfo <- 10
varies <- function(m) apply(m, 2, function(x) {
  x <- x[!is.na(x)]
  u <- unique(x)
  sum(x != u[which.max(tabulate(match(x, u)))]) >= minInfo
})

diffP <- function(tp) as.matrix(datI[datI$pert_time == tp, prot_names]) -
  as.matrix(datI[datI$pert_time == tp, paste0(prot_names, "_0")])
lagP <- function(k) {
  sel <- datI$pert_time == expTimes[k + 1]
  aggData[[k]][datI$label[sel], ] - as.matrix(datI[sel, paste0(prot_names, "_0")])
}
# responses at 24 and 48 hours as used in the direct-effect models, on the samples with lagged measurements
lagRows <- function(k) !is.na(lagP(k)[, 1])
keep <- varies(diffP(6)) & varies(diffP(24)) & varies(diffP(48)) &
  varies(diffP(24)[lagRows(1), , drop = FALSE]) & varies(diffP(48)[lagRows(2), , drop = FALSE]) &
  varies(lagP(1)) & varies(lagP(2))
cat("proteins without variation removed:", sum(!keep), "of", length(keep), "\n")

prot_names <- prot_names[keep]
prot_names_short <- prot_names_short[keep]
aggData <- lapply(aggData, function(a) a[, keep])

# save order and short names
save(file = 'data/order.RData', prot_names_short, drugOrder, prot_names,
     expTimes, nTreatment)
save(datI, prot_names, pert_names, aggData, file = 'data/laggedData.RData')
