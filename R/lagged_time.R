print("start lagged time data preparation")

load("data/prepData.RData")
load("data/protNames.RData")

# remove proteins with less than 11 unique values
abundand <- apply(data[, prot_names], 2, function(x) length(unique(x)) > 10)
prot_names <- prot_names[abundand]

# E[P0|C]
dat0 <- data[data$pert_time == 0, ]
dat0 <- aggregate(dat0[, prot_names], by = list(dat0$protein_plate), FUN = median)

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

# save order and short names
save(file = 'data/order.RData', prot_names_short, drugOrder, prot_names, 
     expTimes, nTreatment)

# add E[P0|C] to data
names(dat0)[-1] <- paste0(names(dat0)[-1], "_0")
row.names(dat0) <- dat0[, 1]
dat0 <- dat0[data$protein_plate, -1]

dat <- data
dat <- cbind(dat, dat0)

# remove M453(ATCC) (no baseline protein expression)
dat <- dat[dat$protein_plate != "M453(ATCC)", ]

# add nois to imputation for computational und statistical stability
add_noise_to_imputation <- function(x){
  if(length(unique(x)) < 3) return(x)
  
  border <- sort(unique(x))[2]
  imp <- which(x < border)
  x_ <- x[-imp]
  y <- x
  y[imp] <- truncnorm::rtruncnorm(length(imp), a = -Inf, b = border, mean = mean(x_), sd = sd(x_))
  y
}

# impute with truncated normal per protein per cell line
datI <- dat
plate_idx <- lapply(unique(dat$protein_plate), function(plate) which(dat$protein_plate == plate))

for(p in prot_names){
  for(plate in plate_idx){
    if(sum(is.na(add_noise_to_imputation(datI[plate, p]))) > 0)stop('nas produced')
    datI[plate, p] <- add_noise_to_imputation(datI[plate, p])
  }
}

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

    # mean over replicates
    colMeans(datI[sel, prot_names])
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

save(datI, prot_names, pert_names, aggData, file = 'data/laggedData.RData')
