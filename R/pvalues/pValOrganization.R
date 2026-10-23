# load ordering of drugs (sort drugs with experiments together)
load('data/order.RData')
nDrugs <- length(drugOrder)
nProt <- length(prot_names_short)
nTreatment <- 122
nTimes <- 3

##### Drug Effects ####
# collect all p values for drug effects on all proteins at all time points
allPvecs <- lapply(1:nProt, function(P){
  pvec <- sapply(c(6, 24, 48), function(t){
    if(file.exists(paste0('results/DrugEffects/', P , '_', t, '.RData'))){
      load(file = paste0('results/DrugEffects/', P , '_', t, '.RData'))
      
      # correction of p values using holm
      #pval.drugs <- p.adjust(pval.drugs, method = 'holm')
      return(pval.drugs)
    }else{
      stop(paste0("results/DrugEffects/", P, "_", t, ".RData not found"))
    }
  })
  pvec
})

# number of significant drug effects on different times
apply(array(p.adjust(unlist(allPvecs), method = "holm") < 0.05, 
            dim = c(nTreatment, nTimes, nProt)), 
      2, sum)

# save p values as array with dimensions (treatment, time, protein)
treatment <- rownames(allPvecs[[1]])
allPvecs <- array(unlist(allPvecs), dim = c(nTreatment, nTimes, nProt))
save(allPvecs, treatment, file = "results/DrugEffects.RData")

##### Protein Network ####
# collect all p values of protein on protein effects
Pval_all <- lapply(c(24, 48), function(t){
  Links <- lapply(1:nProt, function(Prot) {
    path <- paste0('results/ProteinEffects/', Prot , '_', t, '.RData')
    if(file.exists(path)){
      load(file = path)
      # select protein effects (leave drug effects)
      pval <- unname(pval[(length(pval)-nProt + 1):length(pval)])
      links <- data.frame('source' = 1:nProt, 
                          'target' = Prot, 
                          'pvalue' = pval)
    }else{
      stop(paste0(path, " not found"))
    }
    links
  })
  
  do.call(rbind, Links)
})

save(Pval_all, file = "results/proteinNetworkPval.RData")

# save p values minimalistically 
pvalue <- lapply(Pval_all, function(P) P[, "pvalue"])
pvalue <- do.call(cbind, pvalue)

save(pvalue, file = "results/proteinNetworkPval_pvalue.RData")

##### Drug Effects mode 2####
# total causal effects of drugs onto the proteins
# collect all p values for drug effects on all proteins at all time points
# total effects, 6 h has no lagged proteins so direct and total effect coincide
allPvecsTotal <- lapply(1:nProt, function(P){
  sapply(c(6, 24, 48), function(t){
    suffix <- if(t == 6) "" else "_mode2"
    path <- paste0('results/DrugEffects/', P, '_', t, suffix, '.RData')
    if(!file.exists(path)) stop(paste0(path, " not found"))
    load(file = path)
    pval.drugs
  })
})
stopifnot(identical(rownames(allPvecsTotal[[1]]), treatment))
allPvecsTotal <- array(unlist(allPvecsTotal), dim = c(nTreatment, nTimes, nProt))
save(allPvecsTotal, treatment, file = "results/DrugEffects_total.RData")

# number of significant drug effects on different times
apply(array(p.adjust(unlist(allPvecs), method = "holm") < 0.05, 
            dim = c(nTreatment, nTimes, nProt)), 
      2, sum)
