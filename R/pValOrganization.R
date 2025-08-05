# load ordering of drugs (sort drugs with experiments together)
load('data/order.RData')
nDrugs <- length(drugOrder)

##### Drug Effects ####

# collect all p values for drug effects on all proteins at all time points
allPvecs <- lapply(1:length(prot_names_short), function(P){
  pvec <- sapply(c(6, 24, 48), function(t){
    if(file.exists(paste0('results/DrugEffects/', P , '_', t, '.RData'))){
      load(file = paste0('results/DrugEffects/', P , '_', t, '.RData'))
      
      # correction of p values using holm
      pval.drugs <- p.adjust(pval.drugs, method = 'holm')
      return(pval.drugs)
    }else{
      print(paste0(P, '_', t, " not found!"))
      return(rep(1, 122))
    }
  })
  if(!is.null(dim(pvec))) {
    psig <- apply(pvec, 2, function(p) sum(p.adjust(p, method = "holm", 
                                                    n = length(prot_names_short) * 3) < 0.05))
    pvec <- apply(pvec, 1, function(p) min(min(p) * 3, 1))
  }
  list(pvec, psig)
})

# number of significant drug effects on different times
rowSums(sapply(allPvecs, function(P) P[[2]]))

# collect all p values and organize as matrix
allPvecs <- lapply(allPvecs, function(P) P[[1]])
allPvecs <- allPvecs[unlist(lapply(allPvecs, length)) == 122]
allPvecs <- do.call(rbind, allPvecs)

save(allPvecs, file = "results/DrugEffects.RData")

##### Protein Network ####
# collect all p values of protein on protein effects
# correct for FD per target protein
Net <- lapply(c(24, 48), function(t){
  net <- NULL
  targets <- c() # response proteins
  
  for (P in 1:length(prot_names)) {
    path <- paste0('results/ProteinEffects/', P , '_', t, '.RData')
    if(file.exists(path)){
      load(file = path)
      # select protein effects (leave drug effects)
      # apply Benjamini & Hochberg (1995) for fdr
      pval.corr <- p.adjust(pval[(length(pval)-length(prot_names_short) + 1):length(pval)], 
                            method = "BH")
      net <- cbind(net, pval.corr) # add p values to network
      targets <- c(targets, prot_names_short[P]) # add response
    }
  }
  list(net, targets)
})

save(Net, file = "results/proteinNetwork.RData")