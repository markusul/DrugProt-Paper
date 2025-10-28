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
      #pval.drugs <- p.adjust(pval.drugs, method = 'holm')
      return(pval.drugs)
    }else{
      print(paste0(P, '_', t, " not found!"))
      return(rep(1, 122))
    }
  })
  pvec
})

# number of significant drug effects on different times
apply(array(p.adjust(unlist(allPvecs), method = "holm") < 0.05, 
            dim = c(122, 3, length(prot_names_short))), 
      2, sum)

# save p values as array with dimensions (treatment, time, protein)
treatment <- rownames(allPvecs[[1]])
allPvecs <- array(unlist(allPvecs), dim = c(122, 3, length(prot_names_short)))
save(allPvecs, treatment, file = "results/DrugEffects.RData")

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
    }else{
      net <- cbind(net, rep(0, length(prot_names_short)))
      targets <- c(targets, prot_names_short[P])
      print("missing experiment!")
    }
  }
  list(net, targets)
})

save(Net, file = "results/proteinNetwork.RData")

##### Protein Network ####
# collect all p values of protein on protein effects
nProt <- length(prot_names_short)

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
      links <- data.frame('source' = 1:nProt, 
                          'target' = Prot, 
                          'pvalue' = 2)
      print("missing experiment!")
    }
    links
  })
  
  do.call(rbind, Links)
})

save(Pval_all, file = "results/proteinNetworkPval.RData")
