# drug design of models (1) and (2)
# columns are the concentration terms of treatments measured at more than one concentration, followed by the intercepts of all treatments
# for treatments measured at a single concentration the concentration term is a multiple of the intercept and is left out
drugDesign <- function(D){
  conc <- model.matrix(~ -1 + .^2, data = D)
  conc <- conc[, apply(conc, 2, function(x) length(unique(x))) != 1]
  labels <- colnames(conc)
  intercept <- conc
  intercept[conc != 0] <- 1
  colnames(intercept) <- paste0(labels, "_intercept")
  multi <- apply(conc, 2, function(x) length(unique(x[x != 0]))) > 1
  list(design = cbind(conc[, multi, drop = FALSE], intercept),
       labels = labels,
       labels_model = c(labels[multi], labels))
}