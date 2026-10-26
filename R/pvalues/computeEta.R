# Bias factor of the de-sparsified lasso for every predictor (Zhang and Zhang, 2014).
# eta_j = max_{k != j} |Z_j' X_k| / ||Z_j||_2 for a design standardized to ||X_k||_2^2 = n.
# The design and Z are shared by all response proteins at a time point, so eta is a property of the predictor.
# For a drug group (concentration and intercept column) the partner column of the same group is excluded,
# because under the group null both coefficients are zero, and the group value is the larger of the two columns.
# Output results/eta.RData with
#   eta_drug     effect_type, treatment (index as in DrugEffects.RData), time_idx, eta
#   eta_protein  source (index in prot_names), time_idx (1 = 6h->24h, 2 = 24h->48h), eta
# Run from the repository root after getZ.sh.

load("data/laggedData.RData")
load("results/DrugEffects.RData")    # treatment, the label order used everywhere downstream
expTimes <- c(6, 24, 48)

bias_factor <- function(design, Z, labels_model) {
  stopifnot(nrow(Z) == nrow(design), ncol(Z) == ncol(design))
  n <- nrow(design)
  X <- scale(design) * sqrt(n / (n - 1))
  X[is.na(X)] <- 0
  C <- abs(crossprod(Z, X))
  diag(C) <- 0
  # partner column of the same treatment group
  for (g in split(seq_along(labels_model), labels_model)) 
    if (length(g) == 2){
      C[g[1], g[2]] <- C[g[2], g[1]] <- 0  
      eta <- apply(C, 1, max) / sqrt(colSums(Z^2))
    }
  rm(C, X); invisible(gc())
  eta
}

eta_drug <- list()
eta_protein <- list()

for (k in seq_along(expTimes)) {
  tp <- expTimes[k]

  # drug design exactly as in drugInteraction.R
  D <- datI[datI$pert_time == tp, pert_names]
  dd <- drugDesign(D)
  stopifnot(identical(dd$labels, treatment))
  nDrug <- length(dd$labels)
  nCol <- length(dd$labels_model)
  drug_design <- dd$design

  modes <- if (k == 1) "direct" else c("direct", "total")
  for (et in modes) {
    design <- drug_design
    if (et == "direct" && k > 1) {
      protein_design <- aggData[[k - 1]][datI[datI$pert_time == tp, "label"], ]
      protein_design <- protein_design - datI[datI$pert_time == tp, paste0(prot_names, "_0")]
      design <- cbind(design, as.matrix(protein_design))
      design <- design[rowSums(is.na(design)) == 0, ]
    }
    load(if (et == "total") paste0("Z/", tp, "_mode2.RData") else paste0("Z/", tp, ".RData"))

    eta <- bias_factor(design, Z, dd$labels_model)
    eta_group <- sapply(split(eta[1:nCol], factor(dd$labels_model, levels = dd$labels)), max)
    eta_drug[[length(eta_drug) + 1]] <- data.frame(
      effect_type = et, treatment = seq_len(nDrug), time_idx = k, eta = eta_group)
    if (et == "direct" && k > 1) {
      eta_protein[[length(eta_protein) + 1]] <- data.frame(
        source = seq_along(prot_names), time_idx = k - 1, eta = eta[-(1:nCol)])
    }
    message(sprintf("%d h, %s effects: median eta %.2f", tp, et, median(eta)))
  }
}

# at 6 h direct and total effects coincide
eta_drug <- do.call(rbind, eta_drug)
eta_drug <- rbind(eta_drug, transform(eta_drug[eta_drug$time_idx == 1, ], effect_type = "total"))
eta_protein <- do.call(rbind, eta_protein)

save(eta_drug, eta_protein, file = "results/eta.RData")
