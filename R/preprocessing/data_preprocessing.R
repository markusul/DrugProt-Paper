print("start data preprocessing")

# ---------------------------------------------------------------------------
# data_preprocessing.R
#
# Assembles the analysis input of Drug-Prot, one row per sample, holding
# log protein expressions, the drug concentration matrix, the experimental
# annotation and the viability readout.
#
# Inputs, all in data/:
#   01_sample_info.xlsx                sample annotation, db.prottalks.com
#   02_protein_intensity_matrix.xlsx   protein intensities, db.prottalks.com
#   protein_names.csv                  protein identifiers
#   analysis_samples.csv               sample selection and treatment labels
#   ic50_values.csv                    viability readout
#
# Outputs:
#   data/prepData.RData    the analysis input
#   data/drugLookup.RData  drug identifier to drug name
#   data/na_count.RData    fraction of missing values per protein
#
# The protein measurements come from the intensity matrix of Sun et al.
# (Nature 2026, doi:10.1038/s41586-026-11001-9). The three CSV files record
# the conventions of this analysis, which the matrix does not carry: which
# samples are used, the order in which each drug pair is written, the cell
# line naming, the protein column identifiers, and the viability values.
# See the Drug-Prot data deposit for their description.
#
# Two points about the sample selection. The released dataset contains more
# samples than are used here, among them 1117 untreated controls against the
# 840 of this analysis, which are those meeting complete biological
# triplicate and full factorial design criteria. Since the baseline is a
# per-cell-line median over the controls, the selection determines every
# differential expression. And of the 15002 selected samples, 158 are
# technical replicates that appear in the intensity matrix but not in the
# sample annotation, so 14844 enter the analysis.
# ---------------------------------------------------------------------------

library(dplyr)
library(readxl)
library(fastDummies)

# Impute missing values with 80% of minimum
# as we assume missing values are below detection limit
na_imputation <- function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE) * 0.8
  x
}

# ---- Read the inputs -------------------------------------------------------

sampleInfo <- read_excel(
  "data/01_sample_info.xlsx",
  col_types = c("text", "text", "text", "text", "text", "text", "text",
                "text", "numeric", "numeric", "text", "text", "text")
)

n_col <- ncol(read_excel("data/02_protein_intensity_matrix.xlsx", n_max = 0))
PTV1 <- read_excel(
  "data/02_protein_intensity_matrix.xlsx",
  col_types = c("text", "text", rep("numeric", n_col - 2))
)

protein_names    <- read.csv("data/protein_names.csv", stringsAsFactors = FALSE)
analysis_samples <- read.csv("data/analysis_samples.csv", stringsAsFactors = FALSE)
ic50_values      <- read.csv("data/ic50_values.csv", stringsAsFactors = FALSE)

protein_columns <- protein_names$protein[order(protein_names$row)]
stopifnot(nrow(PTV1) == length(protein_columns))

# Guard against a reordering of the intensity matrix. Compare the first gene
# symbol of each row against the corresponding entry of protein_names.csv.
# Protein groups separate symbols by semicolons in the matrix and by periods
# in the identifiers, and 68 rows carry no gene symbol, so those rows are
# excluded from the check rather than counted as failures.
first_gene <- sub("[;,/].*", "", PTV1$geneName)
checkable <- !is.na(first_gene) & nzchar(first_gene)
gene_ok <- mapply(function(g, nm) grepl(g, nm, fixed = TRUE),
                  first_gene[checkable], protein_columns[checkable])

cat("row order check:", sum(gene_ok), "of", sum(checkable), "checkable rows matched\n")
cat("rows without a gene symbol:", sum(!checkable), "\n")

if (mean(gene_ok) < 0.999) {
  stop("row order of the intensity matrix does not match protein_names.csv")
}

# ---- Transpose to one row per sample ---------------------------------------

intensity <- t(as.matrix(PTV1[, -(1:2)]))
colnames(intensity) <- protein_columns
sample_ids <- rownames(intensity)

# Restrict to the samples that entered the analysis and that carry annotation
keep <- sample_ids %in% analysis_samples$Sample_ID & sample_ids %in% sampleInfo$Sample_ID
intensity <- intensity[keep, , drop = FALSE]
sample_ids <- sample_ids[keep]

info <- sampleInfo[match(sample_ids, sampleInfo$Sample_ID), ]
sel  <- analysis_samples[match(sample_ids, analysis_samples$Sample_ID), ]
stopifnot(identical(info$Sample_ID, sample_ids))
stopifnot(identical(sel$Sample_ID, sample_ids))

cat("samples recovered:", nrow(intensity), "of", nrow(analysis_samples), "\n")
lost <- setdiff(analysis_samples$Sample_ID, sample_ids)
if (length(lost)) {
  cat("samples without annotation in the public release:", length(lost), "\n")
}

# ---- Protein data preprocessing --------------------------------------------

data_protein <- as.data.frame(intensity)

# Remove constant columns and keep only HUMAN proteins
data_protein <- data_protein[, apply(data_protein, 2, function(x) length(unique(x))) > 1]
data_protein <- data_protein[, grepl("HUMAN", names(data_protein))]

cat("protein columns after filtering:", ncol(data_protein), "\n")

# Analyze missing values
na_count <- colMeans(is.na(data_protein))
save(na_count, file = "data/na_count.RData")

# Impute and log transform
data_protein <- apply(data_protein, 2, na_imputation)
data_protein <- log(data_protein)

# ---- Drug annotation -------------------------------------------------------

# Taken from analysis_samples.csv rather than from the sample annotation,
# because the order within a drug pair and the cell line naming are
# conventions of this analysis.
pert_id      <- sel$pert_id
pertLabel    <- sel$pertLabel
Anchor_dose  <- sel$Anchor_dose
Library_dose <- sel$Library_dose
pert_time    <- sel$pert_time
cell_line    <- sel$cell_line

combination_idx <- which(pert_id == "")

# Cross-check the doses against the sample annotation, which records them
# per sample. Only combinations are compared, since the two sources use
# different conventions for single drugs.
pub_A <- suppressWarnings(as.numeric(info$Pert_Does1))
pub_L <- suppressWarnings(as.numeric(info$Pert_Does2))
mism <- combination_idx[
  abs(pub_A[combination_idx] - Anchor_dose[combination_idx]) > 1e-8 |
  abs(pub_L[combination_idx] - Library_dose[combination_idx]) > 1e-8
]
if (length(mism)) {
  warning(length(mism), " combination samples disagree on dose between the ",
          "analysis_samples.csv and the sample annotation")
}

# Create drug lookup
# What drug number corresponds to which drug name
drugnames <- unique(data.frame(pert_iname = info$Pert_Name1,
                               pert_id = info$Pert_ID1,
                               stringsAsFactors = FALSE))
drugnames <- drugnames[!is.na(drugnames$pert_id), ]
drug_lookup <- setNames(drugnames$pert_iname, drugnames$pert_id)

# Combination partners must also resolve to a name
partner <- unique(data.frame(pert_iname = info$Pert_Name2,
                             pert_id = info$Pert_ID2,
                             stringsAsFactors = FALSE))
partner <- partner[!is.na(partner$pert_id), ]
extra <- setdiff(partner$pert_id, names(drug_lookup))
if (length(extra)) {
  drug_lookup <- c(drug_lookup,
                   setNames(partner$pert_iname[match(extra, partner$pert_id)], extra))
}
save(drug_lookup, file = "data/drugLookup.RData")

# ---- Drug data preprocessing -----------------------------------------------

drugs <- pert_id
data_drugs <- dummy_cols(data.frame(drug = drugs), remove_first_dummy = FALSE)[, -1]
data_drugs <- data_drugs * 10

# Update drug combinations with actual doses
for (i in combination_idx) {
  drugAB <- paste0("drug_", strsplit(pertLabel[i], " ")[[1]])
  if (all(drugAB %in% colnames(data_drugs))) {
    data_drugs[i, drugAB[1]] <- Anchor_dose[i]
    data_drugs[i, drugAB[2]] <- Library_dose[i]
  }
}

# Remove reference columns
data_drugs <- data_drugs[, !colnames(data_drugs) %in% c("drug_", "drug_no")]

# ---- Additional metadata ---------------------------------------------------

# Both dose fields are 0 for single drugs, since the 10 micromole
# concentration is implicit in the drug indicator. data_preparation.R groups
# on these columns, so the convention matters.
out_A <- Anchor_dose
out_L <- Library_dose
is_single <- pert_id != "" & pert_id != "no"
out_A[is_single] <- 0
out_L[is_single] <- 0

data_additional <- data.frame(
  pert_time     = pert_time,
  protein_plate = cell_line,
  machine       = info$Machine,
  BioRep        = info$BioRep,
  Sample_ID     = info$Sample_ID,
  Anchor_dose   = out_A,
  Library_dose  = out_L,
  stringsAsFactors = FALSE
)

# ---- Response variable -----------------------------------------------------

# One IC50 per cell line, treatment and dose pairing, joined back to samples.
ic_key  <- paste(ic50_values$cell_line, ic50_values$pertLabel,
                 ic50_values$Anchor_dose, ic50_values$Library_dose, sep = "|")
dat_key <- paste(cell_line, pertLabel, Anchor_dose, Library_dose, sep = "|")

data_response <- ic50_values$IC50[match(dat_key, ic_key)]
cat("IC50 values matched:", sum(!is.na(data_response)), "of", length(data_response),
    "samples\n")

# ---- Combine all data ------------------------------------------------------

type <- rep("singleDrug", nrow(data_protein))
type[combination_idx] <- "drugCombination"
type[which(drugs == "no")] <- "noDrug"

data_response[drugs == "no"] <- Inf

data <- cbind(data_protein, data_drugs, data_additional, type, pertLabel,
              IC50 = data_response, NY = NA)
data[data$type == "singleDrug", "IC50"] <- log2(data[data$type == "singleDrug", "IC50"])

save(data, file = "data/prepData.RData")

cat("prepData.RData written:", nrow(data), "rows,", ncol(data), "columns\n")
print(table(data$type))
cat("cell lines:", length(unique(data$protein_plate)), "\n")
cat("distinct treatments:", length(unique(data$pertLabel)), "\n")
cat("distinct IC50 values:", length(unique(data$IC50[is.finite(data$IC50)])), "\n")
