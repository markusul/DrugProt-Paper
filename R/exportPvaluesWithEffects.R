## exportPvaluesWithEffects.R
## ===========================================================================
## Like exportPvalues.R, but each row also carries the effect MAGNITUDE that
## the Shiny app shows:
##   * drug effects  -> dEff  (average treatment effect)  from
##                      results/Coef/drugs/{drug}.RData      [protein x time]
##   * protein edges -> bhat  (de-sparsified coefficient)   from
##                      results/Coef/proteins/{protein}_{hours}.RData
##
## This is "Option 1": coefficients are read straight from the results/ tree,
## NOT from the Parquet store (which only kept the SIGN of bhat). It must
## therefore run where results/ still lives (the cluster), not on shinyapps.
##
## Effects are attached PER TIME POINT (one effect per p-value row), which is
## more complete than the app's min-over-time collapse. Missing/absent
## coefficients are filled with 0, matching app.R (out[is.na(out)] <- 0).
##
## Outputs (out_dir, default data/downloads/):
##   drug_pvalues_effects.csv.gz     (protein, drug, time_h, pvalue, effect)
##   protein_network_effects.csv.gz  (source, target, time_transition, pvalue, effect)
##
## Usage:
##   Rscript R/exportPvaluesWithEffects.R [parquet_dir] [out_dir] [results_dir]
## ===========================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

args        <- commandArgs(trailingOnly = TRUE)
parquet_dir <- if (length(args) >= 1) args[[1]] else "data/parquet"
out_dir     <- if (length(args) >= 2) args[[2]] else "data/downloads"
results_dir <- if (length(args) >= 3) args[[3]] else "results"

stopifnot(dir.exists(parquet_dir), dir.exists(results_dir))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pq <- function(tbl) {
  p <- file.path(parquet_dir, paste0(tbl, ".parquet"))
  if (!file.exists(p)) stop("missing Parquet table: ", p)
  sprintf("read_parquet('%s')", normalizePath(p, winslash = "/"))
}
load_one <- function(path, name) {           # pull one named object from an .RData
  e <- new.env(); load(path, envir = e); get(name, envir = e)
}

con <- dbConnect(duckdb::duckdb())
## (no on.exit: under source() it can close the connection mid-script)

## Pull the dimension/label tables into R once -- all small.
proteins   <- dbGetQuery(con, sprintf("SELECT protein, name FROM %s ORDER BY protein", pq("proteins")))
treatments <- dbGetQuery(con, sprintf("SELECT treatment, display FROM %s ORDER BY treatment", pq("treatments")))
exp_times  <- dbGetQuery(con, sprintf("SELECT idx, hours FROM %s ORDER BY idx", pq("exp_times")))
nProt      <- nrow(proteins)
prot_name  <- setNames(proteins$name, proteins$protein)
drug_name  <- setNames(treatments$display, treatments$treatment)
hours_of   <- setNames(exp_times$hours, exp_times$idx)          # time_idx -> hours

## ===========================================================================
## 1. DRUG EFFECTS  + dEff magnitude
##    Read the p-values from Parquet, attach dEff[protein, time] per row.
##    dEff files are results/Coef/drugs/{treatment}.RData, shape [nProt x nTimes].
## ===========================================================================
message("Drug effects: reading p-values ...")
drug_pv <- dbGetQuery(con, sprintf(
  "SELECT protein, treatment, time_idx, pvalue FROM %s", pq("drug_pvalue")))

## Build a (treatment, protein, time_idx) -> dEff lookup by loading each drug file.
message("Drug effects: attaching dEff magnitudes ...")
drug_pv$effect <- 0                                  # default fill (matches app)
nTimes <- nrow(exp_times)
for (tr in sort(unique(drug_pv$treatment))) {
  f <- file.path(results_dir, "Coef/drugs", paste0(tr, ".RData"))
  if (!file.exists(f)) { message("  missing drug coef: ", f); next }
  dEff <- load_one(f, "dEff")                        # [nProt x nTimes]
  dEff <- matrix(dEff, nrow = nProt)
  stopifnot(nrow(dEff) == nProt, ncol(dEff) == nTimes)
  rows <- which(drug_pv$treatment == tr)
  ## index dEff by this drug's (protein, time_idx) rows; NA -> 0
  vals <- dEff[cbind(drug_pv$protein[rows], drug_pv$time_idx[rows])]
  vals[is.na(vals)] <- 0
  drug_pv$effect[rows] <- vals
}

## Relabel to human-readable and write.
drug_out <- file.path(out_dir, "drug_pvalues_effects.csv.gz")
drug_df  <- data.frame(
  protein = prot_name[as.character(drug_pv$protein)],
  drug    = drug_name[as.character(drug_pv$treatment)],
  time_h  = hours_of[as.character(drug_pv$time_idx)],
  pvalue  = drug_pv$pvalue,
  effect  = drug_pv$effect,
  stringsAsFactors = FALSE
)
drug_df <- drug_df[order(drug_df$protein, drug_df$drug, drug_df$time_h), ]
gz <- gzfile(drug_out, "w")
write.csv(drug_df, gz, row.names = FALSE)
close(gz)
message(sprintf("  wrote %s (%.1f MB, %d rows)",
                basename(drug_out), file.info(drug_out)$size / 1024^2, nrow(drug_df)))
rm(drug_pv, drug_df); gc()

## ===========================================================================
## 2. PROTEIN NETWORK  + bhat magnitude   (the ~58M-row table)
##    Stream from Parquet in target-chunks. For each target protein, load its
##    bhat files (one per transition time) and attach by source index.
##    NOTE: bhat in results/Coef/proteins/{P}_{hours}.RData is ALREADY the
##    protein block, length == nProt, indexed by source protein. No slicing.
## ===========================================================================
message("Protein network: streaming with bhat magnitudes ...")

net_out <- file.path(out_dir, "protein_network_effects.csv.gz")
gz <- gzfile(net_out, "w")
write.csv(data.frame(source = character(), target = character(),
                     time_transition = character(), pvalue = double(),
                     effect = double()),
          gz, row.names = FALSE)            # header only

## network time_idx -> the later hours value (transition t -> t+1)
## expTimes are 6,24,48; transitions are 6->24 (idx 1) and 24->48 (idx 2).
trans_label <- c("1" = "6h->24h", "2" = "24h->48h")
trans_hours <- c("1" = 24, "2" = 48)        # hours used in the bhat filename

chunk <- 300L
target_chunks <- split(seq_len(nProt), ceiling(seq_len(nProt) / chunk))

for (ch in target_chunks) {
  ## pull all p-value rows whose target is in this chunk, in one DuckDB query
  net <- dbGetQuery(con, sprintf(
    "SELECT source, target, time_idx, pvalue FROM %s
     WHERE target BETWEEN %d AND %d", pq("protein_pvalue"), min(ch), max(ch)))
  if (nrow(net) == 0) next

  net$effect <- 0
  ## attach bhat: for each (target, time_idx) load the protein's bhat vector once
  for (tgt in unique(net$target)) {
    for (ti in unique(net$time_idx[net$target == tgt])) {
      hrs <- trans_hours[as.character(ti)]
      f <- file.path(results_dir, "Coef/proteins", paste0(tgt, "_", hrs, ".RData"))
      if (!file.exists(f)) next
      bhat <- as.numeric(load_one(f, "bhat"))
      stopifnot(length(bhat) == nProt)         # guardrail: must be the protein block
      idx <- which(net$target == tgt & net$time_idx == ti)
      vals <- bhat[net$source[idx]]
      vals[is.na(vals)] <- 0
      net$effect[idx] <- vals
    }
  }

  out <- data.frame(
    source          = prot_name[as.character(net$source)],
    target          = prot_name[as.character(net$target)],
    time_transition = trans_label[as.character(net$time_idx)],
    pvalue          = net$pvalue,
    effect          = net$effect,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$source, out$target, out$time_transition), ]
  write.table(out, gz, sep = ",", row.names = FALSE,
              col.names = FALSE, qmethod = "double")
  message(sprintf("  targets %d-%d: %d rows", min(ch), max(ch), nrow(out)))
  rm(net, out); gc()
}
close(gz)
message(sprintf("  wrote %s (%.1f MB)",
                basename(net_out), file.info(net_out)$size / 1024^2))

dbDisconnect(con, shutdown = TRUE)
message("Done. Files written to ", normalizePath(out_dir, winslash = "/"))
