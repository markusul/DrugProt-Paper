## buildDatabase.R
## --------------------------------------------------------------------------
## Collect all Drug-Prot data (drug effects, protein-network p-values,
## drug/protein coefficients, and metadata) into ONE SQLite database that the
## Shiny app can query by protein selection without loading whole arrays.
##
##
## Design notes
##  * Long/tidy tables, not wide arrays, so we can index on protein and pull
##    only the rows a query touches.
##  * The drug-effect array becomes a table indexed on `protein`.
##  * The protein-network matrix becomes a table indexed on (source,target);
##    sentinel p-values (== 2, "no experiment") are dropped, not stored.
##  * Coefficients (drug `dEff`, protein `bhat`) move from thousands of tiny
##    .RData files into two indexed tables.
##  * Metadata vectors go into small key/value-ish tables.
##  * Everything is written in chunks inside a single transaction per table so
##    the build itself never holds a whole array in RAM longer than needed.
##  * Each table is filled via dp_insert(): one dbSendStatement -> dbBind ->
##    dbClearResult cycle per batch. A single result handle must NOT be re-bound
##    in a loop (that raises "Invalid result set" / leaves pending rows), so we
##    bind ONCE per batch with whole column vectors -- which is also faster.
## --------------------------------------------------------------------------

library(DBI)
library(RSQLite)

## ---- configuration --------------------------------------------------------
data_dir <- "data"                                   # where the .RData live
results_dir <- "results"
db_path  <- file.path(data_dir, "drugprot.sqlite")   # output
chunk_proteins <- 3000L   # proteins per insert chunk (tune for RAM)

if (file.exists(db_path)) {
  message("Removing existing ", db_path)
  file.remove(db_path)
}

con <- dbConnect(RSQLite::SQLite(), db_path)
## Faster bulk load; these are session pragmas, safe for a build.
dbExecute(con, "PRAGMA journal_mode = OFF;")
dbExecute(con, "PRAGMA synchronous = OFF;")
dbExecute(con, "PRAGMA cache_size = -200000;")  # ~200 MB page cache

## A tiny helper for loading a single object out of an .RData by name.
load_one <- function(path, name) {
  e <- new.env()
  load(path, envir = e)
  get(name, envir = e)
}

## Insert a batch of rows in one prepared-statement cycle.
## `cols` is a named list of equal-length column vectors, in table-column order.
## Wrap calls in dbBegin/dbCommit yourself when batching many dp_insert()s into
## one transaction (cheaper than a transaction per call).
dp_insert <- function(con, table, cols) {
  n <- length(cols[[1]])
  if (n == 0) return(invisible())
  ph  <- paste(rep("?", length(cols)), collapse = ",")
  sql <- sprintf("INSERT INTO %s (%s) VALUES (%s)",
                 table, paste(names(cols), collapse = ","), ph)
  rs <- dbSendStatement(con, sql)
  on.exit(dbClearResult(rs), add = TRUE)   # always clears, even on error
  dbBind(rs, unname(cols))
  invisible(n)
}

## ===========================================================================
## 1. METADATA
## ===========================================================================
message("== Metadata ==")

load(file.path(data_dir, "order.RData"))        # drugOrder, expTimes, prot_names_short, ...
load(file.path(data_dir, "drugLookup.RData"))   # drug_lookup
load(file.path(results_dir, "Coef/drugs/treatNames.RData"))  # treatNames (named)

nProtein   <- length(prot_names_short)
nTreatment <- 122L
expTimes   <- as.integer(expTimes)
nTimes     <- length(expTimes)

## proteins: 1-based index <-> HGNC short name (the app's prot_names_short order)
dbWriteTable(con, "proteins",
  data.frame(protein = seq_len(nProtein),
             name    = unname(prot_names_short),
             stringsAsFactors = FALSE),
  overwrite = TRUE)

## treatments: 1-based index <-> raw model label <-> display name.
## `treatment` (raw rownames) comes from DrugEffects.RData; load it here.
treatment_raw <- load_one(file.path(results_dir, "DrugEffects.RData"), "treatment")
dbWriteTable(con, "treatments",
  data.frame(treatment = seq_len(nTreatment),
             raw_label = treatment_raw,
             # treatNames is named by display name -> raw label; invert for lookup
             display   = names(treatNames)[match(treatment_raw, unname(treatNames))],
             stringsAsFactors = FALSE),
  overwrite = TRUE)

## expTimes, drugOrder, drug_lookup, and assorted scalars as small tables.
dbWriteTable(con, "exp_times",
  data.frame(idx = seq_along(expTimes), hours = expTimes),
  overwrite = TRUE)

dbWriteTable(con, "drug_order",
  data.frame(ord = seq_along(drugOrder), raw_label = unname(drugOrder)),
  overwrite = TRUE)

dbWriteTable(con, "drug_lookup",
  data.frame(id = names(drug_lookup), name = unname(drug_lookup),
             stringsAsFactors = FALSE),
  overwrite = TRUE)

dbWriteTable(con, "meta",
  data.frame(key   = c("nProtein", "nTreatment", "nTimes"),
             value = c(nProtein,   nTreatment,   nTimes)),
  overwrite = TRUE)

## ===========================================================================
## 2. DRUG EFFECTS  (array -> indexed long table)
##    allPvecs[treatment, time, protein]  ->  rows (protein, treatment, time, pvalue)
## ===========================================================================
message("== Drug-effect p-values ==")

allPvecs <- load_one(file.path(results_dir, "DrugEffects.RData"), "allPvecs")
stopifnot(identical(dim(allPvecs), c(nTreatment, nTimes, nProtein)))

dbExecute(con, "
  CREATE TABLE drug_pvalue (
    protein   INTEGER NOT NULL,
    treatment INTEGER NOT NULL,
    time_idx  INTEGER NOT NULL,
    pvalue    REAL    NOT NULL
  );")

## Insert in protein-chunks. For each chunk we melt the [nTreatment x nTimes x
## |chunk|] sub-array to long form and bind ONCE.
prot_chunks <- split(seq_len(nProtein),
                     ceiling(seq_len(nProtein) / chunk_proteins))

## within one protein's slice, rows run treatment-fastest then time;
## as.numeric() on a [nTreatment x nTimes x |chunk|] array flattens
## column-major (treatment, then time, then protein) -- matching these grids.
tr_grid <- rep(seq_len(nTreatment), times = nTimes)
t_grid  <- rep(seq_len(nTimes),     each  = nTreatment)

for (ch in prot_chunks) {
  n_per_prot    <- nTreatment * nTimes
  protein_col   <- rep(ch, each = n_per_prot)
  treatment_col <- rep(tr_grid, times = length(ch))
  time_col      <- rep(t_grid,  times = length(ch))
  pvalue_col    <- as.numeric(allPvecs[, , ch])   # column-major: matches grids

  dbBegin(con)
  dp_insert(con, "drug_pvalue", list(
    protein   = protein_col,
    treatment = treatment_col,
    time_idx  = time_col,
    pvalue    = pvalue_col))
  dbCommit(con)
  message(sprintf("  drug_pvalue: proteins %d-%d", min(ch), max(ch)))
}
rm(allPvecs); gc()

## ===========================================================================
## 3. PROTEIN NETWORK P-VALUES
##    pvalue[row, time]  where row = (target-1)*nProtein + source  (1-based)
##    -> rows (source, target, time_idx, pvalue), sentinel 2 dropped.
## ===========================================================================
message("== Protein-network p-values ==")

## Load once, write in target-chunks, then free it.
pvalue_mat <- load_one(file.path(results_dir, "proteinNetworkPval_pvalue.RData"),
                       "pvalue")
nNetTimes <- ncol(pvalue_mat)   # = nTimes - 1 (the 24h, 48h transitions)
stopifnot(nrow(pvalue_mat) == nProtein * nProtein)

dbExecute(con, "
  CREATE TABLE protein_pvalue (
    source   INTEGER NOT NULL,
    target   INTEGER NOT NULL,
    time_idx INTEGER NOT NULL,
    pvalue   REAL    NOT NULL
  );")

## Recover (source, target) from the flattened index exactly as app.R does:
##   source = (idx-1) %% nProtein + 1 ;  target = ceiling(idx / nProtein)
## i.e. flat row = (target-1)*nProtein + source. We process targets in chunks,
## accumulate the kept (non-sentinel) rows for the whole chunk, bind ONCE.
target_chunks <- split(seq_len(nProtein),
                       ceiling(seq_len(nProtein) / chunk_proteins))

for (ch in target_chunks) {
  src_acc <- integer(0); tgt_acc <- integer(0)
  tt_acc  <- integer(0); pv_acc  <- numeric(0)

  for (tgt in ch) {
    base <- (tgt - 1L) * nProtein
    rows <- base + seq_len(nProtein)        # the nProtein flattened rows for this target
    src  <- seq_len(nProtein)
    for (tt in seq_len(nNetTimes)) {
      pv   <- pvalue_mat[rows, tt]
      keep <- which(pv < 2)                 # drop "no experiment" sentinels
      if (length(keep)) {
        src_acc <- c(src_acc, src[keep])
        tgt_acc <- c(tgt_acc, rep(tgt, length(keep)))
        tt_acc  <- c(tt_acc,  rep(tt,  length(keep)))
        pv_acc  <- c(pv_acc,  pv[keep])
      }
    }
  }

  dbBegin(con)
  dp_insert(con, "protein_pvalue", list(
    source   = src_acc,
    target   = tgt_acc,
    time_idx = tt_acc,
    pvalue   = pv_acc))
  dbCommit(con)
  message(sprintf("  protein_pvalue: targets %d-%d (%d rows)",
                  min(ch), max(ch), length(pv_acc)))
}
rm(pvalue_mat); gc()

## ===========================================================================
## 4. DRUG COEFFICIENTS  (results/Coef/drugs/{drug}.RData : dEff[protein, time])
## ===========================================================================
message("== Drug coefficients ==")

dbExecute(con, "
  CREATE TABLE drug_coef (
    treatment INTEGER NOT NULL,
    protein   INTEGER NOT NULL,
    time_idx  INTEGER NOT NULL,
    coef      REAL    NOT NULL
  );")

pg <- rep(seq_len(nProtein), times = nTimes)
tg <- rep(seq_len(nTimes),   each  = nProtein)

dbBegin(con)
for (drug in seq_len(nTreatment)) {
  f <- file.path(results_dir, "Coef/drugs", paste0(drug, ".RData"))
  if (!file.exists(f)) { message("  missing drug coef ", f); next }
  dEff <- load_one(f, "dEff")              # nProtein x nTimes
  dEff <- matrix(dEff, nrow = nProtein)
  vals <- as.numeric(dEff)
  keep <- which(!is.na(vals))              # store only present coefficients
  dp_insert(con, "drug_coef", list(
    treatment = rep(drug, length(keep)),
    protein   = pg[keep],
    time_idx  = tg[keep],
    coef      = vals[keep]))
}
dbCommit(con)

## ===========================================================================
## 5. PROTEIN COEFFICIENTS (results/Coef/proteins/{protein}_{hours}.RData : bhat)
##    Used for edge direction (sign of bhat) in the temporal graph.
## ===========================================================================
message("== Protein coefficients ==")

dbExecute(con, "
  CREATE TABLE protein_coef (
    target   INTEGER NOT NULL,   -- the protein whose model this bhat belongs to
    source   INTEGER NOT NULL,   -- index into bhat (predictor protein)
    time_idx INTEGER NOT NULL,   -- time of the target measurement
    coef     REAL    NOT NULL
  );")

## time index here matches expTimes[-1] (24h, 48h) -> time_idx 2..nTimes
for (ti in 2:nTimes) {
  hrs <- expTimes[ti]
  dbBegin(con)
  for (protein in seq_len(nProtein)) {
    f <- file.path(results_dir, "Coef/proteins",
                   paste0(protein, "_", hrs, ".RData"))
    if (!file.exists(f)) next
    bhat <- load_one(f, "bhat")
    bhat <- as.numeric(bhat) > 0
    keep <- which(bhat != 0 & !is.na(bhat))   # only nonzero coefficients matter
    dp_insert(con, "protein_coef", list(
      target   = rep(protein, length(keep)),
      source   = keep,
      time_idx = rep(ti, length(keep)),
      coef     = bhat[keep]))
  }
  dbCommit(con)
  message(sprintf("  protein_coef: time %dh done", hrs))
}

## ===========================================================================
## 6. INDEXES  (built after bulk load = much faster)
## ===========================================================================
message("== Indexes ==")

dbExecute(con, "CREATE INDEX idx_drug_pval_protein  ON drug_pvalue(protein);")
dbExecute(con, "CREATE INDEX idx_net_source         ON protein_pvalue(source, time_idx);")
dbExecute(con, "CREATE INDEX idx_net_target         ON protein_pvalue(target, time_idx);")
dbExecute(con, "CREATE INDEX idx_drug_coef_treat    ON drug_coef(treatment);")
dbExecute(con, "CREATE INDEX idx_prot_coef_target   ON protein_coef(target, time_idx);")
dbExecute(con, "CREATE INDEX idx_proteins_name      ON proteins(name);")

message("== Finalizing ==")
dbExecute(con, "PRAGMA journal_mode = WAL;")   # good for concurrent Shiny reads
dbExecute(con, "VACUUM;")
dbExecute(con, "ANALYZE;")
dbDisconnect(con)

message("Done -> ", db_path)
