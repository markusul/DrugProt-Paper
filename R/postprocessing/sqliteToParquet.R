## sqliteToParquet.R
## --------------------------------------------------------------------------
## Convert each table of data/drugprot.sqlite into a compressed Parquet file
## under data/parquet/, then report sizes so you can compare against the
## 5.4 GB .sqlite. One Parquet file per table (DuckDB queries them together).
##
## Memory note: every table is read fully into R before writing. drug_pvalue
## (~3.9M rows with direct and total effects) and protein_pvalue (~58M rows on
## the real data) are the large
## ones; ~58M rows x 4 cols of doubles is ~1.8 GB in R briefly. If that is too
## much, see the chunked alternative for protein_pvalue at the bottom.
## --------------------------------------------------------------------------

library(DBI)
library(RSQLite)
library(arrow)

db_path  <- "data/drugprot.sqlite"
out_dir  <- "data/parquet"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(RSQLite::SQLite(), db_path)
tables <- dbListTables(con)
cat("Tables:", paste(tables, collapse = ", "), "\n\n")

## Sort the drug tables by effect type first, so a query for one effect type and
## a few proteins (or one treatment) touches only a few row groups.
order_by <- c(drug_pvalue = "effect_type, protein, treatment, time_idx",
              drug_coef   = "effect_type, treatment, protein, time_idx")

for (tbl in tables) {
  sql <- if (tbl %in% names(order_by)) {
    sprintf("SELECT * FROM %s ORDER BY %s", tbl, order_by[[tbl]])
  } else {
    sprintf("SELECT * FROM %s", tbl)
  }
  df  <- dbGetQuery(con, sql)
  out <- file.path(out_dir, paste0(tbl, ".parquet"))
  write_parquet(df, out, compression = "gzip")
  cat(sprintf("  %-16s %10d rows -> %s (%.1f MB)\n",
              tbl, nrow(df), out, file.info(out)$size / 1024^2))
  rm(df); gc()
}

dbDisconnect(con)

## ---- size summary ---------------------------------------------------------
sqlite_mb  <- file.info(db_path)$size / 1024^2
parquet_mb <- sum(file.info(list.files(out_dir, full.names = TRUE))$size) / 1024^2
cat(sprintf("\nSQLite total : %8.1f MB\n", sqlite_mb))
cat(sprintf("Parquet total: %8.1f MB  (%.1f%% of SQLite)\n",
            parquet_mb, 100 * parquet_mb / sqlite_mb))

## --------------------------------------------------------------------------
## If protein_pvalue is too big to hold in RAM at once, replace its loop
## iteration with a chunked writer (uses ParquetFileWriter):
##
##   library(arrow)
##   src <- dbSendQuery(con, "SELECT source,target,time_idx,pvalue
##                              FROM protein_pvalue ORDER BY source")
##   schema <- arrow::schema(source = int32(), target = int32(),
##                           time_idx = int32(), pvalue = float64())
##   writer <- ParquetFileWriter$create(
##     schema, sink = file.path(out_dir, "protein_pvalue.parquet"),
##     properties = ParquetWriterProperties$create(compression = "zstd"))
##   repeat {
##     ch <- dbFetch(src, n = 5e6)
##     if (nrow(ch) == 0) break
##     writer$WriteTable(arrow::as_arrow_table(ch), chunk_size = nrow(ch))
##   }
##   writer$Close(); dbClearResult(src)
##
## ORDER BY source makes a given protein's edges cluster into few row groups,
## so DuckDB can prune them on the source side of your `source IN (...)` filter.
## --------------------------------------------------------------------------
