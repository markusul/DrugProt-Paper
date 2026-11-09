## appendixC.R
## ===========================================================================
## Numbers of S1 Appendix, Section C, the case study on drugs that change the
## abundance of their target. Reads the Parquet store that backs the Shiny
## application (data/parquet, written by sqliteToParquet.R) and writes
## results/C_Results.txt and results/C_table.csv.
##
## The corrections are those of the application. Drug effects: Holm over the
## queried proteins, all 122 treatments and the three time points, direct
## effects at 24 and 48 hours. Network: Benjamini-Hochberg over all tests that
## involve a queried protein as parent or child at both transitions.
##
##   Rscript R/postprocessing/appendixC.R
## ===========================================================================

library(DBI)
library(duckdb)

parquet_dir <- "data/parquet"
pq <- function(tbl) sprintf("read_parquet('%s')", file.path(parquet_dir, paste0(tbl, ".parquet")))
con <- dbConnect(duckdb::duckdb())

alpha <- 0.05
alpha_explore <- 0.15

## ---- metadata ---------------------------------------------------------------
proteins   <- dbGetQuery(con, sprintf("SELECT protein, name FROM %s ORDER BY protein", pq("proteins")))
treatments <- dbGetQuery(con, sprintf("SELECT treatment, display FROM %s ORDER BY treatment", pq("treatments")))
exp_times  <- dbGetQuery(con, sprintf("SELECT idx, hours FROM %s ORDER BY idx", pq("exp_times")))
nProt  <- nrow(proteins)
nTreat <- nrow(treatments)
nTimes <- nrow(exp_times)
hours  <- setNames(exp_times$hours, exp_times$idx)
is_pair <- grepl(":", treatments$display)
stopifnot(nTreat == 122, nTimes == 3)

C <- c("Results of S1 Appendix, Section C, drugs that change the abundance of their target",
       "-----------------------------------",
       paste0("Proteins in the store: ", nProt),
       paste0("Treatments: ", nTreat, ", single drugs: ", sum(!is_pair), ", pairs: ", sum(is_pair)))

## ---- queried proteins and predictions --------------------------------------
targets <- c("DHFR", "TYMS", "DNMT1", "TOP1", "TOP2B", "POLR2A", "FKBP5", "PARP1", "EGFR")
sel <- proteins$protein[match(targets, proteins$name)]
stopifnot(!anyNA(sel))
M <- length(sel)
inSel <- paste(sel, collapse = ",")

pairs <- rbind(
  data.frame(target = "DHFR",   drug = c("Methotrexate", "Pemetrexed Disodium Hydrate"), expected = "up"),
  data.frame(target = "TYMS",   drug = c("5-Fluorouracil", "Pemetrexed Disodium Hydrate"), expected = "up"),
  data.frame(target = "TYMS",   drug = "Capecitabine", expected = "none"),
  data.frame(target = "DNMT1",  drug = c("Decitabine", "Azacytidine"), expected = "down"),
  data.frame(target = "TOP1",   drug = "Irinotecan hydrochloride", expected = "down"),
  data.frame(target = "TOP2B",  drug = c("Doxorubicin Hydrochloride", "Epirubicin Hydrochloride"), expected = "down"),
  data.frame(target = "POLR2A", drug = "Lurbinectedin", expected = "down"),
  data.frame(target = "FKBP5",  drug = c("hydrocortisone", "Medroxyprogesterone Acetate"), expected = "up"),
  data.frame(target = "PARP1",  drug = c("Olaparib", "Niraparib Tosylate", "Talazoparib"), expected = "none"),
  data.frame(target = "EGFR",   drug = c("Afatinib Dimaleate", "Erlotinib Hydrochloride", "Lapatinib Ditosylate Hydrate",
                                         "Neratinib Maleate", "Tucatinib"), expected = "none"))
stopifnot(all(pairs$drug %in% treatments$display))

## ---- drug effects, as the application computes them ------------------------
dp <- dbGetQuery(con, sprintf(
  "SELECT effect_type, protein, treatment, time_idx, pvalue FROM %s WHERE protein IN (%s)",
  pq("drug_pvalue"), inSel))
dc <- dbGetQuery(con, sprintf(
  "SELECT effect_type, protein, treatment, time_idx, coef FROM %s WHERE protein IN (%s)",
  pq("drug_coef"), inSel))
dp <- merge(dp, dc, by = c("effect_type", "protein", "treatment", "time_idx"), all.x = TRUE)
dp$coef[is.na(dp$coef)] <- 0            # absent coefficients are 0, as in the application
dp$target <- proteins$name[dp$protein]
dp$drug   <- treatments$display[dp$treatment]
dp$hours  <- hours[as.character(dp$time_idx)]
dp$pair   <- is_pair[dp$treatment]

d <- dp[dp$effect_type == "direct", ]
stopifnot(nrow(d) == M * nTreat * nTimes)
d$p_holm <- p.adjust(d$pvalue, "holm")   # Holm over the set, all treatments and all time points
d <- d[order(d$target, d$drug, d$hours), ]

fmt <- function(p) ifelse(p >= 0.01, sprintf("%.2f", p), sprintf("%.1e", p))
arrow <- function(coef, p) ifelse(p < alpha, ifelse(coef > 0, "up ", ifelse(coef < 0, "down ", "sign? ")), "")

## the table of Section C
tab <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  r <- d[d$target == pairs$target[i] & d$drug == pairs$drug[i], ]
  r <- r[order(r$hours), ]
  data.frame(target = pairs$target[i], drug = pairs$drug[i], expected = pairs$expected[i],
             p6 = r$p_holm[1], p24 = r$p_holm[2], p48 = r$p_holm[3],
             sign6 = sign(r$coef[1]), sign24 = sign(r$coef[2]), sign48 = sign(r$coef[3]),
             raw6 = r$pvalue[1], raw24 = r$pvalue[2], raw48 = r$pvalue[3])
}))
write.csv(tab, "results/C_table.csv", row.names = FALSE)

C <- c(C, "", "Table, Holm-corrected p-values (direction where p < 0.05), uncorrected p-values and signs",
       "-----------------------------------",
       sprintf("%-7s %-30s %-5s  6h: %s%s  24h: %s%s  48h: %s%s   uncorrected %s (%+d), %s (%+d), %s (%+d)",
               tab$target, tab$drug, tab$expected,
               arrow(tab$sign6, tab$p6), fmt(tab$p6), arrow(tab$sign24, tab$p24), fmt(tab$p24),
               arrow(tab$sign48, tab$p48), fmt(tab$p48),
               fmt(tab$raw6), tab$sign6, fmt(tab$raw24), tab$sign24, fmt(tab$raw48), tab$sign48))

## recovered predictions
rec <- tab[tab$expected != "none" & ((tab$expected == "up" & ((tab$p24 < alpha & tab$sign24 > 0) | (tab$p48 < alpha & tab$sign48 > 0))) |
                                     (tab$expected == "down" & ((tab$p24 < alpha & tab$sign24 < 0) | (tab$p48 < alpha & tab$sign48 < 0)))), ]
pred_targets <- unique(tab$target[tab$expected != "none"])
C <- c(C, "", paste0("Predicted changes: ", length(pred_targets), " targets (", paste(pred_targets, collapse = ", "), ") in ",
                     sum(tab$expected != "none"), " drug and target pairs"),
       paste0("Recovered at 24 or 48 h with the predicted direction: ", length(unique(rec$target)), " targets (",
              paste(unique(rec$target), collapse = ", "), ") in ", nrow(rec), " pairs (",
              paste(paste(rec$target, rec$drug), collapse = ", "), ")"))

## specificity: significant single drugs per target and time point
sig <- d[!d$pair & d$p_holm < alpha, ]
C <- c(C, "", "Significant single drugs per target and time point (Holm over the set), with sign",
       "-----------------------------------")
for (tg in targets) for (h in exp_times$hours) {
  r <- sig[sig$target == tg & sig$hours == h, ]
  r <- r[order(r$p_holm), ]
  C <- c(C, sprintf("%-7s %2dh: %d%s", tg, h, nrow(r),
                    if (nrow(r)) paste0(", ", paste0(r$drug, " (", ifelse(r$coef > 0, "+", "-"), ", ", fmt(r$p_holm), ")", collapse = ", ")) else ""))
}

## set level, as the application reports it: minimum over the set and the time points per treatment
set_level <- tapply(d$p_holm, d$drug, min)
C <- c(C, "", paste0("Set level: significant single drugs ", sum(set_level[!is_pair[match(names(set_level), treatments$display)]] < alpha),
                     " of ", sum(!is_pair), ", significant pairs ", sum(set_level[is_pair[match(names(set_level), treatments$display)]] < alpha),
                     " of ", sum(is_pair)))
r <- sig[sig$target == "POLR2A" & sig$hours == 6, ]
C <- c(C, paste0("POLR2A at 6 h: ", nrow(r), " significant single drugs, ", sum(r$coef > 0), " up and ", sum(r$coef < 0), " down",
                 if (any(r$coef < 0)) paste0(" (", paste(r$drug[r$coef < 0], collapse = ", "), ")") else ""))
r <- sig[sig$target == "PARP1" & sig$hours == 6, ]
C <- c(C, paste0("PARP1 at 6 h: ", nrow(r), " significant single drugs, ", sum(r$coef > 0), " up: ",
                 paste(r$drug[order(r$p_holm)], collapse = ", ")))

## total effects for the predictions that the direct effects miss
tt <- dp[dp$effect_type == "total" & dp$hours > 6, ]
tt$p_holm <- p.adjust(tt$pvalue, "holm")
r <- tt[tt$target %in% c("TOP2B", "POLR2A") & tt$drug %in% c("Doxorubicin Hydrochloride", "Epirubicin Hydrochloride", "Lurbinectedin"), ]
r <- r[r$target == ifelse(r$drug == "Lurbinectedin", "POLR2A", "TOP2B"), ]
r <- r[order(r$target, r$drug, r$hours), ]
C <- c(C, "", "Total effects (Holm over the set at 24 and 48 h) for the two predictions the direct effects miss",
       "-----------------------------------",
       sprintf("%-7s %-26s %2dh: sign %+d, uncorrected %s, Holm %s", r$target, r$drug, r$hours, sign(r$coef), fmt(r$pvalue), fmt(r$p_holm)))

## drug pairs
ps <- d[d$pair & d$p_holm < alpha, ]
ps <- ps[order(ps$p_holm), ]
C <- c(C, "", "Significant drug pairs on the set (Holm over the set), with target, time point and sign",
       "-----------------------------------",
       sprintf("%-7s %2dh: %-55s %s %s", ps$target, ps$hours, ps$drug, ifelse(ps$coef > 0, "+", "-"), fmt(ps$p_holm)))
C <- c(C, paste0("Pairs: ", nrow(ps), " significant, ", sum(ps$hours == 6), " at 6 h, ", sum(ps$coef < 0), " negative, ",
                 sum(ps$target == "POLR2A"), " on POLR2A"))

## ---- network, as the application computes it ---------------------------------
np <- dbGetQuery(con, sprintf(
  "SELECT source, target, time_idx, pvalue FROM %s WHERE source IN (%s) OR target IN (%s)",
  pq("protein_pvalue"), inSel, inSel))
nTests <- (nTimes - 1) * M * (2 * nProt - M)
stopifnot(nrow(np) == nTests)
np$p_bh <- p.adjust(np$pvalue, "BH")
## stored protein coefficients are the positive ones, at the time index of the target
pc <- dbGetQuery(con, sprintf(
  "SELECT target, source, time_idx FROM %s WHERE source IN (%s) OR target IN (%s)",
  pq("protein_coef"), inSel, inSel))
np$positive <- paste(np$source, np$target, np$time_idx + 1) %in% paste(pc$source, pc$target, pc$time_idx)
np$from <- proteins$name[np$source]
np$to   <- proteins$name[np$target]
np$transition <- paste0(hours[as.character(np$time_idx)], "h to ", hours[as.character(np$time_idx + 1)], "h")

e <- np[np$p_bh < alpha, ]
self <- e[e$source == e$target, ]
between <- e[e$source != e$target & e$source %in% sel & e$target %in% sel, ]
C <- c(C, "", paste0("Network at FDR ", alpha, ": ", nrow(e), " edges, ", nrow(self), " self-dependencies of ",
                     length(unique(self$from)), " proteins (", paste(sort(unique(self$from)), collapse = ", "), "), ",
                     nrow(between), " edges between two queried proteins, ", nTests, " tests"))
C <- c(C, sprintf("Self-dependency %-7s %s: %s", self$from, self$transition, fmt(self$p_bh)))
r <- e[e$from == "TOP1" & e$to != "TOP1" & e$time_idx == 2, ]
r <- r[order(r$p_bh), ]
C <- c(C, paste0("TOP1 at 24 h to other proteins at 48 h: ",
                 paste0(r$to, " (", ifelse(r$positive, "+", "-"), ", ", fmt(r$p_bh), ")", collapse = ", ")))

ex <- np[np$p_bh < alpha_explore & np$from == "PARP1" & np$to != "PARP1" & np$time_idx == 1, ]
ex <- ex[order(ex$p_bh), ]
C <- c(C, "", paste0("Exploration at FDR ", alpha_explore, ", PARP1 at 6 h to other proteins at 24 h: ",
                     paste0(ex$to, " (", ifelse(ex$positive, "+", "-"), ", ", fmt(ex$p_bh), ")", collapse = ", ")))

dbDisconnect(con, shutdown = TRUE)

fileConn <- file("results/C_Results.txt")
writeLines(C, fileConn)
close(fileConn)
cat(C, sep = "\n")
