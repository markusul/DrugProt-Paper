print("start combination IC50")

# ---------------------------------------------------------------------------
# combinationIC50.R
#
# Replaces the combination IC50 values in data/ic50_values.csv by the values
# of the original anchored screen of Jaaks et al. (Nature 2022,
# doi:10.1038/s41586-022-04437-2), Figshare doi:10.6084/m9.figshare.16843597.
#
# SYNERGY_XMID_uM is the IC50 of the library drug in the presence of the
# anchor drug at ANCHOR_CONC, in uM. SYNERGY_XMID is the same value on the
# log2 normalised scale of the screen, with -1 to 9 the screening range, so
# log2(SYNERGY_XMID_uM) = SYNERGY_XMID - 9 + log2(LIBRARY_CONC).
#
# Inputs, in data/:
#   Original screen_All tissues_fitted data.csv   Figshare file of Jaaks et al.
#   ic50_values.csv                                current IC50 values
#   01_sample_info.xlsx                            sample annotation, for the drug names
#
# Output:
#   data/ic50_values_jaaks.csv   ic50_values.csv with IC50 in uM for all rows
#
# Run before data_preprocessing.R, which reads data/ic50_values_jaaks.csv.
# ---------------------------------------------------------------------------

jaaks <- read.csv("data/Original screen_All tissues_fitted data.csv", stringsAsFactors = FALSE)
jaaks <- jaaks[jaaks$Tissue == "Breast", ]

ic50 <- read.csv("data/ic50_values.csv", stringsAsFactors = FALSE)

# drug identifier to drug name, from the sample annotation, as in data_preprocessing.R
info <- readxl::read_excel("data/01_sample_info.xlsx", col_types = "text")
drugnames <- unique(rbind(data.frame(id = info$Pert_ID1, name = info$Pert_Name1),
                          data.frame(id = info$Pert_ID2, name = info$Pert_Name2)))
drugnames <- drugnames[!is.na(drugnames$id) & !duplicated(drugnames$id), ]
drug_lookup <- setNames(drugnames$name, drugnames$id)

# ---- name conventions --------------------------------------------------------

# cell lines differ in hyphens and spaces, e.g. Hs-578-T and HS578T
norm_cell <- function(x) toupper(gsub("[^A-Za-z0-9]", "", sub("\\(ATCC\\)", "", x)))

# drugs differ in salt forms, e.g. Lapatinib and Lapatinib Ditosylate Hydrate,
# so the first word is compared, with aliases for the remaining differences
drug_alias <- c("5fluorouracil" = "fluorouracil")
norm_drug <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]", "", sub(" .*", "", x)))
  ifelse(x %in% names(drug_alias), drug_alias[x], x)
}

# ---- Jaaks et al., one value per cell line, anchor, anchor dose and library --

# ANCHOR_CONC is read as text in the original screen, report entries that are not numbers
anchor_conc <- suppressWarnings(as.numeric(jaaks$ANCHOR_CONC))
if (any(is.na(anchor_conc))) {
  cat("non-numeric ANCHOR_CONC entries:\n")
  print(table(jaaks$ANCHOR_CONC[is.na(anchor_conc)], useNA = "ifany"))
}
jaaks$ANCHOR_CONC <- anchor_conc

# SYNERGY_XMID_uM is rounded to four decimals, which loses small IC50 values,
# so the uM scale is computed from the unrounded normalised value
jaaks <- jaaks[is.finite(jaaks$SYNERGY_XMID) & is.finite(jaaks$LIBRARY_CONC) &
               is.finite(jaaks$ANCHOR_CONC), ]
j <- data.frame(cell    = norm_cell(jaaks$CELL_LINE_NAME),
                anchor  = norm_drug(jaaks$ANCHOR_NAME),
                library = norm_drug(jaaks$LIBRARY_NAME),
                conc    = jaaks$ANCHOR_CONC,
                log2_uM = jaaks$SYNERGY_XMID - 9 + log2(jaaks$LIBRARY_CONC),
                xmid    = jaaks$SYNERGY_XMID)

# check of the conversion on values that the rounding leaves intact
big <- jaaks$SYNERGY_XMID_uM > 0.1
cat("conversion check, max abs difference on the log2 scale for IC50 > 0.1 uM:",
    round(max(abs(log2(jaaks$SYNERGY_XMID_uM[big]) - j$log2_uM[big])), 4), "\n")

# replicate measurements are summarised by their median
key_j <- paste(j$cell, j$anchor, j$library, signif(j$conc, 3))
j <- data.frame(key     = unique(key_j),
                log2_uM = tapply(j$log2_uM, key_j, median)[unique(key_j)],
                xmid    = tapply(j$xmid, key_j, median)[unique(key_j)],
                n       = as.vector(table(key_j)[unique(key_j)]))

# ---- combinations of this analysis -------------------------------------------

cb  <- grepl(" ", ic50$pertLabel)
ids <- do.call(rbind, strsplit(ic50$pertLabel[cb], " "))
w <- data.frame(row     = which(cb),
                cell    = norm_cell(ic50$cell_line[cb]),
                anchor  = norm_drug(drug_lookup[ids[, 1]]),
                library = norm_drug(drug_lookup[ids[, 2]]),
                conc    = ic50$Anchor_dose[cb])

m_fwd <- match(paste(w$cell, w$anchor, w$library, signif(w$conc, 3)), j$key)

# combinations in which the first drug of pertLabel is the library drug of the screen
m_swp <- match(paste(w$cell, w$library, w$anchor, signif(w$conc, 3)), j$key)
m <- ifelse(is.na(m_fwd), m_swp, m_fwd)

# ---- diagnostics ---------------------------------------------------------------

cat("combination rows matched:", sum(!is.na(m)), "of", length(m),
    "(", sum(!is.na(m_fwd)), "as anchor and library,",
    sum(is.na(m_fwd) & !is.na(m_swp)), "with the roles swapped )\n")
if (any(is.na(m_fwd) & !is.na(m_swp))) {
  cat("combinations matched with the roles swapped:\n")
  print(unique(w[is.na(m_fwd) & !is.na(m_swp), c("anchor", "library", "conc")]))
}

# unmatched combinations that the screen holds for other cell lines
pair_conc <- sub("^[^ ]+ ", "", j$key)
w_pc <- paste(w$anchor, w$library, signif(w$conc, 3))
w_pc_swp <- paste(w$library, w$anchor, signif(w$conc, 3))
cat("unmatched rows whose combination and dose occur in other cell lines:",
    sum(is.na(m) & (w_pc %in% pair_conc | w_pc_swp %in% pair_conc)), "\n")
cat("cell lines not in the screen:",
    paste(setdiff(unique(w$cell), sub(" .*", "", j$key)), collapse = ", "), "\n")
jdrugs <- unique(unlist(lapply(strsplit(j$key, " "), `[`, 2:3)))
cat("drugs not in the screen:",
    paste(setdiff(unique(c(w$anchor, w$library)), jdrugs), collapse = ", "), "\n")
if (any(is.na(m))) {
  cat("unmatched anchor, library and anchor dose:\n")
  print(unique(w[is.na(m), c("anchor", "library", "conc")]))
}

# what the old file holds, compared on the matched rows, skipped if the file holds no combination values
old <- suppressWarnings(as.numeric(ic50$IC50[w$row]))
ok  <- !is.na(m) & is.finite(old)
if (sum(ok) > 1) for (s in c("xmid", "log2_uM")) {
  d <- old[ok] - j[[s]][m[ok]]
  cat("old value against", s, ": correlation", round(cor(old[ok], j[[s]][m[ok]]), 3),
      ", share equal within 0.01:", round(mean(abs(d) < 0.01), 3), ", differences:\n")
  print(round(quantile(d, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)), 3))
}

# ---- write -------------------------------------------------------------------

# IC50 in uM for all rows, single drugs are already untransformed uM
ic50$IC50[w$row] <- 2^j$log2_uM[m]
write.csv(ic50, "data/ic50_values_jaaks.csv", row.names = FALSE)
cat("written data/ic50_values_jaaks.csv\n")
