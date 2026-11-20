# Short protein names for the figures, the database and the application.
#
# An identifier reads ACC.ACC.ENTRY_HUMAN.GENE.Description, and a protein group repeats each
# field once per member. Every character of the original names other than letters, digits and
# underscores is a period in the identifier, so the hyphen of a gene symbol is lost and is restored here.
short_names <- function(prot_names) {
  # gene symbols with these prefixes carry a hyphen, such as H1-2, H3-3A, HLA-A and MT-CO2
  hyphenated <- c("H1", "H3", "HLA", "MT")
  parts <- strsplit(prot_names, "[.]")
  accession <- sapply(parts, `[`, 1)

  short <- sapply(parts, function(p) {
    human <- grep("_HUMAN", p)
    n     <- length(human)
    rest  <- p[-seq_len(max(human))]   # the gene symbols, followed by the description
    if (n == 1) {
      if (rest[1] %in% hyphenated) return(paste(rest[1], rest[2], sep = "-"))
      return(rest[1])
    }
    # the members of a group can share a gene symbol, then fewer than n symbols are listed
    # and a capitalised word marks the start of the description
    genes <- rest[seq_len(n)]
    word  <- grep("^[A-Z][a-z]", genes)
    if (length(word) > 0) genes <- genes[seq_len(word[1] - 1)]
    paste(genes, collapse = "/")
  })

  # a protein without a gene symbol is named by its UniProt accession
  missing <- short == "nan"
  short[missing] <- accession[missing]

  # proteins that share a gene symbol carry the accession in addition
  shared <- short %in% short[duplicated(short)]
  short[shared] <- paste0(short[shared], " (", accession[shared], ")")

  # the names are unique, and the figures split node names at the underscore
  stopifnot(!anyDuplicated(short), all(nzchar(short)), !grepl("_", short))
  setNames(short, prot_names)
}
