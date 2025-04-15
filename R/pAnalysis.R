library(readr)
pDat <- read_csv("data/p-values_matrix.csv")

pDat <- pDat[pDat$Treatment != "divider", ]

Dcomb <- unique(pDat$Treatment)
Dcomb <- Dcomb[nchar(Dcomb) > 2]

Dcomb
pDat_agg <- lapply(Dcomb, function(D) {
  single <- strsplit(D, ':')[[1]]
  out <- pDat[pDat$Treatment == D, 1:3]
  out$Combination <- -log(out$Coefficient)
  First <- -log(pDat[pDat$Treatment == single[1], 'Coefficient'])
  names(First) <- "First"
  Second <- -log(pDat[pDat$Treatment == single[2], 'Coefficient'])
  names(Second) <- "Second"
  
  cbind(out, First, Second)
})


pDat_agg <- do.call(rbind, pDat_agg)
pDat_agg$min <- pmin(pDat_agg$First, pDat_agg$Second)
pDat_agg$max <- pmax(pDat_agg$First, pDat_agg$Second)
pDat_agg$Protein <- as.factor(pDat_agg$Protein)

pDat_agg$Protein[which(pDat_agg$Coefficient %in% sort(pDat_agg$Coefficient, decreasing = T)[1:4])]

pDat_top <- pDat_agg[pDat_agg$Combination > 9, ]

plot(Combination ~ min, pDat_agg, col = Protein)
text(x = pDat_top$min, y = pDat_top$Combination, labels = pDat_top$Treatment)
title('-log(p-values)')
plot(Combination ~ max, pDat_agg, col = Protein)
text(x = pDat_top$max, y = pDat_top$Combination, labels = pDat_top$Treatment)
title('-log(p-values)')

library(plotly)

fig <- plot_ly(pDat_agg, z = ~Combination, x = ~First, y = ~Second, color = ~Protein, 
               size = 2)
fig %>% add_trace(text = ~Treatment)
