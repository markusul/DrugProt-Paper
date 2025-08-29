load("results/DrugEffects.RData")

# clustering of proteins according to drug effects
d <- dist(allPvecs)
fit <- hclust(d)
plot(fit)
clusters <- cutree(fit, 25)

fit <- kmeans(allPvecs, centers = 10)
fit

res <- umap::umap(allPvecs)
res$layout

plot(res$layout, col = clusters)

