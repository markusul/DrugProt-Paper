library(SDForest)

load("data/prepData.RData")
load("data/protNames.RData")

head(data)

X <- data[, prot_names]

#make UMAP
library(umap)
umap_res <- umap(X, n_neighbors = 50, min_dist = 0.5)

#plot UMAP
umap_res <- data.frame(umap_res$layout)
umap_res$IC50 <- data$IC50
umap_res$protein_plate <- data$protein_plate
umap_res$pert_time <- as.factor(data$pert_time)
umap_res$pertLabel <- data$pertLabel


library(ggplot2)
ggumap1 <- ggplot(umap_res, aes(x = X1, y = X2, color = pertLabel)) + 
    geom_point(size = 0.5) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    ggtitle("UMAP of Protein Expression Data")

ggsave("figures/umap.png", ggumap, width = 10, height = 6)



umap_0 <- umap(X[data$pert_time == 0, ], n_neighbors = 50, min_dist = 0.5)

#plot UMAP
umap_0 <- data.frame(umap_0$layout)
umap_0$IC50 <- data$IC50[data$pert_time == 0]
umap_0$protein_plate <- data$protein_plate[data$pert_time == 0]
umap_0$pert_time <- as.factor(data$pert_time[data$pert_time == 0])
umap_0$pertLabel <- data$pertLabel[data$pert_time == 0]

ggumap0 <- ggplot(umap_0, aes(x = X1, y = X2, color = protein_plate)) + 
    geom_point(size = 0.5) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    ggtitle("UMAP of Protein Expression Data at 0 hours")

ggumap0


umap_6 <- umap(X[data$pert_time == 6, ], n_neighbors = 50, min_dist = 0.5)

#plot UMAP
umap_6 <- data.frame(umap_6$layout)
umap_6$IC50 <- data$IC50[data$pert_time == 6]
umap_6$protein_plate <- data$protein_plate[data$pert_time == 6]
umap_6$pert_time <- as.factor(data$pert_time[data$pert_time == 6])
umap_6$pertLabel <- data$pertLabel[data$pert_time == 6]

ggumap6 <- ggplot(umap_6, aes(x = X1, y = X2, color = IC50)) + 
    geom_point(size = 0.5) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    ggtitle("UMAP of Protein Expression Data at 6 hours")

ggumap6