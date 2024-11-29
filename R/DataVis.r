library(SDForest)
library(ggplot2)
library(umap)

load("data/prepData.RData")
load("data/protNames.RData")


# Visualization of Protein expression
X <- data[, prot_names]
umap_res <- umap(X, n_neighbors = 50, min_dist = 0.5)

#plot UMAP
umap_res <- data.frame(umap_res$layout)
umap_res$IC50 <- data$IC50
umap_res$protein_plate <- data$protein_plate
umap_res$pert_time <- as.factor(data$pert_time)
umap_res$pertLabel <- data$pertLabel

ggumap <- ggplot(umap_res, aes(x = X1, y = X2, color = pertLabel)) + 
    geom_point(size = 0.5) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    ggtitle("UMAP of Protein Expression Data")

ggumap
ggsave("figures/umap.png", ggumap, width = 10, height = 6)


umap_0 <- umap(X[data$pert_time == 0, ], n_neighbors = 50, min_dist = 0.5)
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
        axis.ticks.y=element_blank())
ggumap0
ggsave("figures/umap0.png", ggumap0, width = 5, height = 4)


umap_6 <- umap(X[data$pert_time == 6, ], n_neighbors = 50, min_dist = 0.5)
umap_6 <- data.frame(umap_6$layout)
umap_6$IC50 <- data$IC50[data$pert_time == 6]
umap_6$protein_plate <- data$protein_plate[data$pert_time == 6]
umap_6$pert_time <- as.factor(data$pert_time[data$pert_time == 6])
umap_6$pertLabel <- data$pertLabel[data$pert_time == 6]

ggumap6 <- ggplot(umap_6, aes(x = X1, y = X2, color = IC50)) + 
    geom_point(size = 0.1) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggumap6
ggsave("figures/umap6.png", ggumap6, width = 5, height = 4)

umap_24 <- umap(X[data$pert_time == 24, ], n_neighbors = 50, min_dist = 0.5)
umap_24 <- data.frame(umap_24$layout)
umap_24$IC50 <- data$IC50[data$pert_time == 24]
umap_24$protein_plate <- data$protein_plate[data$pert_time == 24]
umap_24$pert_time <- as.factor(data$pert_time[data$pert_time == 24])
umap_24$pertLabel <- data$pertLabel[data$pert_time == 24]

ggumap24 <- ggplot(umap_24, aes(x = X1, y = X2, color = IC50)) + 
    geom_point(size = 0.1) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggumap24
ggsave("figures/umap24.png", ggumap24, width = 5, height = 4)


umap_48 <- umap(X[data$pert_time == 48, ], n_neighbors = 50, min_dist = 0.5)
umap_48 <- data.frame(umap_48$layout)
umap_48$IC50 <- data$IC50[data$pert_time == 48]
umap_48$protein_plate <- data$protein_plate[data$pert_time == 48]
umap_48$pert_time <- as.factor(data$pert_time[data$pert_time == 48])
umap_48$pertLabel <- data$pertLabel[data$pert_time == 48]

ggumap48 <- ggplot(umap_48, aes(x = X1, y = X2, color = IC50)) + 
    geom_point(size = 0.1) + theme_bw() + 
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
    theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggumap48
ggsave("figures/umap48.png", ggumap48, width = 5, height = 4)


# Information about data
load("data/aggData.RData")
load("data/combData.RData")

# number of protein plates
length(unique(agg_data$protein_plate))

#number of drug combinations
length(unique(comb_data$pertLabel)) - 63

#concentrations of different drugs
table(unlist(data[, pert_names]))

#number of experiments with a given drug
sort(table(data$pertLabel))

#example of concentrations in drug combinations
gg_conc <- ggplot(data, aes(x = `drug_#20`, y = `drug_#32`)) + 
  geom_point(size = 0.7) + theme_bw() + 
  xlab(expression(paste(mu, "mol Olaparib"))) +
  ylab(expression(paste(mu, "mol Lapatinib Ditosylate Hydrate")))
gg_conc

ggsave('figures/conc.png', gg_conc, width = 4, height = 3)

#variation of IC50
length(unique(data$IC50))
length(unique(comb_data$pertLabel)) * 18
dim(unique(data[, c('pertLabel', 'protein_plate')]))

gg_IC50 <- ggplot(data, aes(x = type, y = IC50)) + 
  geom_violin() + theme_bw() + geom_boxplot(width = 0.07, size = 0.5)
gg_IC50

ggsave('figures/IC50.png', gg_IC50, width = 4, height = 3)


#number of proteins
length(prot_names)
