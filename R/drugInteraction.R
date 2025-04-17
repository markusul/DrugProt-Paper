load("data/prepData.RData")
load("data/protNames.RData")

#mean per protein plate at time 0
dat0 <- data[data$pert_time == 0, ]
dat0 <- aggregate(dat0[, prot_names], by = list(dat0$protein_plate), FUN = mean)

# add mean at 0 to data
names(dat0)[-1] <- paste0(names(dat0)[-1], "_0")
row.names(dat0) <- dat0[, 1]
dat0 <- dat0[data$protein_plate, -1]

# add mean at 0 to data
dat <- data
dat <- cbind(dat, dat0)

# remove M453(ATCC) (no baseline protein expression)
dat <- dat[dat$protein_plate != "M453(ATCC)", ]

# add nois to imputation for computational und statistical stability
add_noise_to_imputation <- function(x){
  border <- sort(unique(x))[2]
  imp <- which(x < border)
  x_ <- x[-imp]
  y <- x
  y[imp] <- truncnorm::rtruncnorm(length(imp), a = -Inf, b = border, mean = mean(x_), sd = sd(x_))
  y
}

datI <- dat
plate_idx <- lapply(unique(dat$protein_plate), function(plate) which(dat$protein_plate == plate))

for(p in prot_names){
  for(plate in plate_idx)
    datI[plate, p] <- add_noise_to_imputation(datI[plate, p])
}

# protein of interest
P <- "P08621.P08621.RU17_HUMAN.SNRNP70.U1.small.nuclear.ribonucleoprotein.70.kDa"
t <- 6

# Data for model
Y <- datI[datI$pert_time == t, P] - datI[datI$pert_time == t, paste0(P, "_0")]
D <- datI[datI$pert_time == t, pert_names]
modelDat <- data.frame(Y = Y, D = D)

fit1 <- lm(Y ~ -1 + .^2, modelDat)
col <- c(rep("blue", ncol(D)), rep("red", choose(ncol(D), 2)))
plot(fit1$coefficients, col = col)

# debiased lasso
library(hdi)

## prepare design matrix with interactions
drug_design <- model.matrix(~ -1 + .^2, data = D)

# remove tratements wihtout data
with_data <- apply(drug_design, 2, function(x) length(unique(x))) != 1
drug_design <- drug_design[, with_data]

# add intercept to each treatment
drug_intercept <- drug_design
drug_intercept[drug_design != 0] <- 1
colnames(drug_intercept) <- paste0(colnames(drug_design), "_intercept")

# combine intercept and drug design
drug_design <- cbind(drug_design, drug_intercept)

#labels
single_effects <- rep(c(colnames(D), rep(NA, choose(ncol(D), 2))), 2)
length(single_effects[c(with_data, with_data)])
dim(drug_design)

fit2 <- lasso.proj(x = drug_design, y = Y, ncores = 11, parallel = TRUE)
plot(fit2$pval)
plot(fit2$pval.corr)

fit2

dev.off()
svd(D)
which(is.na(dat[dat$pert_time == t, paste0(P, "_0")]))
plot(drug_design[, '`drug_#56`:`drug_#64`'])
unique(dat$protein_plate)

dat[dat$pert_time == t, ][which(is.na(dat[dat$pert_time == t, paste0(P, "_0")])), "protein_plate"]
apply(D, 2, function(x) cor(x, Y))


plot(Y)
dev.off()


library(pemultinom)
install.packages("pemultinom")

X <- as.matrix(D)
str(X)
fit2 <- debiased_lasso(x = X, y = Y)

debiased_lasso(x = rnorm(100), y = rnorm(100))

library(DDL)
fit <- DDL(X = X, Y = Y, index = 3)


plot(D[1, ])
plot(colMeans(D))

library(hdi)
lasso.proj(x = D, y = Y)



Y <- jitter(Y)
Y <- rnorm(length(Y))
plot(Y)
plot(D[, 16], Y)
fit2 <- boot.lasso.proj(x = D, y = Y, ncores = 20, parallel = TRUE)
plot(fit2$pval)
plot(fit2$pval.corr)

fit2 <- hdi(x = D, y = Y)

library(truncnorm)
?rtruncnorm

x <- dat[dat$protein_plate == dat$protein_plate[1], P]



plot(dat[, 10], datI[, 10], col = as.factor(dat$protein_plate))

dat[, p]
y <- tapply(dat[, p], dat$protein_plate, add_noise_to_imputation)
dim(y)
length(y)
plot(dat[, p], y)
