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

# protein of interest
P <- "P08621.P08621.RU17_HUMAN.SNRNP70.U1.small.nuclear.ribonucleoprotein.70.kDa"
t <- 6

# Data for model
Y <- dat[dat$pert_time == t, P] - dat[dat$pert_time == t, paste0(P, "_0")]
D <- dat[dat$pert_time == t, pert_names]
modelDat <- data.frame(Y = Y, D = D)

fit1 <- lm(Y ~ -1 + .^2, modelDat)
col <- c(rep("blue", ncol(D)), rep("red", choose(ncol(D), 2)))
plot(fit1$coefficients, col = col)

# debiased lasso
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

add_noise_to_imputation <- function(x){
    x <- dat[dat$protein_plate == dat$protein_plate[1], P]
    border <- sort(unique(x))[2]
    imp <- which(x < border)
    x_ <- x[-imp]
    y <- x
    y[imp] <- truncnorm::rtruncnorm(length(imp), a = -Inf, b = border, mean = mean(x_), sd = sd(x_))
    Y
}

plot(x)
points(y, col = 3)

plot(x, y)
