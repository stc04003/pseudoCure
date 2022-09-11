library(paft)
library(geepack)
library(Rcpp)
library(microbenchmark)

sourceCpp(file = "PGEE.cpp")

## Example from ?geese
gendat <- function() {
    id <- gl(50, 4, 200)
    visit <- rep(1:4, 50)
    x1 <- rbinom(200, 1, 0.6) ## within cluster varying binary covariate
    x2 <- runif(200, 0, 1)   ## within cluster varying continuous covariate
    phi <- 1 + 2 * x1         ## true scale model
    ## the true correlation coefficient rho for an ar(1)
    ## correlation structure is 0.667.
    rhomat <- 0.667 ^ outer(1:4, 1:4, function(x, y) abs(x - y))
    chol.u <- chol(rhomat)
    noise <- as.vector(sapply(1:50, function(x) chol.u %*% rnorm(4)))
    e <- sqrt(phi) * noise
    y <- 1 + 3 * x1 - 2 * x2 + e
    dat <- data.frame(y, id, visit, x1, x2)
    dat
}

set.seed(1); dat <- gendat()

fit1 <- geese(y ~ x1 + x2, id = id, data = dat, corstr = "ind")
fit2 <- geese(y ~ x1 + x2, id = id, data = dat, corstr = "ex")
fit3 <- geese(y ~ x1 + x2, id = id, data = dat, corstr = "ar1")

cbind(fit1$beta, fit2$beta, fit3$beta)

## paft returns zero? 
y <- dat$y
X <- model.matrix(~ x1 + x2, data = dat)
id <- as.numeric(dat$id)
weights <- rep(1, length(y))
lambda <- 0
b0 <- rep(0, ncol(X))
pindex <- rep(0, ncol(X))

fit4 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "independence")
fit5 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "exchangeable")
fit6 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "AR1")

cbind(fit4$b1, fit6$b1)
cbind(fit4$b1, fit5$b1, fit6$b1)

## new gee
nt <- tabulate(id)

fit7 <- gee(y, X, b0, nt, rep(1, ncol(X)), "identity", "ind", 0, 1e-6, 1e-7, 50)
fit8 <- gee(y, X, b0, nt, rep(1, ncol(X)), "identity", "ex", 0, 1e-6, 1e-7, 50)
fit9 <- gee(y, X, b0, nt, rep(1, ncol(X)), "identity", "ar1", 0, 1e-6, 1e-7, 50)

fit72 <- gee0(y, X, b0, nt, "identity", "ind", 1e-7, 50)
fit82 <- gee0(y, X, b0, nt, "identity", "ex", 1e-7, 50)
fit92 <- gee0(y, X, b0, nt, "identity", "ar1", 1e-7, 50)

identical(fit7[-4], fit72)
identical(fit8[-4], fit82)
identical(fit9[-4], fit92)

cbind(fit7$b, fit8$b, fit9$b)

microbenchmark(i = gee(y, X, b0, nt, rep(1, ncol(X)), "identity", "ar1", 0, 1e-6, 1e-7, 50),
               ii = gee0(y, X, b0, nt, "identity", "ar1", 1e-7, 50))               

## Compare all
cbind(fit1$beta, fit2$beta, fit3$beta)
cbind(fit4$b1, fit6$b1)
cbind(fit7$b, fit8$b, fit9$b)



sourceCpp(file = "PGEE.cpp")

geeCV(y, X, b0, nt, rep(1, ncol(X)), "identity", "ex", 5, .1, 1e-6, 1e-7, 50)

library(modifiedPGEE)
CVfit(y ~ x1 + x2, id = id, data = dat, fold = 5, lambda.vec = .1)
