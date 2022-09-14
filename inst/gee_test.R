library(geepack)
library(microbenchmark)
library(pseudoCure)

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

## With pseudoCure()
y <- dat$y
X <- model.matrix(~ x1 + x2, data = dat)
id <- as.numeric(dat$id)
weights <- rep(1, length(y))
lambda <- 0
b0 <- coef(glm(y ~ x1 + x2, data = dat))
pindex <- rep(0, ncol(X))
nt <- tabulate(id)

fit4 <- pseudoCure:::gee(y, X, b0, nt, "identity", "ind", 1e-7, 50)
fit5 <- pseudoCure:::gee(y, X, b0, nt, "identity", "ex", 1e-7, 50)
fit6 <- pseudoCure:::gee(y, X, b0, nt, "identity", "ar1", 1e-7, 50)

cbind(fit4$b, fit5$b, fit6$b)

## Our implementations are all faster than geepack's counterparts
microbenchmark(fit1 = geese(y ~ x1 + x2, id = id, data = dat, corstr = "ind"),
               fit2 = geese(y ~ x1 + x2, id = id, data = dat, corstr = "ex"),
               fit3 = geese(y ~ x1 + x2, id = id, data = dat, corstr = "ar1"),
               fit4 = pseudoCure:::gee(y, X, b0, nt, "identity", "ind", 1e-7, 50),
               fit5 = pseudoCure:::gee(y, X, b0, nt, "identity", "ex", 1e-7, 50),
               fit6 = pseudoCure:::gee(y, X, b0, nt, "identity", "ar1", 1e-7, 50))

               
## Compare all; difference < 5e-5
all.equal(cbind(fit1$beta, fit2$beta, fit3$beta),
          cbind(fit4$b, fit5$b, fit6$b))
