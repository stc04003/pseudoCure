library(paft)
library(geepack)
library(Rcpp)
library(microbenchmark)
library(pseudoCure)

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
b0 <- coef(glm(y ~ x1 + x2, data = dat))
pindex <- rep(0, ncol(X))

fit4 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "independence")
fit5 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "exchangeable")
fit6 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, pindex, corstr = "AR1")

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

## ##########################################################################################

lambda <- .5

fit4 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, c(0, 1, 1), corstr = "independence", tol = 1e-7)
fit5 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, c(0, 1, 1), corstr = "exchangeable", tol = 1e-7)
fit6 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, c(0, 1, 1), corstr = "AR1", tol = 1e-7)

fit7 <- gee(y, X, b0, nt, c(1, 0, 0), "identity", "ind", lambda, 1e-6, 1e-7, 50)
fit8 <- gee(y, X, b0, nt, c(1, 0, 0), "identity", "ex", lambda, 1e-6, 1e-7, 50)
fit9 <- gee(y, X, b0, nt, c(1, 0, 0), "identity", "ar1", lambda, 1e-6, 1e-7, 50)

cbind(fit4$b1, fit5$b1, fit6$b1)
cbind(fit7$b, fit8$b, fit9$b)

fit4$E
fit5$E
fit6$E
fit7$E
fit8$E
fit9$E

str(fit5)
str(fit8)

modifiedPGEE::PGEE(y ~ x1 + x2, id = id, data = dat, lambda = .5)$co
modifiedPGEE::PGEE(y ~ x1 + x2, id = id, data = dat, pindex = 1, lambda = .5)$co

fit4 <- paft:::pgeeCpp(X, y, weights, id, lambda, b0, c(0, 1, 1), corstr = "exchangeable",
                       maxit = 100, tol = 1e-7)
fit7 <- gee(y, X, b0, nt, c(1, 0, 0), "identity", "ex", lambda, 1e-6, 1e-7, 99)

str(fit4)
str(fit7)


(q_scad(abs(c(beta_new)), lambda)/(abs(as.vector(beta_new)) + eps))

e
## ##########################################################################################

sourceCpp(file = "PGEE.cpp")

library(modifiedPGEE)

cv1 <- function(lambda) {
    dat0 <- dat
    dat0$id <- rep(sample(unique(dat0$id)), tabulate(dat0$id))
    dat0 <- dat0[order(dat0$id),]      
    tmp <- CVfit(y ~ x1 + x2, id = id, data = dat0, fold = 5, lambda.vec = lambda,
                 eps = 1e-6, tol = 1e-7, maxiter = 50, family = gaussian("identity"))
    c(mean(tmp$cv.raw[[1]]), tmp$cvsd)
}

cv2 <- function(lambda) {
    tmp <- geeCV(y, X, b0, nt, c(1, 0, 0), "identity", "ind", 5, lambda, 1e-6, 1e-7, 50)
    c(mean(tmp), sd(tmp / 5))
}


summary(t(replicate(500, cv1(.8))))
summary(t(replicate(500, cv2(.8))))

miicrobenchmark(cv1(.8), cv2(.8))


tmp <- pseudoCure:::geeCV(y, X, b0, nt, c(1, 0, 0),
                          "identity", "scad", "ind", 5, 1:10 / 10, 1e-6, 1e-7, 50)
