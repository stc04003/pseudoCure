## Function to generate simulated data under the PHMC model
simMC <- function(n) {
  p <- 10
  a <- c(1, 0, -1, 0, 0, 0, 0, 0, 0, 0) # incidence coefs.
  b <- c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0) # latency coefs.
  X <- data.frame(x = matrix(runif(n * p), n))
  X$x.3 <- 1 * (X$x.3 > .5)
  X$x.4 <- 1 * (X$x.4 > .5)
  X[,5:10] <- apply(X[,5:10], 2, qnorm)  
  time <- -3 * exp(-colSums(b * t(X))) * log(runif(n))
  cure.prob <- 1 / (1 + exp(-2 - colSums(a * t(X))))
  Y <- rbinom(n, 1, cure.prob) 
  cen <- rexp(n, .02)
  dat <- NULL  
  dat$Time <- pmin(time / Y, cen)
  dat$Status <- 1 * (dat$Time == time)
  data.frame(dat, X)
}

## Fix seed and generate data
set.seed(1); datMC <- simMC(1000)

## Oracle model with an unpenalized PHMC model
summary(fit1 <- pCure(~ x.1 + x.3, ~ x.1 + x.3, Time, Status, datMC))

## Penalized PHMC model with tuning parameters selected by 10-fold cross validation
## User specifies the range of tuning parameters
summary(update(fit1, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200))
## Auto selection of the range of tuning parameters
summary(update(fit1, lambda1 = "auto", lambda2 = "auto"))

## Penalized PHMC model given tuning parameters
summary(update(fit1, lambda1 = 0.006, lambda2 = 0.022))
