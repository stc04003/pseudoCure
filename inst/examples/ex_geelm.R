gendat <- function() {
  id <- gl(50, 4, 200)
  visit <- rep(1:4, 50)
  x1 <- rbinom(200, 1, 0.6)
  x2 <- runif(200, 0, 1)
  phi <- 1 + 2 * x1
  rhomat <- 0.667^outer(1:4, 1:4, function(x, y) abs(x - y))
  chol.u <- chol(rhomat)
  noise <- as.vector(sapply(1:50, function(x) chol.u %*% rnorm(4)))
  e <- sqrt(phi) * noise
  y <- 1 + 3 * x1 - 2 * x2 + e
  dat <- data.frame(y, id, visit, x1, x2)
  dat
}

set.seed(1); str(dat <- gendat())
geelm(y ~ x1 + x2, id = id, data = dat, corstr = "ar1")
