library(Rcpp)
library(microbenchmark)

sourceCpp(file = "PGEE.cpp")

## ##################################################
## No cluster
## ##################################################
load("test.RData")
y <- test$y
mu <- test$mu
X <- test$X
Rhat <- test$Rhat
N <- test$N
nx <- test$nx
nt <- test$nt
index <- test$index
muEta <- test$muEta

out <- SHM(y, mu, X, Rhat, N, nx, nt, index, muEta)

## ##################################################
## With cluster
## ##################################################

load("test2.RData")
y <- test2$y
mu <- test2$mu
X <- test2$X
Rhat <- test2$Rhat
N <- test2$N
nx <- test2$nx
nt <- test2$nt
index <- test2$index
muEta <- test2$muEta

str(test2)

out1 <- SHM1(y, mu, X, Rhat, N, nx, nt, index, muEta)
out2 <- SHM2(y, mu, X, Rhat, N, nx, nt, index, muEta)
str(out1)
str(out2)

all.equal(out1, out2)
microbenchmark(out1, out2)
