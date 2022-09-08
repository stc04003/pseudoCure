library(Rcpp)
library(microbenchmark)

## Assume family = gaussian for our cure models
## Assume fi = 1 because of the gaussian assumption
sourceCpp(file = "PGEE.cpp")

## ##################################################
## Generate example data
## ##################################################

library(mvtnorm)
n <- 2000; pn <- 10; m <- 4
id.vect <- rep(1:n, each = m) 
X.sigma <- matrix(0,(pn-1),(pn-1))
for (i in 1:(pn-1))
  X.sigma[i,] <- 0.5^(abs((1:(pn-1))-i))  
x.mat <- as.matrix(rmvnorm(n*m, mean = rep(0,(pn-1)), X.sigma))
x.mat <- cbind(rbinom(n*m,1, 0.5), x.mat)
beta.true <- c(2,3,1.5,2,rep(0,6))
sigma2 <- 1
rho <- 0.5
R <- matrix(rho,m,m)+diag(rep(1-rho,m))
SIGMA <- sigma2*R
error <- rmvnorm(n, mean = rep(0,m),SIGMA)
y.temp <- x.mat%*%beta.true
y.vect <- y.temp+as.vector(t(error))    
mydata <- data.frame(id.vect,y.vect,x.mat) 
colnames(mydata) <- c("id","y",paste("x",1:length(beta.true),sep = ""))

formula <- "y ~.-id-1"
data <- mydata
family <- gaussian(link = "identity")
lambda.vec <- seq(0.1,1,0.1)
lambda <- 0.1    

myfit1 <- PGEE(formula = formula, id = id, data = data, na.action = NULL, 
               family = family, corstr = "exchangeable", Mv = NULL, 
               beta_int = c(rep(0,length(beta.true))), R = NULL, scale.fix = TRUE, 
               scale.value = 1, lambda = lambda, pindex = NULL, eps = 10^-6, maxiter = 30, 
               tol = 10^-3, silent = TRUE)

## ##################################################
## Test
## ##################################################

load("testSHEM.RData")
for (i in 1:length(test))
  eval(parse(text = paste(names(test)[i], " <- test$", names(test)[i], sep = "")))

library(modifiedPGEE)
debug(S_H_E_M)
out <- S_H_E_M(N, nt, y, X, nx, family, beta_new, Rhat, fihat, lambda, pindex, eps = 10^-6) 
str(out)

pindex <- rep(0, ncol(X))
out2 <- SHEM(y, X, beta_new, Rhat, nt, pindex, lambda, eps)
str(out2)

sourceCpp(file = "PGEE.cpp")
out <- gee(y, X, beta_new, Rhat, nt, pindex, lambda, eps, 1e-7, 100)
str(out)

e






out1 <- SHM1(y, mu, X, Rhat, N, nx, nt, index, muEta)
out2 <- SHM(y, mu, X, Rhat, N, nx, nt, index, muEta)
str(out1)
str(out2)

all.equal(out1, out2)
microbenchmark(SHM1(y, mu, X, Rhat, N, nx, nt, index, muEta),
               SHM(y, mu, X, Rhat, N, nx, nt, index, muEta))


sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
arma::mat r3(arma::mat A, arma::mat B, arma::mat C) {
  return(A.t() * B * C * C.t() * B * A);
}
// [[Rcpp::export]]
arma::mat r4(arma::mat A, arma::mat B, arma::mat C) {
  arma::mat tmp = A.t() * B * C;
  return(tmp * tmp.t());
}')

n <- 100
A1 <- matrix(runif(n * n), n, n)
A2 <- matrix(runif(n * n), n, n)
A3 <- matrix(runif(n * n), n, n)

all.equal(r3(A1, A2, A3), r4(A1, A2, A3))
microbenchmark(r3(A1, A2, A3), r4(A1, A2, A3))


sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
arma::vec qscad(arma::vec b, double lambda) {
  double a = 3.7;
  arma::vec dif = a * lambda - b;
  dif.elem(find(dif < 0)).zeros();
  arma::vec out(b.n_elem);
  out.fill(lambda);
  arma::uvec ind = find(b > lambda);
  out(ind) = dif(ind) / (a - 1);
  return(out);
}')

r3(1:4)
cumsum(1:3)

qscad(abs(beta_new), lambda)
q_scad(abs(c(beta_new)), lambda)

microbenchmark(qscad(abs(beta_new), lambda),
               q_scad(abs(c(beta_new)), lambda))

sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
Rcpp::List r2(double a) {
  Rcpp::List out(4);
  arma::mat M1(1, 3, arma::fill::zeros);
  arma::mat M2(2, 3, arma::fill::zeros);
  arma::mat M3(3, 3, arma::fill::zeros);
  arma::mat M4(4, 3, arma::fill::zeros);
  out(0) = M1;
  out(1) = M2;
  out(2) = M3;
  out(3) = M4;
  return(out);
}
// [[Rcpp::export]]
arma::vec r3(double a) {
  arma::vec out(10);
  out.fill(a);
  Rcpp::List tmp = r2(a);
  arma::mat m = tmp(1);
  std::cout << m;
  return(out);
}')

r3(3)


q_scad(abs(c(beta_new)), lambda) / (abs(as.vector(beta_new)) + eps)
