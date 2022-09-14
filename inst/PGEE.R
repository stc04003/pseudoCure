library(Rcpp)
library(microbenchmark)

## incidence: gaussian(link = "logit")
## latency: gaussian(link = "cloglog")
## long-term: gaussian(link = "log")
## short-term: gaussian(link = "cloglog")
## Assume fi = 1 because of the gaussian assumption
## Assume equal cluster size
sourceCpp(file = "PGEE.cpp")

## ##################################################
## Generate example data
## ##################################################

set.seed(1)

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
family <- gaussian(link = "log")
lambda.vec <- seq(0.1,1,0.1)
lambda <- 0.1
lambda <- 0


## ##################################################
library(modifiedPGEE)
library(geepack)
library(Rcpp)
sourceCpp(file = "PGEE.cpp")

y <- data$y
X <- as.matrix(data[,3:12])

lambda <- 0

system.time(e1 <- PGEE(y ~ . -id - 1, id = id, data = data, na.action = NULL, 
                       family = gaussian(link = 'log'),
                       corstr = "independence", Mv = NULL,
                       ## corstr = "independence", Mv = NULL, 
                       beta_int = c(rep(0,length(beta.true))), R = NULL, scale.fix = TRUE, 
                       scale.value = 1, lambda = 0, pindex = 1:10,
                       eps = 1e-6, maxiter = 30, 
                       tol = 1e-7, silent = TRUE))
system.time(e2 <- PGEE(y ~ . -id - 1, id = id, data = data, na.action = NULL, 
                       family = gaussian(link = 'log'),
                       corstr = "exchangeable", Mv = NULL,
                       ## corstr = "independence", Mv = NULL, 
                       beta_int = c(rep(0,length(beta.true))), R = NULL, scale.fix = TRUE, 
                       scale.value = 1, lambda = 0, pindex = 1:10,
                       eps = 1e-6, maxiter = 30, 
                       tol = 1e-7, silent = TRUE))
system.time(e3 <- PGEE(y ~ . -id - 1, id = id, data = data, na.action = NULL, 
                       family = gaussian(link = 'log'),
                       corstr = "AR-1", Mv = NULL,
                       ## corstr = "independence", Mv = NULL, 
                       beta_int = c(rep(0,length(beta.true))), R = NULL, scale.fix = TRUE, 
                       scale.value = 1, lambda = 0, pindex = 1:10,
                       eps = 1e-6, maxiter = 30, 
                       tol = 1e-7, silent = TRUE))

system.time(f1 <- geese(y ~ . -id - 1, id = id, data = data, corstr = "ind",
                        family = gaussian(link = 'log'), scale.fix = T, 
                        b = rep(0, 10)))
system.time(f2 <- geese(y ~ . -id - 1, id = id, data = data, corstr = "ex",
                        family = gaussian(link = 'log'), scale.fix = T, 
                        b = rep(0, 10)))
system.time(f3 <- geese(y ~ . -id - 1, id = id, data = data, corstr = "ar1",
                        family = gaussian(link = 'log'), scale.fix = T, 
                        b = rep(0, 10)))

system.time(g1 <- gee(y, X, rep(0, 10), tabulate(data$id), rep(1, 10),
                      "log", "ind", lambda, 1e-6, 1e-7, 30))
system.time(g2 <- gee(y, X, rep(0, 10), tabulate(data$id), rep(1, 10),
                      "log", "ex", lambda, 1e-6, 1e-7, 30))
system.time(g3 <- gee(y, X, rep(0, 10), tabulate(data$id), rep(1, 10),
                      "log", "ar1", lambda, 1e-6, 1e-7, 30))

e1$iter
e2$iter
e3$iter
g1$iter
g2$iter
g3$iter

cbind(drop(e1$coefficients), drop(e2$coefficients), drop(e3$coefficients))
cbind(drop(f1$b), drop(f2$b), drop(f3$b))
cbind(drop(g1$b), drop(g2$b), drop(g3$b))

gaussian()$variance
gaussian('logit')$variance
gaussian('log')$variance
gaussian('cloglog')$variance

e

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
pindex <- rep(0, ncol(X))
out <- gee(y, X, beta_new, nt, pindex, "log", "ex", lambda, 1e-6, 1e-7, 100)
str(out)

system.time(gee(y, X, beta_new, nt, pindex, "log", "ind", lambda, 1e-6, 1e-7, 100))
system.time(gee(y, X, beta_new, nt, pindex, "log", "ex", lambda, 1e-6, 1e-7, 100))

str(gee(y, X, beta_new, nt, pindex, "log", "ind", lambda, 1e-6, 1e-7, 100))
str(gee(y, X, beta_new, nt, pindex, "log", "ex", lambda, 1e-6, 1e-7, 100))


system.time(gee(y, X, beta_new, nt, rep(1, 10), "log", "ex", lambda, 1e-6, 1e-7, 100))

str(gee(y, X, beta_new, nt, rep(1, 10), "log", "ex", 0, 1e-6, 1e-7, 100))

library(geepack)
system.time(f1 <- geese(y ~ X - 1, id = rep(1:2000, nt), family = gaussian("log"),
                        scale.fix = TRUE, b = c(beta_new)))
system.time(f2 <- geese(y ~ X - 1, id = rep(1:2000, nt), family = gaussian("log"),
                        corstr = "ex", scale.fix = TRUE, b = c(beta_new)))
f1$b
f2$b

geese(y ~ X - 1, id = 1:8000, family = gaussian("log"), scale.fix = TRUE, b = c(beta_new))$b

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
// [[Rcpp::export]]
Rcpp::List r3(arma::vec eta) {
	Rcpp::List out(2);
	arma::vec tmp = eta;	
	tmp.elem(find(tmp > 700)).fill(700);
  tmp = exp(tmp) % exp(-exp(tmp));
	tmp.replace(arma::datum::inf, pow(2, 1023));
	out(0) = 1 - exp(-exp(eta));
	out(1) = tmp;
	return out;
}')

r3(1:3)


q_scad(abs(c(beta_new)), lambda) / (abs(as.vector(beta_new)) + eps)



sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
double a1(arma::vec a) {
  arma::mat tmp = a * a.t(); 
  tmp.diag().zeros();
  return(sum(sum(tmp)));
}
// [[Rcpp::export]]
double a2(arma::vec a) {
  double out = 0; 
  int n = a.n_elem;
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
       out += a[i] * a[j];
  }}
  return(2 * out);
}
// [[Rcpp::export]]
double a22(arma::vec a) {
  double out = 0; 
  int n = a.n_elem;
  double a2 = sum(a) * sum(a);
  return(a2 - sum(a % a));
}')

a1(1:50)
a2(1:50)
a22(1:50)

a <- 1:5
microbenchmark(tcrossprod(a), a1(a), a2(a))


res <- drop(y - X %*% beta_new)
res[1:4]
res[5:8]

sum(tcrossprod(res[1:4]))

Rhat[,,1]

sum(sapply(1:N, function(i) {
    tmp <- tcrossprod(res[1:4 + 4 * i - 4])
    diag(tmp) <- 0
    sum(tmp)})) / 12 / N
sum(sapply(1:N, function(i) a1(res[1:4 + 4 * i - 4]))) / 12 / N
sum(sapply(1:N, function(i) a2(res[1:4 + 4 * i - 4]))) / 12 / N

microbenchmark(sum(sapply(1:N, function(i) {
    tmp <- tcrossprod(res[1:4 + 4 * i - 4])
    diag(tmp) <- 0
    sum(tmp)})) / 12 / N,
    sum(apply(matrix(res, 4), 2, a2)) / 12 / N,
    sum(sapply(1:N, function(i) a1(res[1:4 + 4 * i - 4]))) / 12 / N,
    sum(sapply(1:N, function(i) a2(res[1:4 + 4 * i - 4]))) / 12 / N)


sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
double a3(arma::vec a) {
  double out = 0; 
  int n = a.n_elem;
  for (int i = 0; i < n - 1; i++) {
    out += a[i] * a[i + 1];
  }
  return(out);
}
// [[Rcpp::export]]
double a32(arma::vec a) {
  int n = a.n_elem;
  double out = sum(a(span(0, n - 2)) % a(span(1, n - 1)));
  return(out);
}')

a3(1:50)
a32(1:50)

sum(sapply(1:N, function(i) a3(res[1:4 + 4 * i - 4]))) / 3 / N

tmp <- PGEE::mycor_gee2(N, nt, y, X, gaussian(), beta_new, "AR-1", 1, 4, 1, 1, 1)
tmp <- PGEE::mycor_gee2(N, nt, y, X, gaussian(), beta_new, "exchangeable", 1, 4, 1, 1, 1)
tmp$Ehat[,,1]

identical(tmp$Ehat[,,1], tmp$Ehat[,,12])



sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
arma::uvec a4(int a) {
  return(randperm(a));
}
// [[Rcpp::export]]
Rcpp::List a5(arma::vec nt, int nCV) {
  Rcpp::List out(nCV);
  int N = nt.n_elem;
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec cvIndex = shuffle(index);
  arma::vec nTrain(nCV, arma::fill::value(round(N / nCV)));
  if (sum(nTrain) < N) for (int i = 0; i < N - sum(nTrain); i++) nTrain(i)++;
  if (sum(nTrain) > N) for (int i = 0; i < sum(nTrain) - N; i++) nTrain(i)--;
  double offset = 0;
  for (int i = 0; i < nCV; i++) {
    out(i) = sort(cvIndex(span(0 + offset, nTrain(i) - 1 + offset)));
    offset += nTrain(i);
  }
  return(out);
}')


a5(rep(4, 12), 5)
a5(rep(4, 14), 5)
a5(rep(4, 15), 5)
lapply(a5(rep(4, 12), 5), length)
lapply(a5(rep(4, 14), 5), length)
lapply(a5(rep(4, 15), 5), length)

sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
Rcpp::List a6(arma::vec idCV, arma::vec id) {
  Rcpp::List out(2);
  arma::uvec u;
  for (int i; i < idCV.n_elem; i++) {
    u = join_cols(u, find(id == idCV(i)));
  }
  arma::vec dummy(id.n_elem, arma::fill::zeros);
  std::cout << dummy;
  dummy.elem(u).ones();
  std::cout << dummy;
  out(0) = id(u);
  out(1) = id(find(dummy == 0));
  return(out);
}')


a6(2:4, rep(1:10, each = 2))
a7(2:4, rep(1:10, each = 2))

sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
arma::vec a6(arma::vec idCV, arma::vec id) {
  arma::uvec u;
  for (int i; i < idCV.n_elem; i++) {
    u = join_cols(u, find(id == idCV(i)));
  }
  arma::vec dummy(id.n_elem, arma::fill::zeros);
  dummy.elem(u).ones();
  // return(id(u));
  return(id(find(dummy == 0)));
}
// [[Rcpp::export]]
arma::vec a7(arma::vec idCV, arma::vec id) {
  arma::uvec u;
  for (int i; i < idCV.n_elem; i++) {
    u = join_cols(u, find(id == idCV(i)));
  }
  arma::vec dummy(id.n_elem, arma::fill::zeros);
  dummy.elem(u).ones();
  return(id(u));
}')


sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
// [[Rcpp::export]]
arma::mat a1(arma::mat a) {
  return(pinv(a, .01));
}
// [[Rcpp::export]]
arma::mat a2(arma::mat a) {
  arma::mat aa = eye(a.n_rows, a.n_cols);
  return(solve(a, aa, solve_opts::fast));
}')

a1(matrix(1:25, 5))
a2(matrix(1:25, 5))

library(microbenchmark)

microbenchmark(a1(matrix(1:25, 5)), invisible(a2(matrix(1:25, 5))))

