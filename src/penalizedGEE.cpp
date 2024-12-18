#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

arma::mat matvec(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_cols; i++) {
    out.col(i) = x.col(i) % y;
  }
  return out;
}

// Functions to calculate parameter in variance-covariance matrix
// Those functions assume equal cluster size
double ahatEx(arma::vec a,
	      arma::vec nt,
	      arma::vec index) {
  int n = nt.n_elem;
  double out = 0;
  for (int i = 0; i < n; i++) {
    arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
    out += sum(a2) * sum(a2) - sum(a2 % a2);
  }
  return(out / nt(0) / (nt(0) - 1) / n / mean(a % a));
}

double ahatAR1(arma::vec a,
	       arma::vec nt,
	       arma::vec index) {
  int n = nt.n_elem;
  double out = 0;
  for (int i = 0; i < n; i++) {
    arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
    out += sum(a2(span(0, nt(0) - 2)) % a2(span(1, nt(0) - 1)));
  }
  return(out / (nt(0) - 1) / n / mean(a % a));
}

// link functions
// Those functions return a list contains elements 1. linkinv 2. mu.eta
Rcpp::List identityLink (arma::vec eta) {
  Rcpp::List out(2);
  out(0) = eta;
  out(1) = vec(eta.n_elem, arma::fill::ones);
  return out;
}

Rcpp::List logLink (arma::vec eta) {
  Rcpp::List out(2);
  arma::vec eeta = exp(eta);
  eeta.replace(arma::datum::inf, pow(2, 1023));
  out(0) = eeta;
  out(1) = eeta;
  return out;
}

Rcpp::List cloglogLink (arma::vec eta) {
  Rcpp::List out(2);
  arma::vec tmp = eta;	
  tmp.elem(find(tmp > 700)).fill(700);
  tmp = exp(tmp) % exp(-exp(tmp));
  tmp.replace(arma::datum::inf, pow(2, 1023));
  out(0) = 1 - exp(-exp(eta));
  out(1) = tmp;
  return out;
}

Rcpp::List logitLink (arma::vec eta) {
  Rcpp::List out(2);
  arma::vec eeta = exp(eta);
  eeta.replace(arma::datum::inf, pow(2, 1023));
  out(0) = eeta / (eeta + 1);
  out(1) = eeta / (eeta + 1) / (eeta + 1);
  return out;
}


// Penalty functions
arma::vec qscad(arma::vec b, double lambda) {
  double a = 3.7;
  arma::vec dif = a * lambda - b;
  dif.elem(find(dif < 0)).zeros();
  arma::vec out(b.n_elem, arma::fill::value(lambda));
  // out.fill(lambda);
  arma::uvec ind = find(b > lambda);
  out(ind) = dif(ind) / (a - 1);
  return(out);
}

arma::vec qlasso(arma::vec b, double lambda) {
  arma::vec out(b.n_elem, arma::fill::value(lambda));
  out.elem(find(b > lambda)).zeros();
  return(out);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List pgee(arma::vec y,
		arma::mat X,
		arma::vec b0,
		arma::vec nt,
		arma::vec pindex, // 0 means to penalize
		std::string glmlink,
		std::string penalty,
		std::string corstr,
		double lambda,
		double eps,
		double tol,
		int maxit){
  Rcpp::List out(7);
  int N = nt.n_elem;
  int nx = X.n_cols;
  int k = nt(0);
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  for (int j = 0; j <= maxit; j++) {
    arma::vec eta = X * b0;
    arma::vec E1(nx, arma::fill::zeros);
    if (penalty == "scad")
      E1 = qscad(abs(b0), lambda) / (abs(b0) + eps);
    if (penalty == "lasso")
      E1 = qlasso(abs(b0), lambda) / (abs(b0) + eps);
    E1.elem(find(pindex > 0)).zeros();
    arma::mat E = diagmat(E1);
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    Rcpp::List links;
    if (glmlink == "identity") links = identityLink(eta);
    if (glmlink == "log") links = logLink(eta);
    if (glmlink == "cloglog") links = cloglogLink(eta);
    if (glmlink == "logit") links = logitLink(eta);
    arma::vec mu = links(0);
    arma::vec etamu = links(1);
    arma::vec ym = y - mu;
    arma::mat bigD = matvec(X, etamu);
    arma::mat Rhat(k, k, arma::fill::eye);
    double ahat = 0;         
    if (corstr == "ex" && k > 1) {
      ahat = ahatEx(ym / stddev(mu), nt, index);
      Rhat = Rhat * (1 - ahat) + ahat;
    }
    if (corstr == "ar1" && k > 1) {
      ahat = ahatAR1(ym / stddev(mu), nt, index);
      arma::mat Rhat2(k, k, arma::fill::zeros);
      arma::vec tmp(k - 1, arma::fill::value(ahat));
      tmp = cumprod(tmp);
      for (int i = 0; i < k - 1; i++) {
	Rhat2.submat(i + 1, i, k - 1, i) = tmp(span(0, k - i - 2));
      }
      Rhat = Rhat + Rhat2 + Rhat2.t();
    }
    // std::cout << "ahat: " << ahat << "\n";
    for (int i = 0; i < N; i++) {
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat tmp = bigD2.t() * pinv(Rhat);
      S += tmp * ym2;
      H += tmp * bigD2;
      tmp *= ym2;
      M += tmp * tmp.t();
    }
    b1 = b0 + pinv(H + N * E) * (S - N * E * b0);
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(3) = E;
    out(4) = M;
    out(5) = j;
    out(6) = ahat;
    if(min(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter", "alpha");
  return out;
}

// no penality version
// [[Rcpp::export(rng = false)]]
Rcpp::List gee(arma::vec y,
	       arma::mat X,
	       arma::vec b0,
	       arma::vec nt,
	       std::string glmlink,
	       std::string corstr,
	       double tol,
	       int maxit){
  Rcpp::List out(7);
  int N = nt.n_elem;
  int nx = X.n_cols;
  int k = nt(0);
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  for (int j = 0; j <= maxit; j++) {
    arma::vec eta = X * b0;
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    Rcpp::List links;
    if (glmlink == "identity") links = identityLink(eta);
    if (glmlink == "log") links = logLink(eta);
    if (glmlink == "cloglog") links = cloglogLink(eta);
    if (glmlink == "logit") links = logitLink(eta);
    arma::vec mu = links(0);
    arma::vec etamu = links(1);
    arma::vec ym = y - mu;
    arma::mat bigD = matvec(X, etamu);
    arma::mat Rhat(k, k, arma::fill::eye);
    double ahat = 0;
    if (corstr == "ex" && k > 1) {
      if (stddev(mu) > 0)
	ahat = ahatEx(ym / stddev(mu), nt, index);
      else
	ahat = ahatEx(ym, nt, index);
      Rhat = Rhat * (1 - ahat) + ahat;
    }
    if (corstr == "ar1" && k > 1) {
      if (stddev(mu) > 0)
	ahat = ahatAR1(ym / stddev(mu), nt, index);
      else
	ahat = ahatAR1(ym, nt, index);
      arma::vec tmp(k - 1, arma::fill::value(ahat));
      arma::mat Rhat2(k, k, arma::fill::zeros);
      tmp = cumprod(tmp);
      for (int i = 0; i < k - 1; i++) {
	Rhat2.submat(i + 1, i, k - 1, i) = tmp(span(0, k - i - 2));
      }
      Rhat = Rhat + Rhat2 + Rhat2.t();
    }
    for (int i = 0; i < N; i++) {
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat tmp = bigD2.t() * pinv(Rhat);
      S += tmp * ym2;
      H += tmp * bigD2;
      tmp *= ym2;
      M += tmp * tmp.t();
    }
    b1 = b0 + pinv(H) * S;
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(4) = M;
    out(5) = j;
    out(6) = ahat;
    if(min(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  arma::mat E(nx, nx, arma::fill::zeros);
  out(3) = E;
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter", "alpha");
  return out;
}


// Functions need for cross-validation
Rcpp::List getCV(arma::vec index, int nCV) {
  Rcpp::List out(nCV);
  int N = index.n_elem;
  arma::vec cvIndex = shuffle(index);
  arma::vec nTrain(nCV, arma::fill::value(round(N / nCV)));
  if (sum(nTrain) < N) 
    nTrain(span(0, N - sum(nTrain) - 1)) += 1;
  if (sum(nTrain) > N) 
    nTrain(span(0, sum(nTrain) - N - 1)) -= 1;
  double offset = 0;
  for (int i = 0; i < nCV; i++) {
    out(i) = sort(cvIndex(span(0 + offset, nTrain(i) - 1 + offset)));
    offset += nTrain(i);
  }
  return(out);
}

arma::vec getID(arma::vec nt) {
  arma::vec out(sum(nt));
  int offset = 0;
  int n = nt.n_elem;
  for (int i = 0; i < n; i++) {
    out(span(offset, nt(i) - 1 + offset)).fill(i);
    offset += nt(i);
  } 
  return(out);
}

// [[Rcpp::export]]
arma::mat pgeeCV(arma::vec y,
		 arma::mat X,
		 arma::vec b0,
		 arma::vec nt,
		 arma::vec pindex, // 0 means to penalize
		 std::string glmlink,
		 std::string penalty,
		 std::string corstr,
		 int nCV, 
		 arma::vec lambda,
		 double eps,
		 double tol,
		 int maxit){
  arma::mat out(nCV, lambda.n_elem, arma::fill::zeros);
  arma::vec cvm(nCV, arma::fill::zeros);
  arma::vec id = getID(nt);
  Rcpp::List cvList = getCV(unique(id), nCV);
  for (int i = 0; i < nCV; i++) {
    arma::uvec idTest;
    arma::vec idTrain(id.n_elem, arma::fill::zeros);
    arma::vec idCV = cvList(i);
    for (int j = 0; j < (int) idCV.n_elem; j++) {
      idTest = arma::join_cols(idTest, find(id == idCV(j)));
    }
    idTrain.elem(idTest).ones();
    arma::vec yTest = y(idTest);
    arma::mat XTest = X.rows(idTest);
    arma::vec yTrain = y(find(idTrain == 0));
    arma::mat XTrain = X.rows(find(idTrain == 0));
    // assume equal cluster size
    arma::vec ntTrain = nt(span(1, idCV.n_elem));
    for (int j = 0; j < (int) lambda.n_elem; j++) {			
      Rcpp::List tmp = pgee(yTrain, XTrain, b0, ntTrain, pindex, glmlink, penalty, "ind",
			    lambda(j), eps, tol, maxit);
      arma::vec b1 = tmp(0);
      arma::vec eta = XTest * b1;
      Rcpp::List links;
      if (glmlink == "identity") links = identityLink(eta);
      if (glmlink == "log") links = logLink(eta);
      if (glmlink == "cloglog") links = cloglogLink(eta);
      if (glmlink == "logit") links = logitLink(eta);
      arma::vec mu = links(0);
      arma::vec devResids = yTest - mu;    
      out(i, j) = mean(devResids % devResids); // transform to original scale? 
    }
  }
  // out(0) = mean(cvm);
  // out(1) = stddev(cvm / nCV);
  return(out);
}
