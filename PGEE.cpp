// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace arma;

arma::mat matvec(arma::mat x, arma::vec y) {
  arma::mat out(x.n_rows, x.n_cols);
  for (size_t i = 0; i < x.n_cols; i++) {
    out.col(i) = x.col(i) % y;
  }
  return out;
}

//' @noRd
// [[Rcpp::export(rng = false)]]
Rcpp::List SHM(arma::vec y,
	       arma::vec mu,
	       arma::mat X,
	       arma::cube Rhat, 
	       arma::vec nt,
	       arma::vec index){
  int N = nt.n_elem;
  int nx = X.n_cols;
  Rcpp::List out(3);
  arma::mat sumS(nx, 1, arma::fill::zeros);
  arma::mat sumH(nx, nx, arma::fill::zeros);
  arma::mat sumM(nx, nx, arma::fill::zeros);
  arma::vec ym = y - mu;
  // arma::mat bigD = matvec(X, muEta);
  for (int i = 0; i < N; i++) {
    arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
    arma::mat bigD = X.rows(span(index(i), index(i) + nt(i) - 1));
    arma::mat RhatSlice = Rhat.slice(i);
    arma::mat bigV = RhatSlice(span(0, nt(i) - 1), span(0, nt(i) - 1));
    arma::mat tmp = bigD.t() * pinv(bigV);
    sumS += tmp * ym2;
    sumH += tmp * bigD;
    tmp *= ym2;
    sumM += tmp * tmp.t();
  }
  out(0) = sumS;
  out(1) = sumH;
  out(2) = sumM;
  out.names() = Rcpp::CharacterVector::create("sumS", "sumH", "sumM");    
  return out;
}

// [[Rcpp::export(rng = false)]]
arma::vec qscad(arma::vec b, double lambda) {
  double a = 3.7;
  arma::vec dif = a * lambda - b;
  dif.elem(find(dif < 0)).zeros();
  arma::vec out(b.n_elem);
  out.fill(lambda);
  arma::uvec ind = find(b > lambda);
  out(ind) = dif(ind) / (a - 1);
  return(out);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List SHEM(arma::vec y,
		arma::mat X,
		arma::vec b0,
		arma::cube Rhat,
		arma::vec nt,
		arma::vec pindex, // 0 means to penalize
		double lambda,
		double eps){
  Rcpp::List out(4);
  int N = nt.n_elem;
  int nx = X.n_cols;
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec eta = X * b0;
  arma::vec E1 = qscad(abs(b0), lambda) / (abs(b0) + eps);
  E1.elem(find(pindex > 0)).zeros();
  arma::mat E = diagmat(E1);
  arma::mat S(nx, 1, arma::fill::zeros);
  arma::mat H(nx, nx, arma::fill::zeros);
  arma::mat M(nx, nx, arma::fill::zeros);
  arma::vec ym = y - eta;  
  for (int i = 0; i < N; i++) {
    arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
    arma::mat bigD = X.rows(span(index(i), index(i) + nt(i) - 1));
    arma::mat RhatSlice = Rhat.slice(i);
    arma::mat bigV = RhatSlice(span(0, nt(i) - 1), span(0, nt(i) - 1));
    arma::mat tmp = bigD.t() * pinv(bigV);
    S += tmp * ym2;
    H += tmp * bigD;
    tmp *= ym2;
    M += tmp * tmp.t();
  } 
  out(0) = S;
  out(1) = H;
  out(2) = E;
  out(3) = M;
  out.names() = Rcpp::CharacterVector::create("S", "H", "E", "M");    
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List gee(arma::vec y,
	       arma::mat X,
	       arma::vec b0,
	       arma::cube Rhat,
	       arma::vec nt,
	       arma::vec pindex, // 0 means to penalize
	       double lambda, double eps,
	       double tol, int maxit){
  Rcpp::List out(6);
  int N = nt.n_elem;
  int nx = X.n_cols;
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  for (int k = 0; k < maxit; k++) {
    arma::vec eta = X * b0;
    arma::vec E1 = qscad(abs(b0), lambda) / (abs(b0) + eps);
    E1.elem(find(pindex > 0)).zeros();
    arma::mat E = diagmat(E1);
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
    arma::vec ym = y - eta;  
    for (int i = 0; i < N; i++) {
      arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
      arma::mat bigD = X.rows(span(index(i), index(i) + nt(i) - 1));
      arma::mat RhatSlice = Rhat.slice(i);
      arma::mat bigV = RhatSlice(span(0, nt(i) - 1), span(0, nt(i) - 1));
      arma::mat tmp = bigD.t() * pinv(bigV);
      S += tmp * ym2;
      H += tmp * bigD;
      tmp *= ym2;
      M += tmp * tmp.t();
    }
    b1 = b0 + pinv(H + N * E) * (S - N * E * b0);
    out(0) = b1;
    out(1) = S;
    out(2) = H;
    out(3) = E;
    out(4) = M;
    out(5) = k;
    if(min(abs(b1 - b0)) < tol) break;
    b0 = b0;
  }
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter");
  return out;
}

//  
