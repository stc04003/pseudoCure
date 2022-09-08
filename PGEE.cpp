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
Rcpp::List SHM1(arma::vec y, arma::vec mu, arma::mat X,
								arma::cube Rhat, int N, int nx, arma::vec nt,
								arma::vec index, arma::vec muEta){
	Rcpp::List out(3);
	arma::mat sumS = zeros<mat>(nx, 1);
	arma::mat sumH = zeros<mat>(nx, nx);
	arma::mat sumM = zeros<mat>(nx, nx);
	for (int i = 0; i < N; i++) {
		arma::vec ym = zeros<vec>(nt(i));
		arma::mat bigD = zeros<mat>(nt(i), nx);
		for (int j = 0; j < nt(i); j++) {
			ym(j) = y(j + index(i)) -mu(j + index(i));
			for (int k = 0; k < nx; k++) {
	      bigD(j, k) = muEta(j + index(i)) * X(j + index(i), k);
			}
		}
		arma::mat RhatSlice = Rhat.slice(i);
		arma::mat bigV = RhatSlice(span(0, nt(i)-1), span(0, nt(i)-1));
		arma::mat tempS = bigD.t() * pinv(bigV) * ym;
		sumS = sumS + tempS;
		arma::mat SRhat = pinv(bigV);
		arma::mat tempH = bigD.t() * SRhat * bigD;
		sumH = sumH + tempH;
		arma::mat tempM = bigD.t() * SRhat * ym * ym.t() * SRhat * bigD;
		sumM = sumM +tempM;
	}
	out(0) = sumS;
	out(1) = sumH;
	out(2) = sumM;
	out.names() = Rcpp::CharacterVector::create("sumS", "sumH", "sumM");    
	return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List SHM2(arma::vec y,
								arma::vec mu, arma::mat X,
								arma::cube Rhat, int N, int nx, arma::vec nt,
								arma::vec index, arma::vec muEta){
	Rcpp::List out(3);      
	arma::mat sumS(nx, 1, arma::fill::zeros);
	arma::mat sumH(nx, nx, arma::fill::zeros);
	arma::mat sumM(nx, nx, arma::fill::zeros);
	arma::vec ym = y - mu;
	arma::mat bigD = matvec(X, muEta);
	for (int i = 0; i < N; i++) {
		arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
		arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
		arma::mat RhatSlice = Rhat.slice(i);
		arma::mat bigV = RhatSlice(span(0, nt(i) - 1), span(0, nt(i) - 1));
		sumS += bigD2.t() * pinv(bigV) * ym2;
		arma::mat SRhat = pinv(bigV);
		sumH += bigD2.t() * SRhat * bigD2;
		sumM += bigD2.t() * SRhat * ym2 * ym2.t() * SRhat * bigD2;
	}
	out(0) = sumS;
	out(1) = sumH;
	out(2) = sumM;
	out.names() = Rcpp::CharacterVector::create("sumS", "sumH", "sumM");    
	return out;
}
