// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <set>
#include <utility>

using namespace arma;

// [[Rcpp::export(rng = false)]]
double minS(const arma::vec Time,
	    const arma::vec censor) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1));
    r(span(0, i)) += ind1.n_elem;
  }
  return(prod(1 - d / r));
}

// The first N elements are the leave i out min KM
// The last element is the min KM using all data
// [[Rcpp::export(rng = false)]]
arma::vec minSi(const arma::vec Time,
		const arma::vec censor) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  int N = Time.n_elem;
  arma::vec out0(n, arma::fill::ones); // ordered
  arma::vec out1(n, arma::fill::ones); // ordered
  arma::vec out2(N + 1, arma::fill::ones);
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  // calculate d and r
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1));
    r(span(0, i)) += ind1.n_elem;
  }
  // Calculate leave-one-out KM
  out0(span(1, n - 1)) *= 1 - d[0] / (r[0] - 1);
  out1(span(1, n - 1)) *= 1 - d[0] / (r[0] - 1);
  if (r[0] > 1) {
    out1(0) = 1 - (d[0] - 1) / (r[0] - 1);
    out0(0) = 1 - d[0] / (r[0] - 1);    
  }    
  for (int i = 1; i < n - 1; i++) {
    out1(span(0, i - 1)) *= 1 - d[i] / r[i];
    out1(span(i + 1, n - 1)) *= 1 - d[i] / (r[i] - 1);
    out0(span(0, i - 1)) *= 1 - d[i] / r[i];
    out0(span(i + 1, n - 1)) *= 1 - d[i] / (r[i] - 1);
    if (r[i] > 1) {
      out1(i) *= 1 - (d[i] - 1) / (r[i] - 1);
      out0(i) *= 1 - d[i] / (r[i] - 1);
    }
  }
  int i = n - 1;
  out1(span(0, i - 1)) *= 1 - d[i] / r[i];
  out0(span(0, i - 1)) *= 1 - d[i] / r[i];
  if (r[i] > 1) {
    out1(i) *= 1 - (d[i] - 1) / (r[i] - 1);
    out0(i) *= 1 - d[i] / (r[i] - 1);
  }
  // Reassign
  for (int i = 0; i < N; i++) {
    arma::vec tmp;
    if (censor[i] == 0) tmp = out0.elem(find(Time[i] == T0));
    else tmp = out1.elem(find(Time[i] == T0));
    out2[i] = tmp[0];
  }
  out2[N] = prod(1 - d / r); 
  return(out2);
}
