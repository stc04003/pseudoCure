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
  out1(span(1, n - 1)) *= 1 - d[0] / (r[0] - 1);
  if (r[0] > 1) out1(0) = 1 - (d[0] - 1) / (r[0] - 1);
  for (int i = 1; i < n - 1; i++) {
    out1(span(0, i - 1)) *= 1 - d[i] / r[i];
    out1(span(i + 1, n - 1)) *= 1 - d[i] / (r[i] - 1);
    if (r[i] > 1) out1(i) *= 1 - (d[i] - 1) / (r[i] - 1);
  }
  int i = n - 1;
  out1(span(0, i - 1)) *= 1 - d[i] / r[i];
  if (r[i] > 1) out1(i) *= 1 - (d[i] - 1) / (r[i] - 1);
  arma::vec out0 = out1;
  for (int i = 0; i < n - 1; i++) out0(i) *= 1 - 1 / (r[i] - d[i]);
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


// Leave one out at different quantile
// slower than minSi but more useful when calculating multiple S(t) at t0
// [[Rcpp::export(rng = false)]]
arma::vec minSi1(const arma::vec Time,
								 const arma::vec censor) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  int N = Time.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  arma::vec out(N + 1, arma::fill::ones);
  // calculate d and r for all data
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1));
    r(span(0, i)) += ind1.n_elem;
  }
  // calculate d and r for each leave one out object
  arma::mat dmat = repmat(d, 1, N);
  arma::mat rmat = repmat(r, 1, N);
  for (int i = 0; i < N; i++) {
    arma::uvec ind = find(Time[i] == T0, 1);
		// std::cout << "ind: " << ind << "\n";
		int tmp = (int) ind(0);
		if (censor[i] > 0) dmat(tmp, i) -= 1;
		rmat(span(0, tmp), i) -= 1; 
		arma::vec dr = dmat.col(i) / rmat.col(i);
		dr.replace(arma::datum::nan, 0);
    dr.replace(arma::datum::inf, 0);
    out[i] = prod(1 - dr);
  }
  out[N] = prod(1 - d / r);
  return(out);
}

// [[Rcpp::export(rng = false)]]
arma::vec minSi2(const arma::vec Time,
								 const arma::vec censor) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  int N = Time.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  arma::vec out(N + 1, arma::fill::ones);
  // calculate d and r for all data
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1));
    r(span(0, i)) += ind1.n_elem;
  }	
	for (int i = 0; i < N; i++) {
		arma::vec d2 = d;
		arma::vec r2 = r;
    arma::uvec ind = find(Time[i] == T0, 1);
		int tmp = (int) ind(0);
		if (censor[i] > 0) d2(tmp) -= 1;
		r2(span(0, tmp)) -= 1; 
		arma::vec dr = d2 / r2;
		dr.replace(arma::datum::nan, 0);
    dr.replace(arma::datum::inf, 0);
    out[i] = prod(1 - dr);
  }
  out[N] = prod(1 - d / r);
  return(out);
}


// /////////

// [[Rcpp::export(rng = false)]]
arma::mat minSi3(const arma::vec Time,
								 const arma::vec censor,
								 const arma::vec Q) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  int N = Time.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros); 
  arma::mat out(Q.n_elem, N + 1, arma::fill::ones);
  // calculate d and r for all data
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1));
    r(span(0, i)) += ind1.n_elem;
  }
	// Identify quantiles; like findInterval()
	arma::uvec ind1(Q.n_elem, arma::fill::ones);;
	for (int i = 0; i < Q.n_elem; i++) {
		arma::uvec tmp = find(Q(i) >= T0, 1, "last");
		ind1(i) = tmp(0);
	}
  // calculate d and r for each leave one out object
  arma::mat dmat = repmat(d, 1, N);
  arma::mat rmat = repmat(r, 1, N);
  for (int i = 0; i < N; i++) {
    arma::uvec ind2 = find(Time[i] == T0, 1);
		int tmp = (int) ind2(0);
		if (censor[i] > 0) dmat(tmp, i) -= 1;
		rmat(span(0, tmp), i) -= 1; 
		arma::vec dr = dmat.col(i) / rmat.col(i);
		dr.replace(arma::datum::nan, 0);
    dr.replace(arma::datum::inf, 0);
		arma::vec S = cumprod(1 - dr);
		for (int j = 0; j < Q.n_elem; j++) {
			out(j, i) = S(ind1(j));
		}
  }
	arma::vec S = cumprod(1 - d / r); 
  out.col(N) = S(ind1);
  return(out);
}

// [[Rcpp::export(rng = false)]]
arma::mat minSi4(const arma::vec Time,
								 const arma::vec censor,
								 const arma::vec Q) {
  arma::vec T0 = arma::sort(arma::unique(Time));
  int n = T0.n_elem;
  int N = Time.n_elem;
  arma::vec d(n, arma::fill::zeros);
  arma::vec r(n, arma::fill::zeros);
	arma::mat out(Q.n_elem, N + 1, arma::fill::ones);
  // calculate d and r for all data
  for (int i = 0; i < n; i++) {
    arma::uvec ind1 = find(Time == T0[i]);
    d[i] = sum(censor.elem(ind1));
    r(span(0, i)) += ind1.n_elem;
  }
	// Identify quantiles; like findInterval()
	arma::uvec ind1(Q.n_elem, arma::fill::ones);;
	for (int i = 0; i < Q.n_elem; i++) {
		arma::uvec tmp = find(Q(i) >= T0, 1, "last");
		ind1(i) = tmp(0);
	}
	for (int i = 0; i < N; i++) {
		arma::vec d2 = d;
		arma::vec r2 = r;
    arma::uvec ind = find(Time[i] == T0, 1);
		int tmp = (int) ind(0);
		if (censor[i] > 0) d2(tmp) -= 1;
		r2(span(0, tmp)) -= 1; 
		arma::vec dr = d2 / r2;
		dr.replace(arma::datum::nan, 0);
    dr.replace(arma::datum::inf, 0);
		arma::vec S = cumprod(1 - dr);
		for (int j = 0; j < Q.n_elem; j++) {
			out(j, i) = S(ind1(j));
		}
  }
	arma::vec S = cumprod(1 - d / r); 
  out.col(N) = S(ind1);
  return(out);
}
