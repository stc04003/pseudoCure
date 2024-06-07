#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export]]
double fastDabrowska(double tx, double ty,
		     arma::vec x, arma::vec y, arma::vec dx, arma::vec dy, arma::mat L) {
  double out = 1;
  int n = x.size();
  vec r = linspace(n, 1, n);
  vec Sx = cumprod(1 - dx / r);
  vec Sy = cumprod(1 - dy / r);
  uvec rX = find(x <= tx);
  uvec rY = find(y <= ty);
  if (tx != 0 && ty != 0) {
    vec Sxtx = Sx((uvec) rX.tail(1));
    vec Syty = Sy((uvec) rY.tail(1));
    out = Sxtx[0] * Syty[0] * prod(prod(1 - L.submat(0, 0, rX.size() - 1, rY.size() - 1)));
  }
  if (tx == 0 && ty != 0) {
    vec Syty = Sy((uvec) rY.tail(1));
    out = Syty[0];
  }
  if (tx != 0 && ty == 0) {
    vec Sxtx = Sx((uvec) rX.tail(1));
    out = Sxtx[0];
  }
  return(out);
}

//' @noRd
// [[Rcpp::export]]
arma::mat fastL(arma::vec x1, arma::vec y1, arma::vec dx1, arma::vec dy1) {
  int n = x1.size();
  uvec ordX = sort_index(x1);
  uvec ordY = sort_index(y1);
  vec x = x1(ordX);
  vec y = y1(ordY);
  vec dx = dx1(ordX);
  vec dy = dy1(ordY);
  mat L(n, n, fill::zeros);
  for (int u = 0; u < n; u++) {
    if (dx[u] > 0) {
      for (int v = 0; v < n; v++) {
        if (dy[v] > 0) {
          double risk = ((uvec) find(x1 >= x[u] && y1 >= y[v])).size();
          if (risk != 0) {
            double l10 = ((uvec) find(x1 == x[u] && y1 >= y[v] && dx1 > 0)).size() / risk;
            double l01 = ((uvec) find(x1 >= x[u] && y1 == y[v] && dy1 > 0)).size() / risk;
            if (l10 != 1 && l01 != 1) {
              double l11 = ((uvec) find(x1 == x[u] && y1 == y[v] && dx1 > 0 && dy1 > 0)).size() / risk;
              L(u, v) = (l10 * l01 - l11) / (1 - l10) / (1 - l01);
            }
          }
        }
      }
    } 
  } 
  return(L);
}

//' @noRd
// [[Rcpp::export]]
double fastTau2(arma::vec x1, arma::vec y1, arma::vec dx1, arma::vec dy1) {
  int n = x1.size();
  uvec ordX = sort_index(x1);
  uvec ordY = sort_index(y1);
  vec x = x1(ordX);
  vec y = y1(ordY);
  vec dx = dx1(ordX);
  vec dy = dy1(ordY);
  mat L(n, n, fill::zeros);
  L = fastL(x1, y1, dx1, dy1);
  double out = 0;
  mat dab(n + 1, n + 1, fill::zeros);
  vec x0 = linspace(0, n, n + 1);
  vec y0 = linspace(0, n, n + 1);
  for (int i = 0; i < n; i++) {
    x0[i + 1] = x[i];
    y0[i + 1] = y[i];
  }
  for (int i = 0; i < (n + 1); i++) {
    for (int j = 0; j < (n + 1); j++) {
      dab(i, j) = fastDabrowska(x0[i], y0[j], x, y, dx, dy, L);
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (dab(i + 1, j + 1) > 0)
	out = out + dab(i + 1, j + 1) * (dab(i + 1, j + 1) - dab(i + 1, j) - dab(i, j + 1) + dab(i, j));
      // double d1 = dab[i + 1, j + 1];
      // double d2 = dab[i + 1, j];
      // double d3 = dab[i, j + 1];
      // double d4 = dab[i, j];
      // out = out + d1 * (d1 - d2 - d3 + d4);
    }
  }
  return(4 * out - 1);
}

//' @noRd
// [[Rcpp::export]]
arma::vec pseudoTau(arma::vec x1, arma::vec y1, arma::vec dx1, arma::vec dy1) {
  int n = x1.size();
  vec out = zeros<vec>(n);
  vec ind = linspace(0, n - 1, n);
  for (int i = 0; i < n; i++) {
    uvec ind2 = find(i != ind);
    vec x10 = x1(ind2);
    vec dx10 = dx1(ind2);
    vec y10 = y1(ind2);
    vec dy10 = dy1(ind2);
    out[i] = fastTau2(x10, y10, dx10, dy10);
  }
  return(out);
}
