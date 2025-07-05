#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// x, y, dx, and dy are assumed sorted
//' @noRd
// [[Rcpp::export]]
double fastDabrowska(double tx, double ty,
		     const arma::vec& x, const arma::vec& y,
		     const arma::vec& dx, const arma::vec& dy,
		     const arma::mat& L) {
  int n = x.n_elem;
  vec r = linspace<vec>(n, 1, n);
  vec Sx = cumprod(1 - dx / r);
  vec Sy = cumprod(1 - dy / r);
  // Fast index lookup
  uword i = 0, j = 0;
  while (i < n && x[i] <= tx) ++i;
  while (j < n && y[j] <= ty) ++j;
  if (i == 0 && j == 0) return 1.0;
  if (i == 0) return Sy[j - 1];
  if (j == 0) return Sx[i - 1];
  // Compute product of (1 - L(u,v)) over u = 0..i-1, v = 0..j-1
  double prodL = 1.0;
  for (uword u = 0; u < i; ++u) {
    for (uword v = 0; v < j; ++v) {
      prodL *= (1.0 - L(u, v));
    }
  }
  return Sx[i - 1] * Sy[j - 1] * prodL;
}

//' @noRd
// [[Rcpp::export]]
arma::mat fastL(const arma::vec& x1, const arma::vec& y1,
		const arma::vec& dx1, const arma::vec& dy1) {
  int n = x1.n_elem;
  uvec ordX = sort_index(x1);
  uvec ordY = sort_index(y1);
  vec x = x1.elem(ordX);
  vec y = y1.elem(ordY);
  vec dx = dx1.elem(ordX);
  vec dy = dy1.elem(ordY);
  mat L(n, n, fill::zeros);
  for (int u = 0; u < n; ++u) {
    if (dx[u] > 0) {
      double xcut = x[u];
      for (int v = 0; v < n; ++v) {
        if (dy[v] > 0) {
          double ycut = y[v];
          int risk = 0, l10 = 0, l01 = 0, l11 = 0;
          for (int i = 0; i < n; ++i) {
            if (x1[i] >= xcut && y1[i] >= ycut) ++risk;
            if (x1[i] == xcut && y1[i] >= ycut && dx1[i] == 1) ++l10;
            if (x1[i] >= xcut && y1[i] == ycut && dy1[i] == 1) ++l01;
            if (x1[i] == xcut && y1[i] == ycut && dx1[i] == 1 && dy1[i] == 1) ++l11;
          }
          if (risk > 0) {
            double d10 = static_cast<double>(l10) / risk;
            double d01 = static_cast<double>(l01) / risk;
            if (d10 != 1.0 && d01 != 1.0) {
              double d11 = static_cast<double>(l11) / risk;
              L(u, v) = (d10 * d01 - d11) / ((1.0 - d10) * (1.0 - d01));
            }
          }
        }
      }
    }
  }
  return L;
}

// [[Rcpp::export]]
arma::vec fastDabrowska_vec(const arma::vec& tx,
                            const arma::vec& ty,
                            const arma::vec& x,
                            const arma::vec& y,
                            const arma::vec& dx,
                            const arma::vec& dy,
                            const arma::mat& L) {
  int n = x.n_elem;
  int K = tx.n_elem;
  if (ty.n_elem != K)
    stop("tx and ty must have the same length.");
  // Precompute survival curves
  vec r = linspace<vec>(n, 1, n);
  vec Sx = cumprod(1 - dx / r);
  vec Sy = cumprod(1 - dy / r);
  vec out(K, fill::zeros);
  for (int k = 0; k < K; ++k) {
    double t_x = tx[k];
    double t_y = ty[k];
    // Find last i where x[i] <= tx[k]
    uword i = 0;
    while (i < n && x[i] <= t_x) ++i;
    if (i > 0) --i;
    // Find last j where y[j] <= ty[k]
    uword j = 0;
    while (j < n && y[j] <= t_y) ++j;
    if (j > 0) --j;
    // Handle edge cases
    if (t_x == 0 && t_y == 0) {
      out[k] = 1.0;
    } else if (t_x == 0 && t_y != 0) {
      out[k] = Sy[j];
    } else if (t_x != 0 && t_y == 0) {
      out[k] = Sx[i];
    } else {
      // Compute product over L(0:i, 0:j)
      double prodL = 1.0;
      for (uword u = 0; u <= i; ++u) {
        for (uword v = 0; v <= j; ++v) {
          prodL *= (1.0 - L(u, v));
        }
      }
      out[k] = Sx[i] * Sy[j] * prodL;
    }
  }
  return out;
}


//' @noRd
// [[Rcpp::export]]
double fastTau(const arma::vec& x1,
	       const arma::vec& y1,
	       const arma::vec& dx1,
	       const arma::vec& dy1) {
  int n = x1.n_elem;
  // Sort inputs
  uvec ordX = sort_index(x1);
  uvec ordY = sort_index(y1);
  vec x = x1.elem(ordX);
  vec y = y1.elem(ordY);
  vec dx = dx1.elem(ordX);
  vec dy = dy1.elem(ordY);
  // Compute L once
  arma::mat L = fastL(x1, y1, dx1, dy1);
  // Construct (x0, y0) = (0, x[0], ..., x[n-1]), same for y
  vec x0(n + 1), y0(n + 1);
  x0[0] = 0.0;
  y0[0] = 0.0;
  for (int i = 0; i < n; ++i) {
    x0[i + 1] = x[i];
    y0[i + 1] = y[i];
  }
  // Prepare all pairs (tx, ty) for Dabrowska evaluation
  vec tx((n + 1) * (n + 1));
  vec ty((n + 1) * (n + 1));
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= n; ++j) {
      int idx = i * (n + 1) + j;
      tx[idx] = x0[i];
      ty[idx] = y0[j];
    }
  }
  // Vectorized evaluation
  vec dab_vec = fastDabrowska_vec(tx, ty, x, y, dx, dy, L);
  mat dab(dab_vec.memptr(), n + 1, n + 1, false);
  // Vectorized computation of sum
  mat D = dab(span(1, n), span(1, n));
  mat D_left = dab(span(1, n), span(0, n - 1));
  mat D_up = dab(span(0, n - 1), span(1, n));
  mat D_corner = dab(span(0, n - 1), span(0, n - 1));
  mat contrib = D % (D - D_left - D_up + D_corner);
  return 4.0 * accu(contrib) - 1.0;
}


//' @noRd
// [[Rcpp::export]]
arma::vec pseudoTau(arma::vec x1, arma::vec y1,
		    arma::vec dx1, arma::vec dy1) {
  int n = x1.size();
  vec out = zeros<vec>(n);
  vec ind = linspace(0, n - 1, n);
  for (int i = 0; i < n; i++) {
    uvec ind2 = find(i != ind);
    vec x10 = x1(ind2);
    vec dx10 = dx1(ind2);
    vec y10 = y1(ind2);
    vec dy10 = dy1(ind2);
    out[i] = fastTau(x10, y10, dx10, dy10);
  }
  return(out);
}
