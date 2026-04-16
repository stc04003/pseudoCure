#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// x, y, dx, and dy are assumed sorted

uvec group_ids_from_sorted(const vec& z) {
  const uword n = z.n_elem;
  uvec g(n);
  if (n == 0) return g;
  g[0] = 0;
  uword cur = 0;
  for (uword i = 1; i < n; ++i) {
    if (z[i] != z[i - 1]) ++cur;
    g[i] = cur;
  }
  return g;
}

vec unique_from_sorted(const vec& z) {
  if (z.n_elem == 0) return vec();
  std::vector<double> out;
  out.reserve(z.n_elem);
  out.push_back(z[0]);
  for (uword i = 1; i < z.n_elem; ++i) {
    if (z[i] != z[i - 1]) out.push_back(z[i]);
  }
  return vec(out);
}

uvec match_to_unique(const vec& z, const vec& uz) {
  const uword n = z.n_elem;
  uvec idx(n);
  const double* p = uz.memptr();
  const uword m = uz.n_elem;
  for (uword i = 0; i < n; ++i) {
    idx[i] = std::lower_bound(p, p + m, z[i]) - p;
  }
  return idx;
}

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::mat fastL(const arma::vec& x1,
                const arma::vec& y1,
                const arma::vec& dx1,
                const arma::vec& dy1) {
  const uword n = x1.n_elem;
  uvec ordX = sort_index(x1);
  uvec ordY = sort_index(y1);
  vec x  = x1(ordX);
  vec y  = y1(ordY);
  vec dx = dx1(ordX);
  vec dy = dy1(ordY);
  vec ux = unique_from_sorted(x);
  vec uy = unique_from_sorted(y);
  const uword nx = ux.n_elem;
  const uword ny = uy.n_elem;
  uvec gx_sorted = group_ids_from_sorted(x);
  uvec gy_sorted = group_ids_from_sorted(y);
  uvec gx_orig = match_to_unique(x1, ux);
  uvec gy_orig = match_to_unique(y1, uy);
  mat C(nx, ny, fill::zeros);
  mat Cdx(nx, ny, fill::zeros);
  mat Cdy(nx, ny, fill::zeros);
  mat C11(nx, ny, fill::zeros);
  for (uword k = 0; k < n; ++k) {
    uword i = gx_orig[k];
    uword j = gy_orig[k];
    C(i, j) += 1.0;
    if (dx1[k] > 0) Cdx(i, j) += 1.0;
    if (dy1[k] > 0) Cdy(i, j) += 1.0;
    if (dx1[k] > 0 && dy1[k] > 0) C11(i, j) += 1.0;
  }
  mat risk(nx, ny, fill::zeros);
  for (int i = (int)nx - 1; i >= 0; --i) {
    for (int j = (int)ny - 1; j >= 0; --j) {
      double val = C(i, j);
      if (i + 1 < (int)nx) val += risk(i + 1, j);
      if (j + 1 < (int)ny) val += risk(i, j + 1);
      if (i + 1 < (int)nx && j + 1 < (int)ny) val -= risk(i + 1, j + 1);
      risk(i, j) = val;
    }
  }
  mat row_dx(nx, ny, fill::zeros);
  for (uword i = 0; i < nx; ++i) {
    double run = 0.0;
    for (int j = (int)ny - 1; j >= 0; --j) {
      run += Cdx(i, j);
      row_dx(i, j) = run;
    }
  }
  mat col_dy(nx, ny, fill::zeros);
  for (uword j = 0; j < ny; ++j) {
    double run = 0.0;
    for (int i = (int)nx - 1; i >= 0; --i) {
      run += Cdy(i, j);
      col_dy(i, j) = run;
    }
  }
  mat Luniq(nx, ny, fill::zeros);
  for (uword i = 0; i < nx; ++i) {
    for (uword j = 0; j < ny; ++j) {
      double r = risk(i, j);
      if (r <= 0) continue;
      double l10 = row_dx(i, j) / r;
      double l01 = col_dy(i, j) / r;
      if (l10 == 1.0 || l01 == 1.0) continue;
      double l11 = C11(i, j) / r;
      Luniq(i, j) = (l10 * l01 - l11) / ((1.0 - l10) * (1.0 - l01));
    }
  }
  mat L(n, n, fill::zeros);
  for (uword u = 0; u < n; ++u) {
    if (dx[u] <= 0) continue;
    uword i = gx_sorted[u];
    for (uword v = 0; v < n; ++v) {
      if (dy[v] <= 0) continue;
      uword j = gy_sorted[v];
      L(u, v) = Luniq(i, j);
    }
  }
  return L;
}

arma::mat dabrowska_surface(const arma::vec& x,
                            const arma::vec& y,
                            const arma::vec& dx,
                            const arma::vec& dy,
                            const arma::mat& L) {
  const uword n = x.n_elem;
  vec Sx(n), Sy(n);
  for (uword i = 0; i < n; ++i) {
    double r = (double)(n - i);
    Sx[i] = (i == 0 ? 1.0 : Sx[i - 1]) * (1.0 - dx[i] / r);
    Sy[i] = (i == 0 ? 1.0 : Sy[i - 1]) * (1.0 - dy[i] / r);
  }
  // P(i,j) = prod_{a<=i, b<=j} (1 - L(a,b))
  // computed directly, without logs
  mat P(n, n, fill::ones);
  for (uword j = 0; j < n; ++j) {
    double col_prod = 1.0;
    for (uword i = 0; i < n; ++i) {
      col_prod *= (1.0 - L(i, j));
      if (j == 0) {
        P(i, j) = col_prod;
      } else {
        P(i, j) = P(i, j - 1) * col_prod;
      }
    }
  }
  mat D(n + 1, n + 1, fill::ones);
  for (uword i = 0; i < n; ++i) D(i + 1, 0) = Sx[i];
  for (uword j = 0; j < n; ++j) D(0, j + 1) = Sy[j];
  for (uword i = 0; i < n; ++i) {
    for (uword j = 0; j < n; ++j) {
      D(i + 1, j + 1) = Sx[i] * Sy[j] * P(i, j);
    }
  }
  return D;
}

//' @noRd
// [[Rcpp::export(rng = false)]]
double fastTau(const arma::vec& x1,
               const arma::vec& y1,
               const arma::vec& dx1,
               const arma::vec& dy1) {
  const uword n = x1.n_elem;
  uvec ordX = sort_index(x1);
  uvec ordY = sort_index(y1);
  vec x  = x1(ordX);
  vec y  = y1(ordY);
  vec dx = dx1(ordX);
  vec dy = dy1(ordY);
  mat L = fastL(x1, y1, dx1, dy1);
  mat D = dabrowska_surface(x, y, dx, dy, L);
  double out = 0.0;
  for (uword i = 1; i <= n; ++i) {
    for (uword j = 1; j <= n; ++j) {
      double d1 = D(i, j);
      double inc = D(i, j) - D(i, j - 1) - D(i - 1, j) + D(i - 1, j - 1);
      out += d1 * inc;
    }
  }
  return 4.0 * out - 1.0;
}

//' @noRd
// [[Rcpp::export(rng = false)]]
arma::vec pseudoTau(const arma::vec& x1,
                    const arma::vec& y1,
                    const arma::vec& dx1,
                    const arma::vec& dy1) {
  const uword n = x1.n_elem;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    vec x0(n - 1), y0(n - 1), dx0(n - 1), dy0(n - 1);
    if (i > 0) {
      x0.subvec(0, i - 1)  = x1.subvec(0, i - 1);
      y0.subvec(0, i - 1)  = y1.subvec(0, i - 1);
      dx0.subvec(0, i - 1) = dx1.subvec(0, i - 1);
      dy0.subvec(0, i - 1) = dy1.subvec(0, i - 1);
    }
    if (i + 1 < n) {
      x0.subvec(i, n - 2)  = x1.subvec(i + 1, n - 1);
      y0.subvec(i, n - 2)  = y1.subvec(i + 1, n - 1);
      dx0.subvec(i, n - 2) = dx1.subvec(i + 1, n - 1);
      dy0.subvec(i, n - 2) = dy1.subvec(i + 1, n - 1);
    }
    out[i] = fastTau(x0, y0, dx0, dy0);
  }
  return out;
}





















