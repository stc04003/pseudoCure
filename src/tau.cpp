#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Helper: unique values from a sorted vector.
arma::vec unique_from_sorted(const arma::vec& z) {
  if (z.n_elem == 0) return arma::vec();

  std::vector<double> out;
  out.reserve(z.n_elem);
  out.push_back(z[0]);

  for (arma::uword i = 1; i < z.n_elem; ++i) {
    if (z[i] != z[i - 1]) out.push_back(z[i]);
  }

  return arma::vec(out);
}

// Helper: map each z value to its index in sorted unique vector uz.
arma::uvec match_to_unique(const arma::vec& z, const arma::vec& uz) {
  const arma::uword n = z.n_elem;
  arma::uvec idx(n);

  const double* p = uz.memptr();
  const arma::uword m = uz.n_elem;

  for (arma::uword i = 0; i < n; ++i) {
    idx[i] = std::lower_bound(p, p + m, z[i]) - p;
  }

  return idx;
}

// Helper: Kaplan-Meier survival evaluated once per unique sorted time.
// z_sorted and dz_sorted must be sorted by z_sorted.
arma::vec marginal_surv_unique(const arma::vec& z_sorted,
                               const arma::vec& dz_sorted) {
  const arma::uword n = z_sorted.n_elem;
  arma::vec uz = unique_from_sorted(z_sorted);
  const arma::uword m = uz.n_elem;

  arma::vec S(m, arma::fill::ones);
  arma::uword pos = 0;
  double surv = 1.0;

  for (arma::uword k = 0; k < m; ++k) {
    const double t = uz[k];
    const arma::uword start = pos;
    double d = 0.0;

    while (pos < n && z_sorted[pos] == t) {
      if (dz_sorted[pos] > 0) d += 1.0;
      ++pos;
    }

    const double risk = static_cast<double>(n - start);
    if (risk > 0.0) {
      surv *= (1.0 - d / risk);
    }

    S[k] = surv;
  }

  return S;
}

// Internal: Dabrowska L jump matrix on the UNIQUE x-y grid.
arma::mat fastL_unique_internal(const arma::vec& x1,
                                const arma::vec& y1,
                                const arma::vec& dx1,
                                const arma::vec& dy1) {
  const arma::uword n = x1.n_elem;

  arma::uvec ordX = sort_index(x1);
  arma::uvec ordY = sort_index(y1);

  arma::vec x = x1(ordX);
  arma::vec y = y1(ordY);

  arma::vec ux = unique_from_sorted(x);
  arma::vec uy = unique_from_sorted(y);

  const arma::uword nx = ux.n_elem;
  const arma::uword ny = uy.n_elem;

  arma::uvec gx = match_to_unique(x1, ux);
  arma::uvec gy = match_to_unique(y1, uy);

  arma::mat C(nx, ny, arma::fill::zeros);
  arma::mat Cdx(nx, ny, arma::fill::zeros);
  arma::mat Cdy(nx, ny, arma::fill::zeros);
  arma::mat C11(nx, ny, arma::fill::zeros);

  for (arma::uword k = 0; k < n; ++k) {
    const arma::uword i = gx[k];
    const arma::uword j = gy[k];

    C(i, j) += 1.0;
    if (dx1[k] > 0) Cdx(i, j) += 1.0;
    if (dy1[k] > 0) Cdy(i, j) += 1.0;
    if (dx1[k] > 0 && dy1[k] > 0) C11(i, j) += 1.0;
  }

  // risk(i,j) = #{X >= ux[i], Y >= uy[j]}
  arma::mat risk(nx, ny, arma::fill::zeros);
  for (int i = static_cast<int>(nx) - 1; i >= 0; --i) {
    for (int j = static_cast<int>(ny) - 1; j >= 0; --j) {
      double val = C(i, j);
      if (i + 1 < static_cast<int>(nx)) val += risk(i + 1, j);
      if (j + 1 < static_cast<int>(ny)) val += risk(i, j + 1);
      if (i + 1 < static_cast<int>(nx) && j + 1 < static_cast<int>(ny)) {
        val -= risk(i + 1, j + 1);
      }
      risk(i, j) = val;
    }
  }

  // row_dx(i,j) = #{X = ux[i], delta_x = 1, Y >= uy[j]}
  arma::mat row_dx(nx, ny, arma::fill::zeros);
  for (arma::uword i = 0; i < nx; ++i) {
    double run = 0.0;
    for (int j = static_cast<int>(ny) - 1; j >= 0; --j) {
      run += Cdx(i, j);
      row_dx(i, j) = run;
    }
  }

  // col_dy(i,j) = #{Y = uy[j], delta_y = 1, X >= ux[i]}
  arma::mat col_dy(nx, ny, arma::fill::zeros);
  for (arma::uword j = 0; j < ny; ++j) {
    double run = 0.0;
    for (int i = static_cast<int>(nx) - 1; i >= 0; --i) {
      run += Cdy(i, j);
      col_dy(i, j) = run;
    }
  }

  arma::mat L(nx, ny, arma::fill::zeros);
  const double eps = 1e-12;

  for (arma::uword i = 0; i < nx; ++i) {
    for (arma::uword j = 0; j < ny; ++j) {
      const double r = risk(i, j);
      if (r <= 0.0) continue;

      const double l10 = row_dx(i, j) / r;
      const double l01 = col_dy(i, j) / r;
      const double l11 = C11(i, j) / r;
      const double denom = (1.0 - l10) * (1.0 - l01);

      if (denom <= eps) continue;

      const double val = (l10 * l01 - l11) / denom;
      if (std::isfinite(val)) L(i, j) = val;
    }
  }

  return L;
}

// Exported L matrix. This now returns L on the unique-time grid, not the n x n
// observation grid. This is the correct grid for tied data.
// [[Rcpp::export(rng = false)]]
arma::mat fastL(const arma::vec& x1,
                const arma::vec& y1,
                const arma::vec& dx1,
                const arma::vec& dy1) {
  return fastL_unique_internal(x1, y1, dx1, dy1);
}

// Optional helper for debugging from R: returns unique grids and L.
// [[Rcpp::export(rng = false)]]
Rcpp::List fastL_with_grid(const arma::vec& x1,
                           const arma::vec& y1,
                           const arma::vec& dx1,
                           const arma::vec& dy1) {
  arma::uvec ordX = sort_index(x1);
  arma::uvec ordY = sort_index(y1);

  arma::vec ux = unique_from_sorted(x1(ordX));
  arma::vec uy = unique_from_sorted(y1(ordY));
  arma::mat L = fastL_unique_internal(x1, y1, dx1, dy1);

  return Rcpp::List::create(
    Rcpp::Named("x") = ux,
    Rcpp::Named("y") = uy,
    Rcpp::Named("L") = L
  );
}

// Dabrowska survival surface on the UNIQUE x-y grid.
arma::mat dabrowska_surface_unique(const arma::vec& x_sorted,
                                   const arma::vec& y_sorted,
                                   const arma::vec& dx_sorted,
                                   const arma::vec& dy_sorted,
                                   const arma::mat& L) {
  arma::vec Sx = marginal_surv_unique(x_sorted, dx_sorted);
  arma::vec Sy = marginal_surv_unique(y_sorted, dy_sorted);

  const arma::uword nx = Sx.n_elem;
  const arma::uword ny = Sy.n_elem;

  arma::mat P(nx, ny, arma::fill::ones);

  // P(i,j) = prod_{a <= i, b <= j} {1 - L(a,b)}.
  // This is evaluated once per unique time pair, avoiding repeated factors
  // when x or y has ties.
  for (arma::uword j = 0; j < ny; ++j) {
    double col_prod = 1.0;

    for (arma::uword i = 0; i < nx; ++i) {
      col_prod *= (1.0 - L(i, j));

      if (j == 0) {
        P(i, j) = col_prod;
      } else {
        P(i, j) = P(i, j - 1) * col_prod;
      }
    }
  }

  arma::mat D(nx + 1, ny + 1, arma::fill::ones);

  for (arma::uword i = 0; i < nx; ++i) D(i + 1, 0) = Sx[i];
  for (arma::uword j = 0; j < ny; ++j) D(0, j + 1) = Sy[j];

  for (arma::uword i = 0; i < nx; ++i) {
    for (arma::uword j = 0; j < ny; ++j) {
      D(i + 1, j + 1) = Sx[i] * Sy[j] * P(i, j);
    }
  }

  return D;
}

// Optional helper for debugging from R.
// [[Rcpp::export(rng = false)]]
arma::mat fastD(const arma::vec& x1,
                const arma::vec& y1,
                const arma::vec& dx1,
                const arma::vec& dy1) {
  arma::uvec ordX = sort_index(x1);
  arma::uvec ordY = sort_index(y1);

  arma::vec x = x1(ordX);
  arma::vec y = y1(ordY);
  arma::vec dx = dx1(ordX);
  arma::vec dy = dy1(ordY);

  arma::mat L = fastL_unique_internal(x1, y1, dx1, dy1);
  return dabrowska_surface_unique(x, y, dx, dy, L);
}



// Estimate tie probability term:
//   P{(T1_i - T1_j)(T2_i - T2_j) = 0}
// = P(T1_i = T1_j or T2_i = T2_j)
// = sum_s p1(s)^2 + sum_t p2(t)^2 - sum_{s,t} p12(s,t)^2.
//
// Here p1 and p2 are estimated marginal point masses from the KM margins,
// and p12 is estimated by the Dabrowska rectangle mass on the unique grid.
double tie_prob_from_surface(const arma::mat& D) {
  const arma::uword nx = D.n_rows - 1;
  const arma::uword ny = D.n_cols - 1;

  double px2 = 0.0;
  double py2 = 0.0;
  double pxy2 = 0.0;

  // Marginal point masses: p_X(x_i) = S_X(x_{i-1}) - S_X(x_i),
  // where D(0,0)=1 and D(i,0)=S_X(x_i).
  for (arma::uword i = 1; i <= nx; ++i) {
    const double px = D(i - 1, 0) - D(i, 0);
    if (std::isfinite(px)) px2 += px * px;
  }

  // Marginal point masses for Y.
  for (arma::uword j = 1; j <= ny; ++j) {
    const double py = D(0, j - 1) - D(0, j);
    if (std::isfinite(py)) py2 += py * py;
  }

  // Joint point masses on each unique rectangle.
  for (arma::uword i = 1; i <= nx; ++i) {
    for (arma::uword j = 1; j <= ny; ++j) {
      const double pxy = D(i, j) - D(i, j - 1) - D(i - 1, j) + D(i - 1, j - 1);
      if (std::isfinite(pxy)) pxy2 += pxy * pxy;
    }
  }

  double tie_prob = px2 + py2 - pxy2;

  // Numerical guard only; this should be in [0, 1] theoretically.
  if (!std::isfinite(tie_prob)) return NA_REAL;
  if (tie_prob < 0.0 && tie_prob > -1e-10) tie_prob = 0.0;
  if (tie_prob > 1.0 && tie_prob < 1.0 + 1e-10) tie_prob = 1.0;

  return tie_prob;
}

// [[Rcpp::export(rng = false)]]
double fastTau(const arma::vec& x1,
               const arma::vec& y1,
               const arma::vec& dx1,
               const arma::vec& dy1) {
  const arma::uword n = x1.n_elem;
  if (n == 0) return NA_REAL;

  arma::uvec ordX = sort_index(x1);
  arma::uvec ordY = sort_index(y1);

  arma::vec x = x1(ordX);
  arma::vec y = y1(ordY);
  arma::vec dx = dx1(ordX);
  arma::vec dy = dy1(ordY);

  arma::mat L = fastL_unique_internal(x1, y1, dx1, dy1);
  arma::mat D = dabrowska_surface_unique(x, y, dx, dy, L);

  const arma::uword nx = D.n_rows - 1;
  const arma::uword ny = D.n_cols - 1;

  double out = 0.0;

  for (arma::uword i = 1; i <= nx; ++i) {
    for (arma::uword j = 1; j <= ny; ++j) {
      const double d1 = D(i, j);
      const double inc = D(i, j) - D(i, j - 1) - D(i - 1, j) + D(i - 1, j - 1);
      out += d1 * inc;
    }
  }

  const double tau = 4.0 * out - 1.0;
  return tau;
}


// Tie-adjusted version based on equations (3.2) and (3.3) in your screenshot.
// Returns:
//   tau_base   = 4 * P_hat(T1_i > T1_j, T2_i > T2_j) - 1
//   tie_prob   = P_hat{(T1_i - T1_j)(T2_i - T2_j) = 0}
//   tau_tilde  = tau_base + tie_prob
//   gamma      = tau_tilde / (1 - tie_prob)
//
// If you want fastTau() itself to return the tie-adjusted value, use
// fastTauTieAdjusted(...)["tau_tilde"] in R or replace the return value
// in fastTau() by tau + tie_prob.
// [[Rcpp::export(rng = false)]]
Rcpp::List fastTauTieAdjusted(const arma::vec& x1,
                              const arma::vec& y1,
                              const arma::vec& dx1,
                              const arma::vec& dy1) {
  const arma::uword n = x1.n_elem;
  if (n == 0) {
    return Rcpp::List::create(
      Rcpp::Named("tau_base") = NA_REAL,
      Rcpp::Named("tie_prob") = NA_REAL,
      Rcpp::Named("tau_tilde") = NA_REAL,
      Rcpp::Named("gamma") = NA_REAL
    );
  }

  arma::uvec ordX = sort_index(x1);
  arma::uvec ordY = sort_index(y1);

  arma::vec x = x1(ordX);
  arma::vec y = y1(ordY);
  arma::vec dx = dx1(ordX);
  arma::vec dy = dy1(ordY);

  arma::mat L = fastL_unique_internal(x1, y1, dx1, dy1);
  arma::mat D = dabrowska_surface_unique(x, y, dx, dy, L);

  const arma::uword nx = D.n_rows - 1;
  const arma::uword ny = D.n_cols - 1;

  double out = 0.0;

  for (arma::uword i = 1; i <= nx; ++i) {
    for (arma::uword j = 1; j <= ny; ++j) {
      const double d1 = D(i, j);
      const double inc = D(i, j) - D(i, j - 1) - D(i - 1, j) + D(i - 1, j - 1);
      out += d1 * inc;
    }
  }

  const double tau_base = 4.0 * out - 1.0;
  const double tie_prob = tie_prob_from_surface(D);
  const double tau_tilde = tau_base + tie_prob;

  double gamma = NA_REAL;
  if (std::isfinite(tie_prob) && std::abs(1.0 - tie_prob) > 1e-12) {
    gamma = tau_tilde / (1.0 - tie_prob);
  }

  return Rcpp::List::create(
    Rcpp::Named("tau_base") = tau_base,
    Rcpp::Named("tie_prob") = tie_prob,
    Rcpp::Named("tau_tilde") = tau_tilde,
    Rcpp::Named("gamma") = gamma
  );
}

// Diagnostic version: useful for checking whether the estimated surface gives
// nonnegative increments and stable D values.
// [[Rcpp::export(rng = false)]]
Rcpp::List fastTau_diagnostic(const arma::vec& x1,
                              const arma::vec& y1,
                              const arma::vec& dx1,
                              const arma::vec& dy1) {
  arma::mat D = fastD(x1, y1, dx1, dy1);

  const arma::uword nx = D.n_rows - 1;
  const arma::uword ny = D.n_cols - 1;

  double out = 0.0;
  double min_inc = R_PosInf;
  double max_inc = R_NegInf;
  double sum_inc = 0.0;
  double min_D = R_PosInf;
  double max_D = R_NegInf;

  for (arma::uword i = 1; i <= nx; ++i) {
    for (arma::uword j = 1; j <= ny; ++j) {
      const double d1 = D(i, j);
      const double inc = D(i, j) - D(i, j - 1) - D(i - 1, j) + D(i - 1, j - 1);

      out += d1 * inc;
      sum_inc += inc;
      min_inc = std::min(min_inc, inc);
      max_inc = std::max(max_inc, inc);
      min_D = std::min(min_D, d1);
      max_D = std::max(max_D, d1);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("tau") = 4.0 * out - 1.0,
    Rcpp::Named("min_inc") = min_inc,
    Rcpp::Named("max_inc") = max_inc,
    Rcpp::Named("sum_inc") = sum_inc,
    Rcpp::Named("min_D") = min_D,
    Rcpp::Named("max_D") = max_D,
    Rcpp::Named("nx") = nx,
    Rcpp::Named("ny") = ny
  );
}

// [[Rcpp::export(rng = false)]]
arma::vec pseudoTau(const arma::vec& x1,
                    const arma::vec& y1,
                    const arma::vec& dx1,
                    const arma::vec& dy1) {
  const arma::uword n = x1.n_elem;
  arma::vec out(n);

  if (n <= 1) {
    out.fill(NA_REAL);
    return out;
  }

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec x0(n - 1), y0(n - 1), dx0(n - 1), dy0(n - 1);

    if (i > 0) {
      x0.subvec(0, i - 1) = x1.subvec(0, i - 1);
      y0.subvec(0, i - 1) = y1.subvec(0, i - 1);
      dx0.subvec(0, i - 1) = dx1.subvec(0, i - 1);
      dy0.subvec(0, i - 1) = dy1.subvec(0, i - 1);
    }

    if (i + 1 < n) {
      x0.subvec(i, n - 2) = x1.subvec(i + 1, n - 1);
      y0.subvec(i, n - 2) = y1.subvec(i + 1, n - 1);
      dx0.subvec(i, n - 2) = dx1.subvec(i + 1, n - 1);
      dy0.subvec(i, n - 2) = dy1.subvec(i + 1, n - 1);
    }

    out[i] = fastTau(x0, y0, dx0, dy0);
  }

  return out;
}
