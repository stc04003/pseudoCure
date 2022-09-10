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

// assumes equal cluster size
double ahatEx(arma::vec a,
							arma::vec nt,
							arma::vec index) {
	int n = nt.n_elem;
	double out = 0;
	a = a / stddev(a);
	for (int i = 0; i < n; i++) {
		arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
		out += sum(a2) * sum(a2) - sum(a2 % a2);
	}
	return(out / nt(0) / (nt(0) - 1) / n);
}

double ahatAR1(arma::vec a,
							 arma::vec nt,
							 arma::vec index) {
	int n = nt.n_elem;
	double out = 0;
	a = a / stddev(a);
	for (int i = 0; i < n; i++) {
		arma::vec a2 = a(span(index(i), index(i) + nt(i) - 1));
		out += sum(a(span(0, nt(0) - 2)) % a(span(1, nt(0) - 1)));
	}
	return(out / (nt(0) - 1) / n);
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
  arma::vec ym = y - mu; // gaussian()$linkinv
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


// link: 1 = log; 2 = cloglog; 3 = logit
// [[Rcpp::export(rng = false)]]
Rcpp::List SHEM(arma::vec y,
								arma::mat X,
								arma::vec b0,
								arma::cube Rhat,
								arma::vec nt,
								arma::vec pindex, // 0 means to penalize
								std::string glmlink,
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
	Rcpp::List links;
	if (glmlink == "log") links = logLink(eta);
	if (glmlink == "cloglog") links = cloglogLink(eta);
	if (glmlink == "logit") links = logitLink(eta);
	arma::vec mu = links(0);
	arma::vec etamu = links(1);
	arma::vec ym = y - mu;
	arma::mat bigD = matvec(X, etamu);
	for (int i = 0; i < N; i++) {
    arma::vec ym2 = ym(span(index(i), index(i) + nt(i) - 1));
    arma::mat bigD2 = bigD.rows(span(index(i), index(i) + nt(i) - 1));
    arma::mat RhatSlice = Rhat.slice(i);
    arma::mat bigV = RhatSlice(span(0, nt(i) - 1), span(0, nt(i) - 1));
    arma::mat tmp = bigD2.t() * pinv(bigV);
    S += tmp * ym2;
    H += tmp * bigD2;
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
							 arma::vec nt,
							 arma::vec pindex, // 0 means to penalize
							 std::string glmlink,
							 std::string corstr,
							 double lambda, double eps,
							 double tol, int maxit){
  Rcpp::List out(6);
  int N = nt.n_elem;
  int nx = X.n_cols;
	int k = nt(0);
  arma::vec index(N, arma::fill::zeros);
  index(span(1, N - 1)) = cumsum(nt(span(0, N - 2)));
  arma::vec b1 = b0;
  for (int j = 0; j <= maxit; j++) {
    arma::vec eta = X * b0;
    arma::vec E1 = qscad(abs(b0), lambda) / (abs(b0) + eps);
    E1.elem(find(pindex > 0)).zeros();
    arma::mat E = diagmat(E1);
    arma::mat S(nx, 1, arma::fill::zeros);
    arma::mat H(nx, nx, arma::fill::zeros);
    arma::mat M(nx, nx, arma::fill::zeros);
		Rcpp::List links;
		if (glmlink == "log") links = logLink(eta);
		if (glmlink == "cloglog") links = cloglogLink(eta);
		if (glmlink == "logit") links = logitLink(eta);
		arma::vec mu = links(0);
		arma::vec etamu = links(1);
		arma::vec ym = y - mu;
		arma::mat bigD = matvec(X, etamu);
		arma::mat Rhat(k, k, arma::fill::eye);
		double ahat = 0;
		if (corstr == "ex" & k > 1) {
			ahat = ahatEx(ym, nt, index);
			Rhat = Rhat * (1 - ahat) + ahat;
		}
		if (corstr == "ar1" & k > 1) {
			ahat = ahatAR1(ym, nt, index);
			arma::vec tmp(k - 1, arma::fill::value(ahat));
			arma::mat Rhat2(k, k, arma::fill::zeros);
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
      // arma::mat RhatSlice = Rhat.slice(i);
      // arma::mat bigV = RhatSlice(span(0, nt(i) - 1), span(0, nt(i) - 1));
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
    if(min(abs(b1 - b0)) < tol) break;
    b0 = b1;
  }
  out.names() = Rcpp::CharacterVector::create("b", "S", "H", "E", "M", "iter");
  return out;
}
