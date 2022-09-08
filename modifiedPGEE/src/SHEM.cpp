#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//' @noRd
// [[Rcpp::export]]
Rcpp::List SHM_CppArma( arma::vec y, arma::vec mu, arma::mat X,
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
