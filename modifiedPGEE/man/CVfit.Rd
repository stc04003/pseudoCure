\name{CVfit}
\alias{CVfit}
\alias{print.CVfit}

\title{Function to compute cross-validated tuning parameter value}

\description{This function computes cross-validated tuning parameter value for 
longitudinal data with working independence structure.
}

\usage{
CVfit(formula, id, data, family, scale.fix, scale.value, fold, lambda.vec, pindex,  
eps, maxiter, tol)}

\arguments{
\item{formula}{A formula expression in the form of \code{response ~ predictors}.}
\item{id}{A vector for identifying subjects/clusters.}
\item{data}{A data frame which stores the variables in \code{formula} with \code{id} variable.}
\item{scale.fix}{A logical variable; if true, the scale parameter is fixed at the value of \code{scale.value}. 
The default value is \code{TRUE}.}
\item{scale.value}{If \code{scale.fix = TRUE}, this assignes a numeric value to which the scale parameter should be 
fixed. The default value is 1.}
\item{family}{A \code{family} object in \code{\link{PGEE}}.}
\item{fold}{The number of folds used in cross-validation.}
\item{lambda.vec}{A vector of tuning parameters that will be used in the cross-validation.}
\item{pindex}{An index vector showing the parameters which are not subject to penalization. The default value
 is \code{NULL}. However, in case of a model with intercept, the intercept parameter should be never penalized.}
\item{eps}{A numerical value for the epsilon used in minorization-maximization algorithm. The default value is 
\code{10^-6}.}
\item{maxiter}{The number of iterations that is used in the estimation algorithm. The default value is \code{25}.}
\item{tol}{The tolerance level that is used in the estimation algorithm. The default value is \code{10^-3}.} 
}

\value{An object class of \code{CVfit}.}

\references{
Wang, L., Zhou, J., and Qu, A. (2012). Penalized generalized estimating equations 
for high-dimensional longitudinal data analysis. \emph{Biometrics}, \bold{68}, 353--360.
}

\seealso{
\code{\link{PGEE}}
}

