\name{PGEE-package}
\alias{PGEE-package}
%\alias{PGEE}
\docType{package}
\title{Penalized Generalized Estimating Equations}

\description{This package fits penalized generalized estimating equations to longitudinal
data with high-dimensional covariates through accommodating SCAD-penalty function into
generalized estimating equations. 
}

\details{
This package consists of three functions. The function
\code{\link{PGEE}} fits penalized generalized estimating equations to the data. But, before that,
the tuning parameter should be estimated through the function \code{\link{CVfit}}. 
On the other hand, the function \code{\link{MGEE}} fits unpenalized generalized estimating 
equations to the data.
}

\author{
Gul Inan, Jianhui Zhou and Lan Wang

Maintainer: Gul Inan
}

\references{
Wang, L., Zhou, J., and Qu, A. (2012). Penalized generalized estimating equations 
for high-dimensional longitudinal data analysis. \emph{Biometrics}, \bold{68}, 353--360.
}

\keyword{high-dimensional covariates, longitudinal data, marginal models, SCAD-penalty function}

