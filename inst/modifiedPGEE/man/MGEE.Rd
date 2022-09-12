\name{MGEE}
\alias{MGEE}
\alias{MGee}
\alias{mycor_gee1}
\alias{S_H_M}
\alias{print.MGEE}
\alias{print.summary.MGEE}
\alias{summary.MGEE}
\alias{print.MGee}
\alias{print.summary.MGee}
\alias{summary.MGee}

\title{Function to fit generalized estimating equations}

\description{This function fits a generalized estimating equation model to longitudinal data.}

\usage{
MGEE(formula, id, data, na.action = NULL, family = gaussian(link = "identity"), 
corstr = "independence", Mv = NULL, beta_int = NULL, R = NULL, scale.fix = TRUE, 
scale.value = 1, maxiter = 25, tol = 10^-3, silent = TRUE)
}

\arguments{
\item{formula}{A formula expression in the form of \code{response ~ predictors}.}
\item{id}{A vector for identifying subjects/clusters.}
\item{data}{A data frame which stores the variables in \code{formula} with \code{id} variable.}
\item{na.action}{A function to remove missing values from the data. Only \code{na.omit} is allowed here.}
\item{family}{A \code{family} object: a list of functions and expressions for defining \code{link} and 
\code{variance} functions. Families supported in \code{MGEE} are \code{binomial}, \code{gaussian}, \code{gamma} and 
\code{poisson}. The \code{links}, which are not available in \code{gee}, is not available here. The default family 
is \code{gaussian}.} 
\item{corstr}{A character string, which specifies the type of correlation structure. 
Structures supported in \code{MGEE} are \code{"AR-1"},\code{"exchangeable"}, \code{"fixed"}, \code{"independence"},
\code{"stat_M_dep"},\code{"non_stat_M_dep"}, and \code{"unstructured"}. The default \code{corstr} type is 
\code{"independence"}.} 
\item{Mv}{If either \code{"stat_M_dep"}, or \code{"non_stat_M_dep"} is specified in \code{corstr}, then this assigns 
a numeric value for \code{Mv}. Otherwise, the default value is \code{NULL}.}
\item{beta_int}{User specified initial values for regression parameters. The default value is \code{NULL}.}
\item{R}{If \code{corstr = "fixed"} is specified, then \code{R} is a square matrix of dimension maximum cluster 
size containing the user specified correlation. Otherwise, the default value is \code{NULL}.}
\item{scale.fix}{A logical variable; if true, the scale parameter is fixed at the value of \code{scale.value}. 
The default value is \code{TRUE}.}
\item{scale.value}{If \code{scale.fix = TRUE}, this assignes a numeric value to which the scale parameter should be fixed. 
                   The default value is 1.}
\item{maxiter}{The number of iterations that is used in the estimation algorithm. The default value is \code{25}.}
\item{tol}{The tolerance level that is used in the estimation algorithm. The default value is \code{10^-3}.}
\item{silent}{A logical variable; if false, the regression parameter estimates at each iteration are
printed. The default value is \code{TRUE}.}
}

\value{An object class of \code{MGEE} representing the fit.}

\references{

Liang, K.Y. and Zeger, S.L. (1986).
  Longitudinal data analysis using generalized linear models.
  \emph{Biometrika}, \bold{73}, 13--22. 

  
Zeger, S.L. and Liang, K.Y. (1986)
.  Longitudinal data analysis for discrete and continuous outcomes.
  \emph{Biometrics}, \bold{42}, 121--130.

}



\note{
The structures \code{"non_stat_M_dep"} and \code{"unstructured"} are valid only when the data is balanced.
}

\seealso{
\code{\link{PGEE}}
}

