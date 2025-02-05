% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/km.R
\name{km}
\alias{km}
\title{Kaplan-Meier estimate}
\usage{
km(time, status)
}
\arguments{
\item{time}{A numeric vector for the observed survival times.}

\item{status}{A numeric vector for the event indicator;
0 indicates right-censoring and 1 indicates events.}
}
\value{
A data frame with the Kaplan-Meier survival estimates, containing:
  \item{time}{Time points at which the survival probability is estimated.}
  \item{surv}{Estimated survival probability at each time point.}
}
\description{
This function exclusively returns the Kaplan-Meier survival estimate and the corresponding time points.
It does not provide standard errors or any additional outputs
that are typically included with the \code{survfit()} function.
}
\examples{
data(Teeth500)
km(Teeth500$time, Teeth500$event)
}
