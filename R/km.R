#' Kaplan-Meier estimate
#'
#' This function exclusively returns the Kaplan-Meier survival estimate and the corresponding time points.
#' It does not provide standard errors or any additional outputs
#' that are typically included with the \code{survfit()} function.
#'
#' @param time A numeric vector for the observed survival times.
#' @param status A numeric vector for the event indicator;
#' 0 indicates right-censoring and 1 indicates events.
#'
#' @return A data frame with the Kaplan-Meier survival estimates, containing:
#'   \item{time}{Time points at which the survival probability is estimated.}
#'   \item{surv}{Estimated survival probability at each time point.}
#'
#' @importFrom ggplot2 geom_step
#' @example inst/examples/ex_km.R
#' @export
km <- function(time, status) {
  s <- fastKM(time, status) 
  out <- data.frame(time = s[,1], surv = s[,2])
  class(out) <- "pKM"
  attr(out, "events") <- sum(status)
  attr(out, "row.names") <- NULL
  return(out)
}

#' Check class
#' @noRd
is.pKM <- function(x) inherits(x, "pKM")

#' @exportS3Method print pKM
print.pKM <- function(x, ...) {
  if (!is.pKM(x)) stop("Must be a pKM object")
  n <- length(x[[1]])
  qs <- x$time[n - findInterval(1:3 / 4, rev(x$surv)) + 1]
  print(data.frame(n = n, events = attr(x, "events"),
                   Q1 = qs[1], median = qs[2], Q3 = qs[3]), row.names = FALSE)
}

#' @exportS3Method plot pKM
plot.pKM <- function(x, ...) {
  if (!is.pKM(x)) stop("Must be a pKM object")
  ggplot(NULL, aes(x = x$time, y = x$surv)) + geom_step() +
    ylab("Survival probability") + xlab("Time")
}

#' @exportS3Method quantile pKM
quantile.pKM <- function(x, probs = c(0.25, 0.5, 0.75), ...) {
  if (!is.pKM(x)) stop("Must be a pKM object")
  n <- length(x[[1]])
  qs <- x$time[n - findInterval(1 - probs, rev(x$surv)) + 1]
  d <- data.frame(t(qs))
  names(d) <- paste0(100 * probs, "%")  
  print(d, row.names = F)
}
