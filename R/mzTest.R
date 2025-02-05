#' Maller-Zhou test
#'
#' Performs the Maller-Zhou test.
#'
#' @param time A numeric vector for the observed survival times.
#' @param status A numeric vector for the event indicator;
#' 0 indicates right-censoring and 1 indicates events.
#'
#' @return A list containing the Maller-Zhou test results,
#' including the test statistic, p-value, and the number of observed events.
#' 
#' 
#' @example inst/examples/ex_mzTest.R
#' @export
mzTest <- function(time, status) {
  out <- testMZ(time, status) 
  class(out) <- "pMZ"
  attr(out, "n") <- length(time)
  attr(out, "events") <- sum(status)
  return(out)
}

#' Check class
#' @noRd
is.pMZ <- function(x) inherits(x, "pMZ")

#' @exportS3Method print pMZ
print.pMZ <- function(x, ...) {
  if (!is.pMZ(x)) stop("Must be a pMZ object")
  print(data.frame(n = attr(x, "n"),
                   events = attr(x, "events"),
                   statistic = x[[1]],
                   p.value = x[[2]]), row.names = FALSE)
}
