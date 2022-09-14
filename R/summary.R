#' Check class
#' @noRd
is.pCure <- function(x) inherits(x, "pCure")

#' @exportS3Method print pCure
print.pCure <- function(x, cutoff = 1e-3, ...) {
  if (!is.pCure(x)) stop("Must be a pCure object")
  cat("Call: \n")
  dput(x$call)
  if (all(is.null(x$control$lambda1), is.null(x$control$lambda2))) {
    if (x$control$model == "mixture")
      cat("\nFitted with the mixture cure rate model:")
    if (x$control$model == "promotion")
      cat("\nFitted with the promotion time model:")
  } else {
    if (x$control$model == "mixture")
      cat("\nFitted with the penalized mixture cure rate model:")
    if (x$control$model == "promotion")
      cat("\nFitted with the penalized promotion time model:")
  }
  mat1 <- t(format(x$fit1$b, digits = 5))
  mat2 <- t(format(x$fit2$b[-(1:length(x$t0))], digits = 5))
  if (!is.null(x$control$lambda1)) mat1[abs(x$fit1$b) < cutoff] <- "."
  if (!is.null(x$control$lambda2)) mat2[abs(x$fit2$b[-(1:length(x$t0))]) < cutoff] <- "."  
  if (x$control$model == "mixture") cat("\nIncidence component:\n")
  else cat("\nLong-term effect:\n")
  prmatrix(mat1, rowlab = rep("", nrow(mat1)), quote = FALSE)
  if (x$control$model == "mixture") cat("\nLatency component:\n")
  else cat("\nShort-term effect:\n")  
  prmatrix(mat2, rowlab = rep("", nrow(mat2)), quote = FALSE)
}

#' @exportS3Method summary pCure
summary.pCure <- function(object, cutoff = 1e-3,...) {
  if (!is.pCure(object)) stop("Must be a pCure object")
  out <- list(call = object$call, model = object$model, control = object$control)
  cutoff1 <- cutoff2 <- cutoff
  if (is.null(object$control$lambda1)) cutoff1 <- 0
  if (is.null(object$control$lambda2)) cutoff2 <- 0
  out$tab1 <- pvalTab(object$fit1$b, sqrt(diag(object$fit1$vb)))
  out$tab2 <- pvalTab(object$fit2$b, sqrt(diag(object$fit2$vb)))
  out$tab2 <- out$tab2[-(1:length(object$t0)),]
  out$exclude1 <- abs(out$tab1[,1]) < cutoff1
  out$exclude2 <- abs(out$tab2[,1]) < cutoff2 
  class(out) <- "summary.pCure"
  return(out)
}

pvalTab <- function(pe, se, names = NULL) {
  if (is.null(se)) se <- NA
  tab <- cbind(Estimate = pe, StdErr = se,
               z.value = pe / se, p.value = 2 * pnorm(-abs(pe / se)))
  if (!is.null(names)) rownames(tab) <- names
  return(tab)
}

#' @exportS3Method print summary.pCure
print.summary.pCure <- function(x, ...) {
  cat("Call: \n")
  dput(x$call)
  if (all(is.null(x$control$lambda1), is.null(x$control$lambda2))) {
    if (x$control$model == "mixture")
      cat("\nFitted with the mixture cure rate model:")
    if (x$control$model == "promotion")
      cat("\nFitted with the promotion time model:")
  } else {
    if (x$control$model == "mixture")
      cat("\nFitted with the penalized mixture cure rate model:")
    if (x$control$model == "promotion")
      cat("\nFitted with the penalized promotion time model:")
  }

  if (x$control$model == "mixture") cat("\nIncidence component:\n")
  else cat("\nLong-term effect:\n")
  printCoefmat2(x$tab1, exclude = x$exclude1)
  if (x$control$model == "mixture") cat("\nLatency component:\n")
  else cat("\nShort-term effect:\n")
  printCoefmat2(x$tab2, exclude = x$exclude2)
  cat("\n")
}