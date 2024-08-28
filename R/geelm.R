#' Generalized Estimating Equation with Gaussian family
#'
#' Fits a generalized estimating equation (GEE) model with
#' Gaussian family with different link functions.
#' The \code{geelm} function also supports LASSO or SCAD
#' regularization.
#' 
#'
#' @param formula  A formula object starting with \code{~} for the model formula.
#' @param data An optional data frame that contains the covariates and response variables.
#' @param subset An optional logical vector specifying a subset of
#' observations to be used in the fitting process.
#' @param id A vector which identifies the clusters.
#' If not specified, each observation is treated as its own cluster.
#' @param link A character string specifying the model link function. Available options are
#' \code{"identity"}, \code{"log"}, \code{"cloglog"}, and \code{"logit"}.
#' @param corstr A character string specifying the correlation structure.
#' Available options are \code{"independence"}, \code{"exchangeable"}, and \code{"ar1"}.
#' @param lambda An option for specifying the tuning parameter used in penalization.
#' When this is unspecified or has a \code{NULL} value,
#' penalization will not be applied and \code{pCure()} will uses all covariates
#' specified in the formulas.
#' Alternatively, this can be specified as a vector numeric vector of non-negative values
#' or "auto" for auto selection.
#' @param nfolds An optional integer value specifying the number of folds.
#' The default value is 5.
#' @param nlambda An optional integer value specifying the number of tuning parameters to try
#' if \code{lambda = "auto"}.
#' @param exclude A binary numerical vector specifying which variables to exclude in variable selection.
#' The length of \code{exclude} must match with the number of covariates.
#' A value of 1 means to exclude in the variable selection.
#' @param penalty A character string specifying the penalty function.
#' The available options are \code{"lasso"} and \code{"scad"}.
#' @param binit A optional numerical vector for the initial value.
#' A zero vector is used when not specified. 
#' @param tol A positive numerical value specifying the absolute
#' error tolerance in root search. Default at \code{1e-7}.
#' @param maxit A positive integer specifying the maximum number of iteration.
#' Default at 100.
#'
#' @importFrom rlang f_rhs
#' @export
geelm <- function(formula, data, subset, id, 
                 link = c("identity", "log", "cloglog", "logit"),
                 corstr = c("independence", "exchangeable", "ar1"), 
                 lambda, nfolds = 5, nlambda = 100, exclude, 
                 penalty = c("lasso", "scad"),
                 binit,
                 tol = 1e-7, maxit = 100) {  
  link <- match.arg(link)
  corstr <- match.arg(corstr)
  penalty <- match.arg(penalty)
  fExcl <- deparse(substitute(id))
  if (missing(formula))
    stop("A 'formula' needs to be specified.")
  if (!missing(lambda) && !is.character(lambda) && any(lambda < 0))
    stop("Positive tuning parameters is required.")
  if (missing(lambda)) lambda <- 0
  if (missing(data)) data <- environment(formula)
  if (!missing(subset)) {
    sSubset <- substitute(subset)
    subIdx <- eval(sSubset, data, parent.frame())
    if (!is.logical(subIdx)) stop("'subset' must be logical")
    subIdx <- subIdx & !is.na(subIdx)
    data <- data[subIdx, ]
  }
  if (is.character(lambda) && lambda != "auto")    
    stop("Only 'auto' is allowed when 'lambda' is a character string.")
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[match(c("formula", "data", "id"), names(mf), 0L)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  id <- model.extract(mf, id)  
  y <- as.numeric(model.response(mf, "numeric"))
  xlevel <- .getXlevels(attr(mf, "terms"), mf)
  mm <- stats::model.matrix(formula, data = mf)
  if (any(as.character(f_rhs(formula)) == "."))
    mm <- mm[,!sapply(colnames(mm), grepl, pattern = fExcl)]
  n <- nrow(mm); p <- ncol(mm)
  if (missing(binit)) binit <- rep(0, p)
  if (missing(id) || is.null(id)) id <- 1:nrow(mm)
  if (missing(lambda)) lambda <- 0
  if (missing(exclude)) exclude <- rep(0, p)
  if (all(lambda == 0))
    fit <- gee(y, mm, binit, as.numeric(table(id)), link, corstr, tol, maxit)
  if (is.character(lambda) && lambda == "auto") {
    trys <- exp(-7:10)
    lambda.max <- max(trys)
    for (i in 1:length(trys)) {
      tmp <- pgee(y, mm, binit, as.numeric(table(id)), exclude, link,
                  penalty, corstr, trys[i], 1e-6, tol, maxit)
      if (all(abs(tmp$b) < 1e-6)) lambda.max <- trys[i]
      if (lambda.max < max(trys)) break
    }
    lambda <- exp(seq(log(1e-04), log(lambda.max), length.out = nlambda))
  }
  if (length(lambda) == 1 && lambda > 0)
    fit <- pgee(y, mm, binit, as.numeric(table(id)), exclude, link,
                penalty, corstr, lambda, 1e-6, tol, maxit)
  if (length(lambda) > 1) {
    cv.raw <- pgeeCV(y, mm, binit, as.numeric(table(id)), exclude, link,
                     penalty, corstr, nfolds, lambda, 1e-6, tol, maxit)
    cv.mean <- colMeans(cv.raw)
    cv.sd <- apply(cv.raw, 2, sd) / nfolds
    cv.1se <- cv.mean + cv.sd
    lambda.min <- lambda[which.min(cv.mean)]
    lambda.1se <- lambda[which.max(which(cv.mean < cv.1se[which.min(cv.mean)]))]
    fit <- pgee(y, mm, binit, as.numeric(table(id)), exclude, link,
                     penalty, "ind", lambda.min, 1e-6, tol, maxit)
    fit$cv.raw <- cv.raw
    fit$lambda.min <- lambda.min
    fit$lambda.1se <- lambda.1se    
  }
  fit$b <- drop(fit$b)
  names(fit$b) <- colnames(mm)
  fit$vb <- with(fit, ginv(H + n * E) %*% M %*% ginv(H + n * E))
  fit$fitted <- drop(exp(mm %*% fit$b))
  fit$resid <- drop(y - fit$fitted)
  fit$call <- match.call()
  fit$xlevel <- xlevel
  fit$link <- link
  fit$corstr <- corstr
  fit$exclude <- exclude
  fit$lambda <- lambda  
  fit <- fit[order(names(fit))]
  class(fit) <- "geelm"
  return(fit)  
}


## S3 methods

#' Check class
#' @noRd
is.geelm <- function(x) inherits(x, "geelm")

#' @exportS3Method print geelm
print.geelm <- function(x, cutoff = 1e-3, ...) {
  if (!is.geelm(x)) stop("Must be a geelm object")
  cat("Call: \n")
  dput(x$call)
  cat("\nMean Model:\n")
  cat(" Mean Link:                ", x$link, "\n")
  cat(" Variance to Mean Relation: gaussian \n")
  mat1 <- t(format(x$b, digits = 5))
  if (!is.null(x$lambda)) mat1[abs(x$b) < cutoff] <- "."
  cat("\n Coefficients:\n")
  prmatrix(mat1, rowlab = rep("", nrow(mat1)), quote = FALSE)
}

#' @exportS3Method summary geelm
summary.geelm <- function(object, cutoff = 1e-3,...) {
  if (!is.geelm(object)) stop("Must be a geelm object")
  out <- list(call = object$call)
  if (is.null(object$lambda)) cutoff <- 0
  out$tab <- pvalTab(object$b, sqrt(diag(object$vb)))
  out$exclude <- abs(out$tab[,1]) < cutoff
  out$lambda.min <- object$lambda.min
  out$lambda.1se <- object$lambda.1se
  out$link <- object$link
  class(out) <- "summary.geelm"
  return(out)
}

#' @exportS3Method print summary.geelm
print.summary.geelm <- function(x, ...) {
  cat("Call: \n")
  dput(x$call)
  cat("\nMean Model:\n")
  cat(" Mean Link:                ", x$link, "\n")
  cat(" Variance to Mean Relation: gaussian \n")
  cat("\n Coefficients:\n")    
  if (!is.null(x$lambda.min))
    cat("Tuning parameters selected by minimum prediction error:", x$lambda.min, "\n\n")
  printCoefmat2(x$tab, exclude = x$exclude)
  cat("\n")
}


#' @exportS3Method coef geelm
coef.geelm <- function(object, ...) {
  object$b
}

#' @exportS3Method vcov geelm
vcov.geelm <- function(object, ...) {
  object$vb
}

#' @exportS3Method residuals geelm
residuals.geelm <- function(object, ...) {
  object$resid
}

#' @exportS3Method fitted geelm
fitted.geelm <- function(object, ...) {
  object$fitted
}

