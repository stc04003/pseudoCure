#' Cure Rate Model with pseudo-observation approach
#'
#' Fits either a mixture cure model or a bounded cumulative hazard (promotion time) model
#' with pseudo-observation approach.
#'
#' @param formula1 A formula object starting with \code{~} for the model formula.
#' This specifies the covariates in the incidence component and the long-term component under
#' the mixture cure model and the bounded cumulative model, respectively.
#' @param formula2 A formula object starting with \code{~} for the model formula.
#' This specifies the covariates in the latency component and the short-term component under
#' the mixture cure model and the bounded cumulative model, respectively.
#' @param time A numeric vector for the observed survival times.
#' @param status A numeric vector for the event indicator;
#' 0 indicates right-censoring and 1 indicates events.
#' @param data An optional data frame that contains the covariates and response variables
#' (\code{time} and \code{event}).
#' @param subset An optional logical vector specifying a subset of
#' observations to be used in the fitting process.
#' @param model A character string specifying the underlying model.
#' The available functional form are \code{"mixture"} and \code{"promotion"}
#' correspond to the mixture cure model and the bounded cumulative model, respectively.
#' @param t0 A vector of times, where the pseudo-observations are constructed.
#' When not specified, the default values are the 10, 20, ..., 90, 95 percentiles of
#' uncensored event times.
#' @param lambda1,lambda2  An option for specifying the tuning parameter used in penalization.
#' When this is unspecified or has a \code{NULL} value,
#' penalization will not be applied and \code{pCure()} will uses all covariates
#' specified in the formulas.
#' Alternatively, this can be specified as a vector numeric vector of non-negative values
#' or "auto" for auto selection.
#' @param penalty1,penalty2 A character string specifying the penalty function.
#' The available options are code{"scad"} and code{"lasso"}.
#' @param exclude1,exclude2 A character string specifying variables to
#' exclude in variable selection.
#' @param nfolds An optional integer value specifying the number of folds.
#' The default value is 5. 
#' @param control A list of control parameters. See detail.
#'
#' @importFrom stats model.frame model.matrix model.extract 
#' @importFrom stats .getXlevels formula pnorm quantile runif sd symnum
#'
#' @example inst/examples/ex_pCure.R
#' @export
pCure <- function(formula1, formula2, time, status, data, subset, t0, 
                  model = c("mixture", "promotion"), nfolds = 5,
                  lambda1 = NULL, exclude1 = NULL, penalty1 = c("scad", "lasso"), 
                  lambda2 = NULL, exclude2 = NULL, penalty2 = c("scad", "lasso"),
                  control = list()) {
    model <- match.arg(model)
    penalty1 <- match.arg(penalty1)
    penalty2 <- match.arg(penalty2)
    ## Checks and define control
    if (is.null.missing(formula1) & is.null.missing(formula2))
      stop("At least one 'formula' need to be specified.")
    if (missing(time)) stop("Argument 'time' is required.")
    if (missing(status)) stop("Argument 'status' is required.")
    if (!is.null(lambda1) && !is.character(lambda1) && any(lambda1 < 0))
        stop("Positive tuning parameters ('lambda1') is required.")
    if (!is.null(lambda2) && !is.character(lambda2) && any(lambda2 < 0))
        stop("Positive tuning parameters ('lambda2') is required.")
    if (missing(data)) {
      if (is.null.missing(formula1)) data <- environment(formula1)
      else data <- environment(formula2)
    }
    if (!missing(subset)) {
        sSubset <- substitute(subset)
        subIdx <- eval(sSubset, data, parent.frame())
        if (!is.logical(subIdx)) stop("'subset' must be logical")
        subIdx <- subIdx & !is.na(subIdx)
        data <- data[subIdx, ]
    }
    ctrl <- pCure.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    ctrl$model <- model
    ctrl$lambda1 <- lambda1
    ctrl$lambda2 <- lambda2
    ctrl$penalty1 <- penalty1
    ctrl$penalty2 <- penalty2
    ctrl$exclude1 <- exclude1
    ctrl$exclude2 <- exclude2
    ctrl$nfolds <- nfolds
    if (is.null.missing(formula1)) ctrl$formula1 <- NULL
    else ctrl$formula1 <- formula1
    if (is.null.missing(formula2)) ctrl$formula2 <- NULL
    else ctrl$formula2 <- formula2
    if (is.character(lambda1) && lambda1 != "auto")    
        stop("Only 'auto' is allowed when 'lambda1' is a character string.")
    if (is.character(lambda2) && lambda2 != "auto")    
      stop("Only 'auto' is allowed when 'lambda2' is a character string.")
    mf <- match.call(expand.dots = FALSE)
    if (is.null.missing(formula1)) 
      mf <- mf[match(c("formula2", "data", "time", "status"), names(mf), 0L)]
    else
      mf <- mf[match(c("formula1", "data", "time", "status"), names(mf), 0L)]
    mf$data <- data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    xlevel1 <- xlevel2 <- .getXlevels(attr(mf, "terms"), mf)
    time <- as.numeric(model.extract(mf, time))
    status <- as.numeric(model.extract(mf, status))
    if (is.null.missing(formula1)) {   
      mm1 <- NULL
      mm2 <- stats::model.matrix(formula2, data = mf)
      mm2 <- mm2[, colnames(mm2) != "(Intercept)"]   
    } else {
      mm1 <- mm2 <- stats::model.matrix(formula1, data = mf)
      if (is.null.missing(formula2)) mm2 <- NULL
    }
    if (!is.null.missing(formula1) && !is.null.missing(formula2)) {
      mf <- match.call(expand.dots = FALSE)
      mf <- mf[match(c("formula2", "data"), names(mf), 0L)]
      mf$data <- data
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval(mf, parent.frame())
      xlevel2 <- .getXlevels(attr(mf, "terms"), mf)
      mm2 <- stats::model.matrix(formula2, data = mf)
      mm2 <- mm2[, colnames(mm2) != "(Intercept)"]   
    }    
    tmax <- max(time[status > 0])
    if (missing(t0)) t0 <- quantile(time[status > 0], c(1:9 / 10, .95))
    ## auto choose lambda; still under development
    if (is.character(lambda1) || is.character(lambda2))
        tmp <- auto.lambda(mm1, mm2, time, status, t0, ctrl)
    if (is.character(lambda1) && lambda1 == "auto") ctrl$lambda1 <- tmp$lambda1
    if (is.character(lambda2) && lambda2 == "auto") ctrl$lambda2 <- tmp$lambda2
    ## Fit models
    if (model == "mixture")
        out <- fitPHMC(mm1, mm2, time, status, t0, ctrl)
    if (model == "promotion")
        out <- fitPHPH(mm1, mm2, time, status, t0, ctrl)
    out$t0 <- t0
    out$call <- match.call()
    out$xlevel1 <- xlevel1
    out$xlevel2 <- xlevel2
    out$control <- ctrl
    out <- out[order(names(out))]
    class(out) <- "pCure"
    return(out)
}

#' Package options for pseudoCure
#'
#' This function provides the fitting options for the \code{pCure()} function.
#'
#' @param binit1 Initial value for the first component.
#' A zero vector will be used if not specified.
#' @param binit2 Initial value for the second component
#' A zero vector will be used if not specified.
#' @param corstr A character string specifying the correlation structure.
#' The following are permitted: \code{"independence"}, \code{"exchangeable"},
#' and \code{"ar1"}.
#' @param nlambda1,nlambda2 An integer value specifying the number of lambda.
#' This is only evoked when \code{lambda1 = "auto"} or \code{lambda2 = "auto"}.
#' @param eps A positive numerical value used to ...
#' @param tol A positive numerical value specifying the absolute error tolerance in GEE algorithms.
#' @param maxit An integer value specifying the maximum number of iteration.
#' 
#' @seealso \code{\link{pCure}}
#' @export
pCure.control <- function(binit1 = NULL, binit2 = NULL,
                          corstr = c("independence", "exchangeable", "ar1"),
                          nlambda1 = 100, nlambda2 = 100,
                          eps = 1e-6, tol = 1e-7, maxit = 100) {
    corstr <- match.arg(corstr)
    list(binit1 = binit1, binit2 = binit2, corstr = corstr,
         nlambda1 = 100, nlambda2 = 100,
         eps = eps, tol = tol, maxit = maxit)
}


#' Auto select lambda to apply CV
#'
#' Still under development
#' 
#' @noRd
auto.lambda <- function(mm1, mm2, time, status, t0, ctrl) {
    trys <- exp(-7:10)
    lambda1.max <- lambda2.max <- max(trys)
    for (i in 1:length(trys)) {
        ctrl$lambda1 <- ctrl$lambda2 <- trys[i]
        if (ctrl$model == "mixture")
            tmp <- fitPHMC(mm1, mm2, time, status, t0, ctrl)
        if (ctrl$model == "promotion")
            tmp <- fitPHPH(mm1, mm2, time, status, t0, ctrl)
        if (all(abs(tmp$fit1$b) < 1e-6)) lambda1.max <- trys[i]
        if (all(abs(tmp$fit2$b) < 1e-6)) lambda2.max <- trys[i]
        if (max(lambda1.max, lambda2.max) < max(trys)) break
    }
    list(lambda1 = exp(seq(log(1e-04), log(lambda1.max), length.out = ctrl$nlambda1)),
         lambda2 = exp(seq(log(1e-04), log(lambda2.max), length.out = ctrl$nlambda2)))
}

is.null.missing <- function(x) missing(x) || is.null(x)
