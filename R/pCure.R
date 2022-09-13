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
#' When this is unspecified or has a code{NULL} value,
#' penalization will not be applied and the \code{pCure} fit will uses all covariates
#' specified in the formulas.
#' Alternatively, this can be specified as a numeric vector of non-negative values.
#' @param penality1,penality2 A character string specifying the penality function.
#' The available options are code{"scad"} and code{"lasso"}.
#' @param exclude1,exclude2 A character string specifying 
#' @param control A list of control parameters. See detail.
#'
#' @importFrom stats model.frame model.matrix
#'
#' @export
pCure <- function(formula1, formula2, time, status, data, subset, t0, 
                  model = c("mixture", "promotion"),
                  lambda1 = NULL, lambda2 = NULL,
                  penalty1 = c("scad", "lasso"), penalty2 = c("scad", "lasso"),
                  exclude1 = NULL, exclude2 = NULL, 
                  control = list()) {
    model <- match.arg(model)
    penalty1 <- match.arg(penalty1)
    penalty2 <- match.arg(penalty2)
    if (missing(formula1)) stop("Argument 'formula' is required.")
    if (missing(time)) stop("Argument 'time' is required.")
    if (missing(status)) stop("Argument 'status' is required.")
    if (missing(data)) data <- environment(formula)
    if (!missing(subset)) {
        sSubset <- substitute(subset)
        subIdx <- eval(sSubset, data, parent.frame())
        if (!is.logical(subIdx)) stop("'subset' must be logical")
        subIdx <- subIdx & !is.na(subIdx)
        data <- data[subIdx, ]
    }
    mf <- match.call(expand.dots = FALSE)    
    mf <- mf[match(c("formula1", "data", "time", "status"), names(mf), 0L)]
    mf$data <- data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    xlevel1 <- xlevel2 <- .getXlevels(attr(mf, "terms"), mf)
    time <- as.numeric(model.extract(mf, time))
    status <- as.numeric(model.extract(mf, status))
    mm1 <- mm2 <- stats::model.matrix(formula1, data = mf)
    if (!missing(formula2)) {
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
    ctrl <- pCure.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (model == "mixture")
        out <- fitPHMC(mm1, mm2, time, status, t0,
                       lambda1, penalty1, exclude1, lambda2, penalty2, exclude2, ctrl)
    if (model == "promotion")
        out <- fitPHPH(mm1, mm2, time, status, t0,
                       lambda1, penalty1, exclude1, lambda2, penalty2, exclude2, ctrl)
    out$t0 <- t0
    out$call <- match.call()
    out$xlevel1 <- xlevel1
    out$xlevel2 <- xlevel2
    class(out) <- "pCure"
    out <- out[order(names(out))]
    return(out)
}


#' Pakcage options for pseudoCure
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
#' @param eps A positive numerical value used to ...
#' @param tol A positive numerical value specifying the absolute error tolerance in GEE algorithms.
#' @param maxit An integer value specifying the maximum number of iteration.
#' 
#' @seealso \code{\link{pCure}}
#' @export
pCure.control <- function(binit1 = NULL, binit2 = NULL,
                          corstr = c("independence", "exchangeable", "ar1"),
                          eps = 1e-6, tol = 1e-7, maxit = 100) {
    corstr <- match.arg(corstr)
    list(binit1 = binit1, binit2 = binit2, corstr = corstr,
         eps = eps, tol = tol, maxit = maxit)
}
