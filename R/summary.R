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
      cat("\nFitted with the mixture cure model")
    if (x$control$model == "promotion")
      cat("\nFitted with the promotion time model")
  } else {
    if (x$control$model == "mixture")
      cat("\nFitted with the penalized mixture cure model")
    if (x$control$model == "promotion")
      cat("\nFitted with the penalized promotion time model")
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
summary.pCure <- function(object, cutoff = 1e-3, digits = max(3L, getOption("digits") - 2L), ...) {
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
  out$lambda1.min <- object$fit1$lambda1.min
  out$lambda2.min <- object$fit2$lambda2.min
  out$lambda1.1se <- object$fit1$lambda1.1se
  out$lambda2.1se <- object$fit2$lambda2.1se
  out$digits <- digits
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
      cat("\n\nFitted with the mixture cure model")
    if (x$control$model == "promotion")
      cat("\n\nFitted with the promotion time model")
  } else {
    if (x$control$model == "mixture")
      cat("\n\nFitted with the penalized mixture cure model")
    if (x$control$model == "promotion")
      cat("\n\nFitted with the penalized promotion time model")
  }
  if (x$control$model == "mixture") cat("\nIncidence component:\n")
  else cat("\nLong-term effect:\n")
  if (!is.null(x$lambda1.min))
    cat("Tuning parameters selected by minimum prediction error:", x$lambda1.min, "\n\n")
  ## cat("\nTuning parameters selected by 1se rule: ", x$lambda1.1se, "\n\n")
  printCoefmat2(x$tab1, exclude = x$exclude1, digits = x$digits)
  if (x$control$model == "mixture") cat("\n\nLatency component:\n")
  else cat("\n\nShort-term effect:\n")
  if (!is.null(x$lambda2.min))
    cat("Tuning parameters selected by minimum prediction error:", x$lambda2.min, "\n\n")
  ## cat("\nTuning parameters selected by 1se rule: ", x$lambda2.1se, "\n\n")
  printCoefmat2(x$tab2, exclude = x$exclude2, digits = x$digits)
  cat("\n")
}


#' @exportS3Method coef pCure
coef.pCure <- function(object, part = "both", ...) {
  out <- list(object$fit1$b, object$fit2$b)
  if (object$control$model == "mixture") {
    part <- match.arg(part, c("both", "incidence", "latency"))
    names(out) <- c("incidence", "latency")
  }
  if (object$control$model == "promotion") {
    part <- match.arg(part, c("both", "long", "short"))
    names(out) <- c("long", "short")
  }
  if (part == "both") return(out)
  else return(out[[part]])
}

#' @exportS3Method vcov pCure
vcov.pCure <- function(object, part = "both", ...) {
  out <- list(object$fit1$vb, object$fit2$vb)
  colnames(out[[1]]) <- rownames(out[[1]]) <- names(object$fit1$b)
  colnames(out[[2]]) <- rownames(out[[2]]) <- names(object$fit2$b)
  if (object$control$model == "mixture") {
    part <- match.arg(part, c("both", "incidence", "latency"))
    names(out) <- c("incidence", "latency")
  }
  if (object$control$model == "promotion") {
    part <- match.arg(part, c("both", "long", "short"))
    names(out) <- c("long", "short")
  }
  if (part == "both") return(out)
  else return(out[[part]])
}


#' @exportS3Method residuals pCure
residuals.pCure <- function(object, part = "both", ...) {
  out <- list(object$fit1$resid, object$fit2$resid)
  if (object$control$model == "mixture") {
    part <- match.arg(part, c("both", "incidence", "latency"))
    names(out) <- c("incidence", "latency")
  }
  if (object$control$model == "promotion") {
    part <- match.arg(part, c("both", "long", "short"))
    names(out) <- c("long", "short")
  }
  if (part == "both") return(out)
  else return(out[[part]])
}

#' @exportS3Method fitted pCure
fitted.pCure <- function(object, part = "both", ...) {
  out <- list(object$fit1$fitted, object$fit2$fitted)
  if (object$control$model == "mixture") {
    part <- match.arg(part, c("both", "incidence", "latency"))
    names(out) <- c("incidence", "latency")
  }
  if (object$control$model == "promotion") {
    part <- match.arg(part, c("both", "long", "short"))
    names(out) <- c("long", "short")
  }
  if (part == "both") return(out)
  else return(out[[part]])
}

#' @import utils 
utils::globalVariables(c("vari", "lam"))

#' Plot method for 'pCure' objects 
#'
#' @param x An object of class 'pCure', usually returned by the 'pCure()' function.
#' @param part A character string specifies which component of the cure model to plot.
#' The default is "both", which plots both the incidence and latency components if a
#' mixture cure model was fitted,
#' or both the long- and short-term effects if a promotion time model was fitted.
#' @param type A character string specifying the type of plot to generate.
#' Available options are "residuals," "cv," and "trace,"
#' which correspond to the pseudo-residual plot, cross-validation plot,
#' and trace plot for different values of the tuning parameter, respectively.
#' @param ... Other arguments for future extension.
#' 
#' @importFrom ggplot2 ggplot geom_point geom_line  geom_errorbar ggtitle scale_x_reverse geom_vline
#' @importFrom ggplot2 xlab ylab aes aes_string
#' @importFrom rlang call_args
#' @importFrom ggpubr ggarrange
#'
#' @exportS3Method plot pCure
plot.pCure <- function(x, part = "both", type = c("residuals", "cv", "trace"),...) {
  if (!is.pCure(x)) stop("Must be a pCure object")
  type <- match.arg(type, c("residuals", "cv", "trace"))
  if (type == "residuals") {
    dat <- data.frame(fitted = unlist(fitted(x)),
                      resid = unlist(resid(x)))
    tmp <- lapply(fitted(x), length)
    dat$Comp <- rep(names(tmp), unlist(tmp))
    if (x$control$model == "mixture") 
      part <- match.arg(part, c("both", "incidence", "latency"))
    if (x$control$model == "promotion") 
      part <- match.arg(part, c("both", "long", "short"))
    if (part != "both") dat <- dat[dat$Comp == part,] # subset(dat, Comp == part)
    if (length(unique(dat$Comp)) == 1)
      p <- ggplot(dat, aes(x = fitted,y = resid))
    else 
      p <- ggplot(dat, aes_string(x = "fitted", y = "resid", color = "Comp"))
    p <- p + geom_point() + xlab("Fitted values") + ylab("Residuals")
  }
  else {
    if (is.null(x$control$lambda1) & is.null(x$control$lambda2)) 
      stop("No tuning parameters for penalization have been detected.")
  }
  if (type == "cv") {    
    d1 <- data.frame(lam = x$control$lambda1,
                     mean = colMeans(x$fit1$cv.raw),
                     sd = apply(x$fit1$cv.raw, 2, sd))
    d2 <- data.frame(lam = x$control$lambda2,
                     mean = colMeans(x$fit2$cv.raw),
                     sd = apply(x$fit2$cv.raw, 2, sd))
    if (nrow(d1) > 0) {
      if (call_args(x$call)$lambda1 == "auto") {
        d1$lam <- log(d1$lam)
        xlab <- expression(log(lambda))
      } else
        xlab <- expression(lambda)
      title <- ifelse(x$control$model == "mixture", "Incidence component", "Long-term effect")
      p1 <- ggplot(d1, aes(x = lam, y = mean)) +
        geom_point() + 
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .002) +
        geom_point(data = d1[which.min(d1$mean),], color = "red") +
        geom_errorbar(data = d1[which.min(d1$mean),], aes(ymin = mean - sd, ymax = mean + sd), width = .002, color = "red") +
        ylab("Prediction error") + ggtitle(title) + xlab(xlab)
    }
    if (nrow(d2) > 0) {
      if (call_args(x$call)$lambda2 == "auto") {
        d2$lam <- log(d2$lam)
        xlab <- expression(log(lambda))
      } else
        xlab <- expression(lambda)
      title <- ifelse(x$control$model == "mixture", "Latency component", "Short-term effect")      
      p2 <- ggplot(d2, aes(x = lam, y = mean)) +
        geom_point() + 
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .002) +
        geom_point(data = d2[which.min(d2$mean),], color = "red") +
        geom_errorbar(data = d2[which.min(d2$mean),], aes(ymin = mean - sd, ymax = mean + sd), width = .002, color = "red") +
        ylab("Prediction error") + ggtitle(title) + xlab(xlab)
    }
    p <- ggarrange(p1, p2, ncol = 1)    
  }
  if (type == "trace") {
    coef1 <- sapply(x$control$lambda1, function(i) update(x, lambda1 = i, formula2 = NULL, lambda2 = NULL)$fit1$b)
    coef2 <- sapply(x$control$lambda2, function(i) update(x, lambda2 = i, formula1 = NULL, lambda1 = NULL)$fit2$b)
    coef1 <- coef1[rownames(coef1) != "(Intercept)",]
    coef2 <- coef2[rownames(coef2) %in% rownames(coef1),]
    d1 <- data.frame(coef = c(coef1),
                     vari = rownames(coef1),
                     lam = rep(x$control$lambda1, each = nrow(coef1)))
    d2 <- data.frame(coef = c(coef2),
                     vari = rownames(coef2), 
                     lam = rep(x$control$lambda2, each = nrow(coef2)))
    if (nrow(d1) > 0) {
      coef10 <- t(ifelse(coef1 < 1e-3, 0, coef1))
      if (call_args(x$call)$lambda1 == "auto") {
        keep1 <- !duplicated(coef10)
        d1 <- subset(d1, lam %in% x$control$lambda1[keep1])
      }
      xint <- unique(d1$lam[d1$lam == x$fit1$lambda1.min])
      if (call_args(x$call)$lambda1 == "auto") {
        d1$lam <- log(d1$lam)
        xlab <- expression(log(lambda))
        xint <- log(xint)
      } else
        xlab <- expression(lambda)
      title <- ifelse(x$control$model == "mixture", "Incidence component", "Long-term effect")
      p1 <- ggplot(d1, aes(x = lam, y = coef, color = vari)) +
        geom_line() + scale_x_reverse() +
        geom_vline(xintercept = xint, linetype="dotted") + 
        ylab("Coefficients") + ggtitle(title) + xlab(xlab)
    }
    if (nrow(d2) > 0) {
      coef20 <- t(ifelse(coef2 < 1e-3, 0, coef2))
      if (call_args(x$call)$lambda2 == "auto") {
        keep2 <- !duplicated(coef20)
        d2 <- subset(d2, lam %in% x$control$lambda2[keep2])
      }
      xint <- unique(d2$lam[d2$lam == x$fit2$lambda2.min])
      if (call_args(x$call)$lambda2 == "auto") {
        d2$lam <- log(d2$lam)
        xlab <- expression(log(lambda))
        xint <- log(xint)
      } else
        xlab <- expression(lambda)
      title <- ifelse(x$control$model == "mixture", "Latency component", "Short-term effect")      
      p2 <- ggplot(d2, aes(x = lam, y = coef, color = vari)) +
        geom_line() + scale_x_reverse() +
        geom_vline(xintercept = xint, linetype="dotted") + 
        ylab("Coefficients") + ggtitle(title) + xlab(xlab)
    }
    p <- ggarrange(p1, p2, ncol = 1)    
  }
  p
}




