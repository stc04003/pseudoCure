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
      cat("\nFitted with the mixture cure rate model")
    if (x$control$model == "promotion")
      cat("\nFitted with the promotion time model")
  } else {
    if (x$control$model == "mixture")
      cat("\nFitted with the penalized mixture cure rate model")
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
  out$lambda1.min <- object$fit1$lambda1.min
  out$lambda2.min <- object$fit2$lambda2.min
  out$lambda1.1se <- object$fit1$lambda1.1se
  out$lambda2.1se <- object$fit2$lambda2.1se
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
      cat("\nFitted with the mixture cure rate model")
    if (x$control$model == "promotion")
      cat("\nFitted with the promotion time model")
  } else {
    if (x$control$model == "mixture")
      cat("\nFitted with the penalized mixture cure rate model")
    if (x$control$model == "promotion")
      cat("\nFitted with the penalized promotion time model")
  }
  if (x$control$model == "mixture") cat("\nIncidence component:\n")
  else cat("\nLong-term effect:\n")
  if (!is.null(x$lambda1.min))
    cat("Tuning parameters selected by minimum prediction error:", x$lambda1.min, "\n\n")
  ## cat("\nTuning parameters selected by 1se rule: ", x$lambda1.1se, "\n\n")
  printCoefmat2(x$tab1, exclude = x$exclude1)
  if (x$control$model == "mixture") cat("\nLatency component:\n")
  else cat("\nShort-term effect:\n")
  if (!is.null(x$lambda2.min))
    cat("Tuning parameters selected by minimum prediction error:", x$lambda2.min, "\n\n")
  ## cat("\nTuning parameters selected by 1se rule: ", x$lambda2.1se, "\n\n")
  printCoefmat2(x$tab2, exclude = x$exclude2)
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


#' @importFrom ggplot2 ggplot geom_point geom_line  geom_errorbar ggtitle scale_x_reverse xlab ylab aes geom_vline
#' @importFrom rlang call_args
#' @importFrom ggpubr ggarrange
#' 
#' @exportS3Method plot pCure
plot.pCure <- function(object, part = "both", type = c("residuals", "cv", "trace"),...) {
  type <- match.arg(type, c("residuals", "cv", "trace"))
  if (type == "residuals") {
    dat <- data.frame(fitted = unlist(fitted(object)),
                      resid = unlist(resid(object)))
    tmp <- lapply(fitted(object), length)
    dat$Component <- rep(names(tmp), unlist(tmp))
    if (object$control$model == "mixture") 
      part <- match.arg(part, c("both", "incidence", "latency"))
    if (object$control$model == "promotion") 
      part <- match.arg(part, c("both", "long", "short"))
    if (part != "both") dat <- subset(dat, Component == part)
    if (length(unique(dat$Component)) == 1)
      p <- ggplot(dat, aes(x = fitted,y = resid))
    else 
      p <- ggplot(dat, aes(x = fitted,y = resid, color = Component))
    p <- p + geom_point() + xlab("Fitted values") + ylab("Residuals")
  }
  else {
    if (is.null(object$control$lambda1) & is.null(object$control$lambda2)) 
      stop("No tuning parameters for penalization have been detected.")
  }
  if (type == "cv") {    
    d1 <- data.frame(lambda = object$control$lambda1,
                     mean = colMeans(object$fit1$cv.raw),
                     sd = apply(object$fit1$cv.raw, 2, sd))
    d2 <- data.frame(lambda = object$control$lambda2,
                     mean = colMeans(object$fit2$cv.raw),
                     sd = apply(object$fit2$cv.raw, 2, sd))
    if (nrow(d1) > 0) {
      if (call_args(object$call)$lambda1 == "auto") {
        d1$lambda <- log(d1$lambda)
        xlab <- expression(log(lambda))
      } else
        xlab <- expression(lambda)
      title <- ifelse(object$control$model == "mixture", "Incidence component", "Long-term effect")
      p1 <- ggplot(d1, aes(x = lambda, y = mean)) +
        geom_point() + 
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .002) +
        geom_point(data = d1[which.min(d1$mean),], color = "red") +
        geom_errorbar(data = d1[which.min(d1$mean),], aes(ymin = mean - sd, ymax = mean + sd), width = .002, color = "red") +
        ylab("Prediction error") + ggtitle(title) + xlab(xlab)
    }

    if (nrow(d2) > 0) {
      if (call_args(object$call)$lambda2 == "auto") {
        d2$lambda <- log(d2$lambda)
        xlab <- expression(log(lambda))
      } else
        xlab <- expression(lambda)
      title <- ifelse(object$control$model == "mixture", "Latency component", "Short-term effect")      
      p2 <- ggplot(d2, aes(x = lambda, y = mean)) +
        geom_point() + 
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .002) +
        geom_point(data = d2[which.min(d2$mean),], color = "red") +
        geom_errorbar(data = d2[which.min(d2$mean),], aes(ymin = mean - sd, ymax = mean + sd), width = .002, color = "red") +
        ylab("Prediction error") + ggtitle(title) + xlab(xlab)
    }
    p <- ggarrange(p1, p2, ncol = 1)    
  }
  if (type == "trace") {
    coef1 <- sapply(object$control$lambda1, function(i) update(object, lambda1 = i, formula2 = NULL, lambda2 = NULL)$fit1$b)
    coef2 <- sapply(object$control$lambda2, function(i) update(object, lambda2 = i, formula1 = NULL, lambda1 = NULL)$fit2$b)
    coef1 <- coef1[rownames(coef1) != "(Intercept)",]
    coef2 <- coef2[rownames(coef2) %in% rownames(coef1),]
    d1 <- data.frame(coef = c(coef1),
                     variable = rownames(coef1),
                     lambda = rep(object$control$lambda1, each = nrow(coef1)))
    d2 <- data.frame(coef = c(coef2),
                     variable = rownames(coef2), 
                     lambda = rep(object$control$lambda2, each = nrow(coef2)))
    if (nrow(d1) > 0) {
      coef10 <- t(ifelse(coef1 < 1e-3, 0, coef1))
      if (call_args(object$call)$lambda1 == "auto") {
        keep1 <- !duplicated(coef10)
        d1 <- subset(d1, lambda %in% object$control$lambda1[keep1])
      }
      xint <- unique(d1$lambda[d1$lambda == object$fit1$lambda1.min])
      if (call_args(object$call)$lambda1 == "auto") {
        d1$lambda <- log(d1$lambda)
        xlab <- expression(log(lambda))
        xint <- log(xint)
      } else
        xlab <- expression(lambda)
      title <- ifelse(object$control$model == "mixture", "Incidence component", "Long-term effect")
      p1 <- ggplot(d1, aes(x = lambda, y = coef, color = variable)) +
        geom_line() + scale_x_reverse() +
        geom_vline(xintercept = xint, linetype="dotted") + 
        ylab("Coefficients") + ggtitle(title) + xlab(xlab)
    }
    if (nrow(d2) > 0) {
      coef20 <- t(ifelse(coef2 < 1e-3, 0, coef2))
      if (call_args(object$call)$lambda2 == "auto") {
        keep2 <- !duplicated(coef20)
        d2 <- subset(d2, lambda %in% object$control$lambda2[keep2])
      }
      xint <- unique(d2$lambda[d2$lambda == object$fit2$lambda2.min])
      if (call_args(object$call)$lambda2 == "auto") {
        d2$lambda <- log(d2$lambda)
        xlab <- expression(log(lambda))
        xint <- log(xint)
      } else
        xlab <- expression(lambda)
      title <- ifelse(object$control$model == "mixture", "Latency component", "Short-term effect")      
      p2 <- ggplot(d2, aes(x = lambda, y = coef, color = variable)) +
        geom_line() + scale_x_reverse() +
        geom_vline(xintercept = xint, linetype="dotted") + 
        ylab("Coefficients") + ggtitle(title) + xlab(xlab)
    }
    p <- ggarrange(p1, p2, ncol = 1)    
  }
  p
}




