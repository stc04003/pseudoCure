#' @importFrom MASS ginv
fitPHPH <- function(X1, X2, time, status, t0, control) {
    tmax <- max(time[status > 0])
    KMs <- pseudoKM1(time, status, c(t0, tmax))
    n <- length(time)
    S <- KMs[1:length(t0), n + 1]
    Si <- KMs[1:length(t0), 1:n]
    cure <- KMs[length(t0) + 1, n + 1]
    curei <- KMs[length(t0) + 1, 1:n]
    thetai <- n * (-log(cure)) - (n - 1) * (-log(curei))
    if (is.null(control$binit1)) control$binit1 <- rep(0, ncol(X1))
    fit1 <- gee(thetai, X1, control$binit1, rep(1, n), "log", "ind", control$tol, control$maxit)    
    Fi <- n * (log(S) / log(cure)) - c((n - 1) * t(log(t(Si)) / log(curei)))
    nt <- rep(length(t0), length(time))
    X22 <- cbind(model.matrix(~ 0 + as.factor(rep(t0, length(time)))),
                 X2[rep(1:nrow(X2), each = length(t0)),])
    colnames(X22) <- c(paste("t", seq_along(t0), sep = ""), colnames(X2))
    if (is.null(control$binit2)) control$binit2 <- rep(0, ncol(X22))
    fit2 <- gee(Fi, X22, control$binit2, nt, "cloglog",
                control$corstr, control$tol, control$maxit)
    fit1$b <- drop(fit1$b)
    fit2$b <- drop(fit2$b)
    names(fit1$b) <- colnames(X1)
    names(fit2$b) <- colnames(X22)
    fit1$vb <- with(fit1, ginv(H) %*% M %*% ginv(H))
    fit2$vb <- with(fit2, ginv(H) %*% M %*% ginv(H))
    list(fit1 = fit1, fit2 = fit2)
}

