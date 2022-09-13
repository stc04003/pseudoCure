fitPHMC <- function(X1, X2, time, status, t0,
                    lambda1, penalty1, exclude1, lambda2, penalty2, exclude2, 
                    control) {
    tmax <- max(time[status > 0])
    KMs <- pseudoKM1(time, status, c(t0, tmax))
    n <- length(time)
    S <- KMs[1:length(t0), n + 1]
    Si <- KMs[1:length(t0), 1:n]
    cure <- KMs[length(t0) + 1, n + 1]
    curei <- KMs[length(t0) + 1, 1:n]
    thetai <- n * (1 - cure) - (n - 1) * (1 - curei)
    if (is.null(control$binit1)) control$binit1 <- rep(0, ncol(X1))
    if (is.null(lambda1))
        fit1 <- gee(thetai, X1, control$binit1, rep(1, n), "logit", "ind", control$tol, control$maxit)
    else
        fit1 <- pgee(thetai, X1, control$binit1, rep(1, n), exclude1,
                     "logit", penalty1, "ind", lambda1, control$eps, control$tol, control$maxit)
    Si <- n * (S - cure) / (1 - cure) - c((n - 1) * t((t(Si) - curei) / (1 - curei)))
    Fi <- 1 - Si
    nt <- rep(length(t0), length(time))
    X22 <- cbind(model.matrix(~ 0 + as.factor(rep(t0, length(time)))),
                 X2[rep(1:nrow(X2), each = length(t0)),])
    colnames(X22) <- c(paste("t", seq_along(t0), sep = ""), colnames(X2))
    if (is.null(control$binit2)) control$binit2 <- rep(0, ncol(X22))
    if (is.null(lambda2))
        fit2 <- gee(Fi, X22, control$binit2, nt, "cloglog",
                    control$corstr, control$tol, control$maxit)
    else
        fit2 <- pgee(Fi, X22, control$binit2, nt, exclude2, "cloglog",
                     penalty2, control$corstr, lambda2, control$eps, control$tol, control$maxit)
    fit1$b <- drop(fit1$b)
    fit2$b <- drop(fit2$b)
    names(fit1$b) <- colnames(X1)
    names(fit2$b) <- colnames(X22)
    fit1$vb <- with(fit1, ginv(H + n * E) %*% M %*% ginv(H + n * E))
    fit2$vb <- with(fit2, ginv(H + n * E) %*% M %*% ginv(H + n * E))
    list(fit1 = fit1, fit2 = fit2)
}

