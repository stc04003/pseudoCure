fitPHMC <- function(X1, X2, time, status, t0, control) {
    tmax <- max(time[status > 0])
    KMs <- pseudoKM1(time, status, c(t0, tmax))
    n <- length(time)
    S <- KMs[1:length(t0), n + 1]
    Si <- KMs[1:length(t0), 1:n]
    cure <- KMs[length(t0) + 1, n + 1]
    curei <- KMs[length(t0) + 1, 1:n]
    if (is.null(control$binit1)) control$binit1 <- rep(0, ncol(X1))
    fit1 <- gee(curei, X1, control$binit1, rep(1, n), "logit", "ind", control$tol, control$maxit)
    Si <- n * (S - cure) / (1 - cure) - c((n - 1) * (Si - curei) / (1 - curei))
    Fi <- 1 - Si
    nt <- rep(length(t0), length(time))
    X2 <- cbind(model.matrix(~ 0 + as.factor(rep(t0, length(time)))), X2)
    if (is.null(control$binit2)) control$binit2 <- rep(0, ncol(X2))
    fit2 <- gee(Fi, X2, control$binit2, nt, "cloglog",
                control$corstr, control$tol, control$maxit)
    list(fit1 = fit1, fit2 = fit2)
}

