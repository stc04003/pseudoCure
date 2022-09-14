fitPHMC <- function(X1, X2, time, status, t0, control) {
    tmax <- max(time[status > 0])
    KMs <- pseudoKM1(time, status, c(t0, tmax))
    n <- length(time)
    S <- KMs[1:length(t0), n + 1]
    Si <- KMs[1:length(t0), 1:n]
    cure <- KMs[length(t0) + 1, n + 1]
    curei <- KMs[length(t0) + 1, 1:n]
    thetai <- n * (1 - cure) - (n - 1) * (1 - curei)
    if (is.null(control$binit1)) control$binit1 <- rep(0, ncol(X1))
    if (is.null(control$exclude1)) control$exclude1 <- rep(0, ncol(X1))
    if (is.null(control$lambda1))
        fit1 <- gee(thetai, X1, control$binit1, rep(1, n), "logit", "ind", control$tol, control$maxit)
    else if (length(control$lambda1) == 1)
        fit1 <- pgee(thetai, X1, control$binit1, rep(1, n), control$exclude1,
                     "logit", control$penalty1, "ind",
                     control$lambda1, control$eps, control$tol, control$maxit)
    else {
        cv.raw1 <- pgeeCV(thetai, X1, control$binit1, rep(1, n), control$exclude1,
                          "logit", control$penalty1, "ind", control$nfolds, 
                          control$lambda1, control$eps, control$tol, control$maxit)
        cv.mean <- colMeans(cv.raw1)
        cv.sd <- apply(cv.raw1, 2, sd) / control$nfolds
        cv.1se <- cv.mean + cv.sd
        lambda1.min <- control$lambda1[which.min(cv.mean)]
        lambda1.1se <- control$lambda1[which.max(which(cv.mean < cv.1se[which.min(cv.mean)]))]
        fit1 <- pgee(thetai, X1, control$binit1, rep(1, n), control$exclude1,
                     "logit", control$penalty1, "ind",
                     lambda1.1se, control$eps, control$tol, control$maxit)
        fit1$cv.raw <- cv.raw1
        fit1$lambda1.min <- lambda1.min
        fit1$lambda1.1se <- lambda1.1se
    }    
    Si <- n * (S - cure) / (1 - cure) - c((n - 1) * t((t(Si) - curei) / (1 - curei)))
    Fi <- 1 - Si
    nt <- rep(length(t0), length(time))
    X22 <- cbind(model.matrix(~ 0 + as.factor(rep(t0, length(time)))),
                 X2[rep(1:nrow(X2), each = length(t0)),])
    colnames(X22) <- c(paste("t", seq_along(t0), sep = ""), colnames(X2))
    if (is.null(control$binit2)) control$binit2 <- rep(0, ncol(X22))
    if (is.null(control$exclude2)) control$exclude2 <- rep(0, ncol(X22))
    if (is.null(control$lambda2))
        fit2 <- gee(Fi, X22, control$binit2, nt, "cloglog",
                    control$corstr, control$tol, control$maxit)
    else if (length(control$lambda2) == 1)
        fit2 <- pgee(Fi, X22, control$binit2, nt, control$exclude2, "cloglog",
                     control$penalty2, control$corstr, control$lambda2,
                     control$eps, control$tol, control$maxit)
    else {
        cv.raw2 <- pgeeCV(Fi, X22, control$binit2, nt, control$exclude2,
                          "cloglog", control$penalty2, control$corstr, control$nfolds, 
                          control$lambda2, control$eps, control$tol, control$maxit)
        cv.mean <- colMeans(cv.raw2)
        cv.sd <- apply(cv.raw2, 2, sd) / control$nfolds
        cv.1se <- cv.mean + cv.sd
        lambda2.min <- control$lambda2[which.min(cv.mean)]
        lambda2.1se <- control$lambda2[which.max(which(cv.mean < cv.1se[which.min(cv.mean)]))]
        fit2 <- pgee(Fi, X22, control$binit2, nt, control$exclude2,
                     "cloglog", control$penalty2, control$corstr,
                     lambda2.1se, control$eps, control$tol, control$maxit)
        fit2$cv.raw <- cv.raw2
        fit2$lambda1.min <- lambda2.min
        fit2$lambda1.1se <- lambda2.1se
    }
    fit1$b <- drop(fit1$b)
    fit2$b <- drop(fit2$b)
    names(fit1$b) <- colnames(X1)
    names(fit2$b) <- colnames(X22)
    fit1$vb <- with(fit1, ginv(H + n * E) %*% M %*% ginv(H + n * E))
    fit2$vb <- with(fit2, ginv(H + n * E) %*% M %*% ginv(H + n * E))
    list(fit1 = fit1, fit2 = fit2)
}

