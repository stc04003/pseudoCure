#' @importFrom MASS ginv
fitPHPH <- function(X1, X2, time, status, t0, control) {
    tmax <- max(time[status > 0])
    KMs <- pseudoKM1(time, status, c(t0, tmax))
    n <- length(time)
    S <- KMs[1:length(t0), n + 1]
    Si <- KMs[1:length(t0), 1:n]
    cure <- KMs[length(t0) + 1, n + 1]
    curei <- KMs[length(t0) + 1, 1:n]
    fit1 <- fit2 <- NULL
    if (!is.null(X1)) fit1 <- fitPHPH1(X1, time, status, control, cure, curei)
    if (!is.null(X2)) fit2 <- fitPHPH2(X2, time, status, t0, control, cure, curei, S, Si)    
    list(fit1 = fit1, fit2 = fit2)
}

fitPHPH1 <- function(X1, time, status, control,
                     cure = NULL, curei = NULL) {
  n <- length(time)
  if (any(is.null(cure), is.null(curei))) {
    KMs <- pseudoKM(time, status)
    cure <- KMs[n + 1]
    curei <- KMs[1:n]
  }
  thetai <- n * (-log(cure)) - (n - 1) * (-log(curei))
  if (is.null(control$binit1)) control$binit1 <- rep(0, ncol(X1))
  if (is.null(control$exclude1)) control$exclude1 <- rep(0, ncol(X1))
  if (is.null(control$lambda1)) 
    fit1 <- gee(thetai, X1, control$binit1, rep(1, n), "log", "ind", control$tol, control$maxit)
  else if (length(control$lambda1) == 1)
    fit1 <- pgee(thetai, X1, control$binit1, rep(1, n), control$exclude1,
                 "log", control$penalty1, "ind", control$lambda1,
                 control$eps, control$tol, control$maxit)
  else {
    cv.raw1 <- pgeeCV(thetai, X1, control$binit1, rep(1, n), control$exclude1,
                      "log", control$penalty1, "ind", control$nfolds, 
                      control$lambda1, control$eps, control$tol, control$maxit)
    cv.mean <- colMeans(cv.raw1)
    cv.sd <- apply(cv.raw1, 2, sd) / control$nfolds
    cv.1se <- cv.mean + cv.sd
    lambda1.min <- control$lambda1[which.min(cv.mean)]
    lambda1.1se <- control$lambda1[which.max(which(cv.mean < cv.1se[which.min(cv.mean)]))]
    fit1 <- pgee(thetai, X1, control$binit1, rep(1, n), control$exclude1,
                 "log", control$penalty1, "ind",
                 lambda1.1se, control$eps, control$tol, control$maxit)        
    fit1$cv.raw <- cv.raw1
    fit1$lambda1.min <- lambda1.min
    fit1$lambda1.1se <- lambda1.1se
  }
  fit1$b <- drop(fit1$b)
  names(fit1$b) <- colnames(X1)
  fit1$vb <- with(fit1, ginv(H + n * E) %*% M %*% ginv(H + n * E))
  fit1$fitted <- drop(exp(X1 %*% fit1$b))
  fit1$resid <- drop(thetai - fit1$fitted)
  return(fit1)
}

fitPHPH2 <- function(X2, time, status, t0, control,
                     cure = NULL, curei = NULL, S = NULL, Si = NULL) {
  n <- length(time)
  tmax <- max(time[status > 0])
  if (any(is.null(cure), is.null(curei), is.null(S), is.null(Si))) {
    KMs <- pseudoKM1(time, status, c(t0, tmax))
    S <- KMs[1:length(t0), n + 1]
    Si <- KMs[1:length(t0), 1:n]
    cure <- KMs[length(t0) + 1, n + 1]
    curei <- KMs[length(t0) + 1, 1:n]         
  }
  Fi <- n * (log(S) / log(cure)) - c((n - 1) * t(log(t(Si)) / log(curei)))
  nt <- rep(length(t0), length(time))
  X22 <- cbind(model.matrix(~ 0 + as.factor(rep(t0, length(time)))),
               X2[rep(1:nrow(X2), each = length(t0)),])
  colnames(X22) <- c(paste("t", seq_along(t0), sep = ""), colnames(X2))
  if (is.null(control$binit2)) control$binit2 <- runif(ncol(X22), -.001, .001)
  ## control$binit2 <- rep(0, ncol(X22))
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
    fit2$lambda2.min <- lambda2.min
    fit2$lambda2.1se <- lambda2.1se        
  }
  fit2$b <- drop(fit2$b)
  names(fit2$b) <- colnames(X22)
  fit2$vb <- with(fit2, ginv(H + n * E) %*% M %*% ginv(H + n * E))
  fit2$fitted <- drop(1 - exp(-exp(X22 %*% fit2$b)))
  fit2$resid <- -drop(Fi - fit2$fitted)
  return(fit2)
}
