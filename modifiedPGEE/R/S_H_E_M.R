S_H_E_M <- function(N,nt,y,X,nx,family,beta_new,Rhat,fihat,lambda,pindex,eps=10^-6) {
    aindex <- cumsum(nt)
    index <- c(0, aindex[-length(aindex)])    
    eta <- X %*% beta_new
    mu <- family$linkinv(eta)
    muEta <- family$mu.eta(eta)
    ## This is E on Wang et al.(2012)
    E1 <- diag(q_scad(abs(c(beta_new)), lambda) / (abs(as.vector(beta_new)) + eps))
    
    if (is.null(pindex) == TRUE) {
        E <- E1
    } else
        if (is.null(pindex) != TRUE) {
            E1[, pindex] <- 0
            E <- E1
        }
    ##  FamilyVar <- family$variance(mu)
    ##  FamilyMu <- family$mu.eta(eta)    
    ##  sum201 <- matrix(0, nx, 1)      #gradient:S(sum201)
    ##  sum301 <- matrix(0, nx, nx)     #naive variance:H(sum301)
    ##  sum401 <- matrix(0, nx, nx)     #a component for robust variance:M(sum401)
    InnerMatrix <- SHM_CppArma(y, mu, X, Rhat, N, nx, nt, index, muEta)
    S <- fihat * InnerMatrix$sumS
    H <- fihat * InnerMatrix$sumH
    E <- E
    M <- fihat * InnerMatrix$sumM  
    return(list(
        "S" = S,
        "H" = H,
        "E" = E,
        "M" = M
    ))
}


