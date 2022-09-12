CVfit <- function(formula, id, data, family, scale.fix, scale.value, fold, 
                  lambda.vec, pindex, eps, maxiter, tol) {
  call_cv <- match.call()
  mf <- match.call(expand.dots=FALSE)
  mf$family <- mf$link <- mf$varfun <-
    mf$scale.fix <- mf$scale.value <-
      mf$fold <- mf$lambda.vec <-mf$pindex <-
        mf$eps <- mf$maxiter <- mf$tol <- NULL
  if(is.null(mf$id)) mf$id <- as.name("id")
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  ## X <- model.matrix(Terms, mf, contrasts)
  X <- model.matrix(Terms, mf)
  id <- model.extract(mf, id)
  if(missing(family)) family=gaussian(link="identity")
  ## if (is.character(family)) family <- get(family)
  ## if (is.function(family))  family <- family()
  if(missing(pindex)) pindex <- NULL
  if(missing(scale.fix))   scale.fix <- TRUE
  if(missing(scale.value)) scale.value <- 1
  if(missing(eps)) eps <- 10^-6
  if(missing(maxiter)) maxiter <- 30
  if(missing(tol)) tol <- 10^-3
  N <- length(unique(id))
  K <- dim(X)[2]-1
  nx <- dim(X)[2]
  nt <- dim(X)[1]/N
  nt <- rep(nt,N)
  ## pay attention to this part  
  lam.min <- -1
  cv.min <- Inf
  cv.vect <- list()
  for (j in 1:length(lambda.vec)) {
    ## get one of the lambda's
    lam.temp <- lambda.vec[j]
    ## initial value for cv.values is 0.
    cv.value <- 0
    tmp <- NULL
    for (k in 1:fold) {
      ## select the index that will be omitted.
      index.cv <- ((k-1)*nt[1]*(N/fold)+1):(k*nt[1]*(N/fold))
      ## training part#
      y.train <- y[-index.cv]
      if (colnames(X)[1]== "(Intercept)")
        x.train <- X[-index.cv,-1]
      else
        x.train <- X[-index.cv,]
      id.train <- id[-index.cv]
      ## compute beta.train
      data.train=data.frame("id"=id.train,"y"=y.train, x.train)
      mm<-match.call(expand.dots = FALSE)
      mm$formula<-formula
      mm$id<-id.train
      mm$data<-data.train
      mm$na.action<-"na.omit"
      mm$family<-family
      mm$corstr<-"independence"
      mm$Mv<-NULL
      ## mm$beta_int <- rep(0,nx)  ##NULL##
      mm$beta_int <- NULL
      mm$R<-NULL
      mm$scale.fix<-scale.fix
      mm$scale.value<-scale.value
      mm$lambda<-lam.temp
      mm$pindex<-pindex
      mm$eps<-eps
      mm$maxiter<-maxiter
      mm$tol<-tol
      mm$silent<-TRUE
      mm$lambda.vec<-NULL
      mm$fold<-NULL
      mm[[1]]<-as.name("PGEE")
      beta.train <- tryCatch(
        suppressWarnings(eval(mm, parent.frame())$coefficients),
        error = function(e) {
          mm$beta_int <- geese(as.formula(formula), data = data.train, 
                               id = id, family = family$family, mean.link = family$link)$beta
          eval(mm, parent.frame())$coefficients
        })
      y.cv<-y[index.cv]
      x.cv<-X[index.cv,]
      id.cv<-id[index.cv]
      yy=y.cv
      eta=x.cv%*%beta.train 
      mu=family$linkinv(eta)
      ##family$dev.resids gives the square of the residuals
      ## cv.value<-cv.value+sum((family$dev.resids(yy,mu,wt=1)))
      tmp <- cbind(tmp, family$dev.resids(yy,mu,wt=1))
    } #k
    cv.vect[[j]] <- tmp
    ## cv.vect<-c(cv.vect, cv.value)
    ## if(cv.value<cv.min) {
    ##     lam.min<-lam.temp
    ##     cv.min<-cv.value
    ## }      
  } #j
  ## cvsd <- unlist(lapply(cv.vect, function(x)
  ##     sqrt(mean(apply(x, 2, function(y) mean((y - mean(y))^2))))))
  cvsd <- unlist(lapply(cv.vect, function(x) sd(colMeans(x) / fold)))
  cvm <- unlist(lapply(cv.vect, mean))
  ind <- which.max(which(cvm < (cvm + cvsd)[which.min(cvm)]))
  out<-list()
  attr(out, "class") <- c("CVfit")
  out$fold=fold
  out$lam.vect <- lambda.vec
  out$cv.vect <- cvm
  out$cv.raw <- cv.vect
  out$cvsd <- cvsd
  out$lam.opt <- lambda.vec[which.min(cvm)]
  out$cv.min <- cvm[which.min(cvm)]
  out$lam.opt.1se <- lambda.vec[ind]
  out$cv.min.1se <- cvm[ind] ## cvm[which.min(cvm)]
  out$call <- call_cv
  out   
}
