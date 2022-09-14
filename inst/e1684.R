## #######################################################################
## Load and prepare the Melanoma data 
## #######################################################################
library(survival)
library(microbenchmark)
library(Rcpp)
library(pseudoCure)

data(e1684, package = "smcure")
e1684 <- na.omit(e1684)
e1684$SEX2 <- as.factor(e1684$SEX)
e1684$TRT2 <- as.factor(e1684$TRT)

## #######################################################################
## EM-based approach for mixture cure model
## #######################################################################
library(smcure)
smcure(Surv(FAILTIME, FAILCENS) ~ AGE + TRT + SEX, 
       ~ AGE + TRT + SEX,
       data = e1684, model = "ph", nboot = 100)

e1684$TRT.SEX <- e1684$TRT * e1684$SEX
coxPHMC <- smcure(Surv(FAILTIME, FAILCENS) ~ AGE + TRT + SEX + TRT.SEX, 
                  ~ AGE + TRT + SEX + TRT.SEX, 
                  data = e1684, model = "ph", nboot = 100)

newDat <- expand.grid(TRT = 0:1, SEX = 0:1)
newDat$AGE <- with(e1684, tapply(AGE, interaction(TRT, SEX), mean))
newDat$cure.EM <- 1 - 1 / (1 + exp(-model.matrix(~ AGE + TRT * SEX, data = newDat) %*% coxPHMC$b)) 
newDat


## #######################################################################
## Pseudo-observation approach for mixture cure model
## #######################################################################

head(e1684)

foo <- pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684)
foo1 <- pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684, model = 'p')


str(pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684))
str(pCure(~ AGE + TRT * SEX2, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684))


## Incidence component
n <- nrow(e1684)
tmax <- max(e1684$FAILTIME[e1684$FAILCENS > 0])
KM <- survfit(Surv(FAILTIME, FAILCENS) ~ 1, data = e1684)
(cure <- min(KM$surv) )

e1684.inc <- e1684
e1684.inc$id <- 1:n
e1684.inc$curei <- n * (1 - cure) - sapply(1:n, function(k) {
  cure.deletei <- min(update(KM, subset = -k)$surv)
  (n - 1) * (1 - cure.deletei)
})
head(e1684.inc)

library(geepack)
fit.inc <- geese(curei ~ AGE + TRT * SEX, data = e1684.inc,
                 jack = TRUE, scale.fix = TRUE, family = gaussian,
                 mean.link = "logit")

summary(fit.inc)

## Pseudo-residual
e1684.inc$resid <- e1684.inc$curei - 
  1 / (exp(-model.matrix(~ AGE + TRT * SEX, e1684.inc) %*% fit.inc$beta) + 1)
ggplot(e1684.inc, aes(x = interaction(TRT, SEX), y = resid, fill = interaction(TRT, SEX))) +
  geom_boxplot() + ylab("Residual") +
  scale_fill_discrete(labels = c("Control/Male", "Control/Female", 
                                 "Treatment/Male", "Treatment/Female")) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")

## Cure rate
newDat$cure.ps <- 1 - 1 / (1 + exp(-model.matrix(~ AGE + TRT * SEX, data = newDat) %*% fit.inc$b)) 
newDat

## Latency component
t0 <- quantile(e1684$FAILTIME[e1684$FAILCENS > 0], c(1:9 / 10, .95))
S <- KM$surv[findInterval(t0, KM$time)]
e1684.lat <- e1684[rep(1:n, each = length(t0)),]
e1684.lat$id <- rep(1:n, each = length(t0))
e1684.lat$Si <- n * (S - cure) / (1 - cure) - c(sapply(1:n, function(k) { 
  KM.reduce <- update(KM, subset = -k)
  S.deletei <- KM.reduce$surv[findInterval(t0, KM.reduce$time)]  
  cure.deletei <- min(KM.reduce$surv)
  (n - 1) * (S.deletei - cure.deletei) / (1 - cure.deletei)
}))
e1684.lat$Fi <- 1 - e1684.lat$Si
e1684.lat$t <- as.factor(rep(t0, n))
rownames(e1684.lat) <- NULL
head(e1684.lat)

fit.lat <- geese(Fi ~ 0 + t + AGE + TRT * SEX, data = e1684.lat, id = id, 
                jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "cloglog")
summary(fit.lat)

## Pseudo-residual
e1684.lat$resid <- 
  e1684.lat$Si - exp(-exp(model.matrix(~ 0 + t + AGE + TRT * SEX, e1684.lat) %*% fit.lat$beta))
levels(e1684.lat$t) <- paste0("t = ", unique(e1684.lat$t))
ggplot(subset(e1684.lat, t0 %in% t0[2:5 * 2]), 
       aes(x = interaction(TRT, SEX), y = resid, fill = interaction(TRT, SEX))) +
  facet_wrap(~ t) + geom_boxplot() + ylab("Residual") +
  scale_fill_discrete(labels = c("Control/Male", "Control/Female", 
                                 "Treatment/Male", "Treatment/Female")) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")

fit.inc$b
drop(foo$fit1$b)

fit.lat$b
drop(foo$fit2$b)

max(abs(fit.inc$b - drop(foo$fit1$b)))
max(abs(fit.lat$b - drop(foo$fit2$b)))

e
## #######################################################################
## Pseudo-observation approach for Bounded cumulative hazard model
## #######################################################################
## Long-term effect

e1684.long <- e1684
e1684.long$id <- 1:n
e1684.long$thetai <- n * (-log(cure)) - sapply(1:n, function(k) {
  cure.deletei <- min(update(KM, subset = -k)$surv)
  (n - 1) * (-log(cure.deletei))
})
head(e1684.long)

fit.long <- geese(thetai ~ AGE + TRT * SEX, data = e1684.long,
                  jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "log")
summary(fit.long)

## Pseudo residual
e1684.long$resid <- e1684.long$thetai - 
  exp(model.matrix(~ AGE + TRT * SEX, e1684.long) %*% fit.long$beta)
ggplot(e1684.long, aes(x = interaction(TRT, SEX), y = resid, fill = interaction(TRT, SEX))) +
  geom_boxplot() + ylab("Residual") +
  scale_fill_discrete(labels = c("Control/Male", "Control/Female", 
                                 "Treatment/Male", "Treatment/Female")) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")

## Short-term effect

e1684.short <- e1684[rep(1:n, each = length(t0)),]
e1684.short$id <- rep(1:n, each = length(t0))
e1684.short$Fi <- n * (log(S) / log(cure)) - c(sapply(1:n, function(k) {
  KM.reduce <- update(KM, subset = -k)
  S.deletei <- KM.reduce$surv[findInterval(t0, KM.reduce$time)]  
  cure.deletei <- min(KM.reduce$surv)
  (n - 1) * log(S.deletei) / log(cure.deletei) 
}))
e1684.short$t <- as.factor(rep(t0, n))
rownames(e1684.short) <- NULL
head(e1684.short)

fit.short <- geese(Fi ~ 0 + t + AGE + TRT * SEX, data = e1684.short,
             jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "cloglog")
summary(fit.short)

## Pseudo residual 
e1684.short$resid <- 1 - e1684.short$Fi - exp(-exp(
  model.matrix(~ 0 + t + AGE + TRT * SEX, e1684.short) %*% fit.short$beta))
levels(e1684.short$t) <- paste0("t = ", unique(e1684.short$t))
ggplot(subset(e1684.short, t0 %in% t0[2:5 * 2]), 
       aes(x = interaction(TRT, SEX), y = resid, fill = interaction(TRT, SEX))) +
  facet_wrap(~ t) + geom_boxplot() + ylab("Residual") +
  scale_fill_discrete(labels = c("Control/Male", "Control/Female", 
                                 "Treatment/Male", "Treatment/Female")) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")



fit.long$b
drop(foo1$fit1$b)

fit.short$b
drop(foo1$fit2$b)

max(abs(fit.long$b - drop(foo1$fit1$b)))
max(abs(fit.short$b - drop(foo1$fit2$b)))

library(MASS)
with(foo1$fit1, ginv(H) %*% M %*% ginv(H))
fit.long$vbeta




## #######################################################################
## Prints
## #######################################################################

head(e1684)

(foo1 <- pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684))
(foo2 <- pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684, model = 'p'))
(foo3 <- pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684, lambda1 = .02, lambda2 = .01))
(foo4 <- pCure(~ AGE + TRT * SEX, ~ AGE + TRT * SEX, FAILTIME, FAILCENS, e1684, model = 'p', lambda1 = .02, lambda2 = .02))

summary(foo1)
summary(foo2)
summary(foo3)
summary(foo4)

str(foo3)
## debug(pseudoCure:::fitPHMC)

pgee(thetai, X1, control$binit1, rep(1, n), control$exclude1, "logit", 
     control$penalty1, "ind", control$lambda1, control$eps, control$tol,  control$maxit)

x <- foo3

pvalTab <- function(pe, se, names = NULL) {
  if (is.null(se)) se <- NA
  tab <- cbind(Estimate = pe, StdErr = se,
               z.value = pe / se, p.value = 2 * pnorm(-abs(pe / se)))
  if (!is.null(names)) rownames(tab) <- names
  return(tab)
}

tab <- pvalTab(foo1$fit1$b, sqrt(diag(foo1$fit1$vb)))
tab2 <- tab
tab2[3,] <- "."

tab
tab2

printCoefmat(as.data.frame(tab), P.values = TRUE,
             has.Pvalue = TRUE, signif.legend = FALSE)

printCoefmat(as.data.frame(tab2), P.values = TRUE,
             has.Pvalue = TRUE, signif.legend = FALSE)

printCoefmat2(as.data.frame(tab), P.values = TRUE,
              has.Pvalue = TRUE, signif.legend = FALSE)



printCoefmat2 <- function(x, digits = max(3L, getOption("digits") - 2L),
                          signif.stars = getOption("show.signif.stars"), 
                          signif.legend = signif.stars,
                          dig.tst = max(1L, min(5L, digits - 1L)),
                          cs.ind = 1:k, tst.ind = k + 1, zap.ind = integer(), 
                          P.values = NULL,
                          has.Pvalue = nc >= 4L &&
                            length(cn <- colnames(x)) && substr(cn[nc], 1L, 3L)
                          %in% c("Pr(", "p-v"), eps.Pvalue = .Machine$double.eps, 
                          na.print = "NA", quote = FALSE, right = TRUE, ...) 
{
  if (is.null(d <- dim(x)) || length(d) != 2L) 
    stop("'x' must be coefficient matrix/data frame")
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  }
  else if (P.values && !has.Pvalue) 
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  }
  else xm <- data.matrix(x)
  k <- nc - has.Pvalue - (if (missing(tst.ind)) 
                            1
                          else length(tst.ind))
  if (!missing(cs.ind) && length(cs.ind) > k) 
    stop("wrong k / cs.ind")
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0])) 
                      floor(log10(range(acs[acs != 0], finite = TRUE)))
                    else 0
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - 
                                                    digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) 
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), 
                            digits = digits)
  if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) 
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) 
           ok[, -nc]
         else ok
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") 
    x1 <- chartr(dec, ".", x1)
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L, 
                                                                    digits - 1L))
  }
  if (any(ina)) 
    Cf[ina] <- na.print
  if (any(inan <- is.nan(xm))) 
    Cf[inan] <- "NaN"
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst, 
                                 eps = eps.Pvalue)
      signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                         symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
      }
    }
    else signif.stars <- FALSE
  }
  else signif.stars <- FALSE
  print.default(Cf, quote = quote, right = right, na.print = na.print, ...)
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, 
                                                       "legend"))) 
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w + 
                                                   4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}
