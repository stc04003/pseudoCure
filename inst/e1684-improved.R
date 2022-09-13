## ###########################################################################
## Load and prepare the Melanoma data 
## If latency component is involved, call minSi4(), otherwise, call minSi()
## If short-term component is involved, call minSi4(), otherwise, call minSi()
## minSi4() -> pseudoKM1
## minSi() -> pseudoKM
## ###########################################################################

library(survival)
library(microbenchmark)
library(Rcpp)
library(geepack)
library(MASS)

sourceCpp("RcppCodes.cpp")
sourceCpp("PGEE.cpp")

data(e1684, package = "smcure")
head(e1684)
e1684 <- na.omit(e1684)
n <- nrow(e1684)
tmax <- max(e1684$FAILTIME[e1684$FAILCENS > 0])
t0 <- quantile(e1684$FAILTIME[e1684$FAILCENS > 0], c(1:9 / 10, .95))

## ####################################
## Incidence component
## ####################################

KMs <- minSi(e1684$FAILTIME, e1684$FAILCENS)
cure <- KMs[length(KMs)]

e1684.inc <- e1684
e1684.inc$id <- 1:n
e1684.inc$curei <- n * (1 - cure) - (n - 1) * (1 - KMs[1:n])
head(e1684.inc)

fit.inc <- geese(curei ~ AGE + TRT * SEX, data = e1684.inc,
                 jack = TRUE, scale.fix = TRUE, family = gaussian,
                 mean.link = "logit")

summary(fit.inc)



## If we know we are doing latent compoent
KMs2 <- minSi4(e1684$FAILTIME, e1684$FAILCENS, c(t0, max(e1684$FAILTIME)))
e1684.inc <- e1684
e1684.inc$id <- 1:n
e1684.inc$curei <- n * (1 - cure) - (n - 1) * (1 - KMs2[length(t0) + 1, 1:n])
head(e1684.inc)

fit.inc <- geese(curei ~ AGE + TRT * SEX, data = e1684.inc,
                jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "logit")

summary(fit.inc)

## Use pgee
fit.inc2 <- gee(y = e1684.inc$curei,
                X = model.matrix(~ AGE + TRT * SEX, data = e1684.inc),
                b0 = rep(0, 5), nt = rep(1, 284), pindex = rep(1, 5),
                glmlink = "logit", corstr = "ind",
                lambda = 0, eps = 1e-6, tol = 1e-7, maxit = 30)

fit.inc$beta
fit.inc2$b

with(fit.inc2, ginv(H + 284 * E) %*% M %*% ginv(H + 284 * E))
fit.inc$vbeta
     
## Pseudo-residual
e1684.inc$resid <- e1684.inc$curei - 
  1 / (exp(-model.matrix(~ AGE + TRT * SEX, e1684.inc) %*% fit.inc$beta) + 1)
ggplot(e1684.inc, aes(x = interaction(TRT, SEX), y = resid,
                      fill = interaction(TRT, SEX))) +
  geom_boxplot() + ylab("Residual") +
  scale_fill_discrete(labels = c("Control/Male", "Control/Female", 
                                 "Treatment/Male", "Treatment/Female")) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "bottom")

## Cure rate
newDat$cure.ps <-
  1 - 1 / (1 + exp(-model.matrix(~ AGE + TRT * SEX, data = newDat) %*% fit.inc$b)) 
newDat

## ####################################
## Latency component
## ####################################

S <- KMs2[1:length(t0), n + 1]
Si <- KMs2[1:length(t0), 1:n]
cure <- KMs2[length(t0) + 1, n + 1]
curei <- KMs2[length(t0) + 1, 1:n]

e1684.lat <- e1684[rep(1:n, each = length(t0)),]
e1684.lat$id <- rep(1:n, each = length(t0))
## e1684.lat$Si <- n * (S - cure) / (1 - cure) - c((n - 1) * (Si - curei) / (1 - curei))
e1684.lat$Si <- n * (S - cure) / (1 - cure) - c((n - 1) * t((t(Si) - curei) / (1 - curei)))
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


str(fit.lat2 <- gee0(e1684.lat$Fi, model.matrix(~ 0 + t + AGE + TRT * SEX, data = e1684.lat),
                     rep(0, length(fit.lat$b)), rep(length(t0), nrow(e1684)),
                     "cloglog", "ind", 1e-6, 50))
fit.lat$b
fit.lat2$b

## ####################################
## Long-term effect
## ####################################

e1684.long <- e1684
e1684.long$id <- 1:n
e1684.long$thetai <- n * (-log(cure)) - (n - 1) * (-log(curei))
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


## ####################################
## Short-term effect
## ####################################

e1684.short <- e1684[rep(1:n, each = length(t0)),]
e1684.short$id <- rep(1:n, each = length(t0))
e1684.short$Fi <- n * (log(S) / log(cure)) - c((n - 1) * log(Si) / log(curei))
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


library(intsurv)
cox_cure

head(e1684)

lapply(pCure(~ AGE + SEX, time = FAILTIME, status = FAILCENS, data = e1684), head)
lapply(pCure(~ AGE + SEX, ~ AGE, time = FAILTIME, status = FAILCENS, data = e1684), head)
lapply(pCure(~ AGE + SEX, ~ AGE - 1, time = FAILTIME, status = FAILCENS, data = e1684), head)

pCure <- function(formula1, formula2, time, status, data, subset,
                  model = c("mixture", "promotion")) {
    model <- match.arg(model)
    if (missing(formula1)) stop("Argument 'formula' is required.")
    if (missing(time)) stop("Argument 'time' is required.")
    if (missing(status)) stop("Argument 'status' is required.")
    if (missing(data)) data <- environment(formula)
    if (!missing(subset)) {
        sSubset <- substitute(subset)
        subIdx <- eval(sSubset, data, parent.frame())
        if (!is.logical(subIdx)) stop("'subset' must be logical")
        subIdx <- subIdx & !is.na(subIdx)
        data <- data[subIdx, ]
    }
    mf <- match.call(expand.dots = FALSE)    
    mf <- mf[match(c("formula1", "data", "time", "status"), names(mf), 0L)]
    mf$data <- data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    time <- as.numeric(model.extract(mf, time))
    status <- as.numeric(model.extract(mf, status))
    mm1 <- mm2 <- stats::model.matrix(formula1, data = mf)
    if (!missing(formula2)) {
        mf <- match.call(expand.dots = FALSE)
        mf <- mf[match(c("formula2", "data"), names(mf), 0L)]
        mf$data <- data
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- quote(stats::model.frame)
        mf <- eval(mf, parent.frame())
        mm2 <- stats::model.matrix(formula2, data = mf)
    }
    list(mm1 = mm1, mm2 = mm2, time = time, status = status)
}

debug(pCure)
