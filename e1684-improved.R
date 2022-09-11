## ###########################################################################
## Load and prepare the Melanoma data 
## If latency component is involved, call minSi4(), otherwise, call minSi()
## If short-term component is involved, call minSi4(), otherwise, call minSi()
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
e1684.lat$Si <- n * (S - cure) / (1 - cure) - c((n - 1) * (Si - curei) / (1 - curei))
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
