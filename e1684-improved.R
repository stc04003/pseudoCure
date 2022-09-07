## #######################################################################
## Load and prepare the Melanoma data 
## #######################################################################
library(survival)
library(microbenchmark)
library(Rcpp)
library(geepack)

sourceCpp("RcppCodes.cpp")

data(e1684, package = "smcure")
head(e1684)
e1684 <- na.omit(e1684)

## ####################################
## Incidence component
## ####################################
n <- nrow(e1684)
tmax <- max(e1684$FAILTIME[e1684$FAILCENS > 0])

KMs <- minSi(e1684$FAILTIME, e1684$FAILCENS)
cure <- KMs[length(KMs)]

e1684.inc <- e1684
e1684.inc$id <- 1:n
e1684.inc$curei <- n * (1 - cure) - (n - 1) * (1 - KMs[1:n])
head(e1684.inc)

fit.inc <- geese(curei ~ AGE + TRT * SEX, data = e1684.inc,
                jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "logit")

summary(fit.inc)


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

head(drop(minSi(e1684$FAILTIME, e1684$FAILCENS)))
head(drop(minSi1(e1684$FAILTIME, e1684$FAILCENS)))
tail(drop(minSi(e1684$FAILTIME, e1684$FAILCENS)))
tail(drop(minSi1(e1684$FAILTIME, e1684$FAILCENS)))
            
t0 <- quantile(e1684$FAILTIME[e1684$FAILCENS > 0], c(1:9 / 10, .95))

KM <- survfit(Surv(FAILTIME, FAILCENS) ~ 1, data = e1684)
S <- KM$surv[findInterval(t0, KM$time)]

KM$surv[findInterval(t0, KM$time)]

sapply(findInterval(t0, KM$time), function(e) 
  drop(minSi1(e1684$FAILTIME, e1684$FAILCENS, e))[n + 1])


Si <- sapply(1:nrow(e1684), function(k) { 
  KM.reduce <- update(KM, subset = -k)
  findInterval(t0, KM.reduce$time)
})

e

## #######################################################################
## Try
## #######################################################################

n <- 1000
tt <- rexp(n)
dd <- sample(0:1, n, T)
## dd <- rep(1, n)
tt <- sort(tt)

KM <- survfit(Surv(tt, dd) ~ 1)
all.equal(minS(tt, dd), min(KM$surv))
all.equal(minS(tt, dd), drop(minSi(tt, dd))[n + 1])
all.equal(drop(minSi(tt, dd))[1:n], sapply(1:n, function(k) min(update(KM, subset = -k)$surv)))

tt <- round(tt, 1)

KM <- survfit(Surv(tt, dd) ~ 1)
all.equal(minS(tt, dd), min(KM$surv))
all.equal(minS(tt, dd), drop(minSi(tt, dd))[n + 1])
all.equal(drop(minSi(tt, dd))[1:n], sapply(1:n, function(k) min(update(KM, subset = -k)$surv)))

drop(minSi(tt, dd))[1:n]
sapply(1:n, function(k) min(update(KM, subset = -k)$surv))

cbind(tt, dd)

microbenchmark(drop(minS(tt, dd)), drop(minSi(tt, dd)))

microbenchmark(drop(minS(tt, dd)),
               drop(minSi(tt, dd)),
               sapply(1:n, function(k) min(update(KM, subset = -k)$surv)))

drop(minSi(tt, dd)) 
drop(minSi2(tt, dd))
all.equal(drop(minSi(tt, dd)), drop(minSi2(tt, dd)))

microbenchmark(drop(minSi(tt, dd)), drop(minSi2(tt, dd)))


sourceCpp("RcppCodes.cpp")

b <- 6
yy <- sort(rexp(b))
## cc <- rep(1, 6)
cc <- sample(0:1, b, T)

drop(minSi(yy, cc))
drop(minSi1(yy, cc))
