## #######################################################################
## Load and prepare the Melanoma data 
## #######################################################################

data(e1684, package = "smcure")
head(e1684)
e1684 <- na.omit(e1684)

## Checking for cure

library(survival)
library(survminer)
ggsurvplot(survfit(Surv(FAILTIME, FAILCENS) ~ TRT + SEX, data = e1684),
           legend.title = "", 
           legend.labs = c("Control/Male", "Control/Female", "Treatment/Male", "Treatment/Female")) +  
  labs(title = "The melanoma data")

npcure::testmz(FAILTIME, FAILCENS, e1684)

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
                jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "logit")

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
