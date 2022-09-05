## #######################################################################
## Load and prepare the Dental data 
## #######################################################################

data(Teeth, package = "MST")
head(Teeth[,1:10])
Teeth$x15 <- 1 * !(Teeth$x15 == "No Crown")
Teeth$x16 <- 1 * !(Teeth$x16 == "No Endo Therapy")
Teeth$x17 <- 1 * !(Teeth$x17 == "No Implant")
Teeth$x18 <- 1 * !(Teeth$x18 == "No Bridge Pontic")
Teeth$x19 <- 1 * !(Teeth$x19 == "Missing")
Teeth$x20 <- 1 * !(Teeth$x20 == "Not Filled")
Teeth$x21 <- 1 * !(Teeth$x21 == "Not Decayed")
Teeth$x49 <- 1 * !(Teeth$x49 == "Male")
Teeth$x50 <- 1 * !(Teeth$x50 == "No Diabetes")
Teeth$x51 <- 1 * !(Teeth$x51 == "Never Had Tobacco")
Teeth$x52 <- 1 * Teeth$molar

## Extract first event
Teeth1 <- do.call(rbind, lapply(split(Teeth, Teeth$id),
                                function(d) d[which.min(d$time),]))

## Test for cure
npcure::testmz(time, event, Teeth1)


## #######################################################################
## Pseudo-observation approach for mixture cure model
## #######################################################################
n <- nrow(Teeth1)
tmax <- max(Teeth1$time[Teeth1$event > 0])
library(survival)
KM <- survfit(Surv(time, event) ~ 1, data = Teeth1)
(cure <- min(KM$surv))

## Incidence component
Teeth.inc <- Teeth1
Teeth.inc$id <- 1:n
Teeth.inc$curei <- n * (1 - cure) - sapply(1:n, function(k) {
  cure.deletei <- min(update(KM, subset = -k)$surv)
  (n - 1) * (1 - cure.deletei)
})
head(Teeth.inc[,1:10])

fn <- ~ x52 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 
  x11 + x12 + x13 + x14 + x15 + x16 + 
  x20 + x21 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
  x30 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + 
  x40 + x42 + x43 + x44 + x45 + x46 + x48 + x49 + x50 + x51

library(geepack)
fit.inc <- geese(as.formula(paste("curei ~", fn)[2]),
                 data = Teeth.inc, jack = TRUE, scale.fix = TRUE,
                 family = gaussian, mean.link = "logit")

summary(fit.inc)

## Variable selection
## to install modifiedPGEE from source code
## install.packages("modifiedPGEE_0.0.1.tar.gz", repo = NULL, type = "source") 
library(modifiedPGEE)
packageVersion("modifiedPGEE")

## 5 fold cross-validation
Teeth.inc$y <- Teeth.inc$curei
system.time(
  cv <- CVfit(as.formula(paste("y ~", fn)[2]), id = id, 
              data = subset(Teeth.inc, select = c("y", "id", all.vars(fn))),
              family = gaussian(link = "logit"), scale.fix = TRUE,
              fold = 5, lambda.vec = 1:45 / 50, pindex = 1:2)
)

fit.inc2 <- PGEE(as.formula(paste("y ~", fn)[2]), id = id, data = Teeth.inc, 
                 family = gaussian(link = "logit"), scale.fix = TRUE,
                 beta_int = fit.inc$beta, lambda = cv$lam.opt, pindex = 1:2)

fit.inc2$coef[abs(fit.inc2$coef) > 1e-3]

## Latency component
t0 <- quantile(Teeth1$time[Teeth1$event > 0], c(1:9 / 10, .95))
S <- KM$surv[findInterval(t0, KM$time)]
Teeth.lat <- data.frame(
  id = rep(1:n, each = length(t0)), 
  Si = n * (S - cure) / (1 - cure) - c(sapply(1:n, function(k) { 
    KM.reduce <- update(KM, subset = -k)
    S.deletei <- KM.reduce$surv[findInterval(t0, KM.reduce$time)]  
    cure.deletei <- min(KM.reduce$surv)
    (n - 1) * (S.deletei - cure.deletei) / (1 - cure.deletei)
  })),
  Teeth1[rep(1:n, each = length(t0)),],
  t = kronecker(rep(1, n), diag(length(t0))))
Teeth.lat$Fi <- 1 - Teeth.lat$Si
rownames(Teeth.lat) <- NULL
head(Teeth.lat[,1:10])

fn2 <- ~ 0 + t.1 + t.2 + t.3 + t.4 + t.5 + t.6 + t.7 + t.8 + t.9 + t.10 +
  x52 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 
  x11 + x12 + x13 + x14 + x15 + x16 + 
  x20 + x21 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
  x30 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + 
  x40 + x42 + x43 + x44 + x45 + x46 + x48 + x49 + x50 + x51

fit.lat <- geese(as.formula(paste("Fi ~", fn2)[2]), data = Teeth.lat, id = id, 
                jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "cloglog")
summary(fit.lat)

## Variable selection
Teeth.lat$y <- Teeth.lat$Fi
system.time(
  cv <- CVfit(as.formula(paste("y ~", fn2)[2]), id = id, 
              data = subset(Teeth.lat, select = c("y", "id", all.vars(fn2))),
              family = gaussian(link = "cloglog"), scale.fix = TRUE,
              fold = 5, lambda.vec = 1:45 / 50, pindex = 1:11)
)

fit.lat2 <- PGEE(as.formula(paste("y ~", fn2)[2]), id = id, data = Teeth.lat, 
                 family = gaussian(link = "cloglog"), scale.fix = TRUE,
                 beta_int = fit.lat$beta, lambda = cv$lam.opt, pindex = 1:11)

fit.lat2$coef[abs(fit.lat2$coef) > 1e-3]

## #######################################################################
## With the intsurv package
## #######################################################################

library(intsurv)
system.time(
  fit.intsurv <- cox_cure_net(fn, fn, data = Teeth1,
                              time = time, event = event,
                              surv_nlambda = 10, cure_nlambda = 10,
                              surv_alpha = 1, cure_alpha = 1,
                              surv_l1_penalty_factor = c(-1e-10, rep(43)),
                              cure_l1_penalty_factor = c(1, -1e10, rep(43))))
lapply(coef(fit.intsurv), function(x) x[x != 0])

## #######################################################################
## Pseudo-observation approach for Bounded cumulative hazard model
## #######################################################################
## Long-term effect

Teeth.long <- Teeth1
Teeth.long$id <- 1:n
Teeth.long$thetai <- n * (-log(cure)) - sapply(1:n, function(k) {
  cure.deletei <- min(update(KM, subset = -k)$surv)
  (n - 1) * (-log(cure.deletei))
})
head(Teeth.long[,1:10])
fit.long <- geese(as.formula(paste("thetai ~", fn)[2]), data = Teeth.long,
                  jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "log")
summary(fit.long)

## Variable selection
Teeth.long$y <- Teeth.long$thetai
system.time(
  cv <- CVfit(as.formula(paste("y ~", fn)[2]), id = id, 
              data = subset(Teeth.long, select = c("y", "id", all.vars(fn))),
              family = gaussian(link = "log"), scale.fix = TRUE,
              fold = 5, lambda.vec = 1:45 / 50, pindex = 1:2)
)

fit.long2 <- PGEE(as.formula(paste("y ~", fn)[2]), id = id, data = Teeth.long, 
                 family = gaussian(link = "log"), scale.fix = TRUE,
                 beta_int = fit.long$beta, lambda = cv$lam.opt, pindex = 1:2)

fit.long2$coef[abs(fit.long2$coef) > 1e-3]

## Short-term effect

Teeth.short <- data.frame(
  id = rep(1:n, each = length(t0)), 
  Fi = n * (log(S) / log(cure)) - c(sapply(1:n, function(k) {
    KM.reduce <- update(KM, subset = -k)
    S.deletei <- KM.reduce$surv[findInterval(t0, KM.reduce$time)]  
    cure.deletei <- min(KM.reduce$surv)
    (n - 1) * log(S.deletei) / log(cure.deletei) 
  })),
  Teeth1[rep(1:n, each = length(t0)), ], 
  t = kronecker(rep(1, n), diag(length(t0))))
rownames(Teeth.short) <- NULL
head(Teeth.short[,1:10])

fit.short <- geese(as.formula(paste("Fi ~", fn2)[2]), data = Teeth.short,id = id, 
             jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "cloglog")
summary(fit.short)

## Variable selection
Teeth.short$y <- Teeth.short$Fi
system.time(
  cv <- CVfit(as.formula(paste("y ~", fn2)[2]), id = id, 
              data = subset(Teeth.short, select = c("y", "id", all.vars(fn2))),
              family = gaussian(link = "cloglog"), scale.fix = TRUE,
              fold = 5, lambda.vec = 1:45 / 50, pindex = 1:11)
)

fit.short2 <- PGEE(as.formula(paste("y ~", fn2)[2]), id = id, data = Teeth.short, 
                 family = gaussian(link = "cloglog"), scale.fix = TRUE,
                 beta_int = fit.short$beta, lambda = cv$lam.opt, pindex = 1:11)

fit.short2$coef[abs(fit.short2$coef) > 1e-3]
