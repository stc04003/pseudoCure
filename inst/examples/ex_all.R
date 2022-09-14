## ##################################################################################
## All examples used in testing pCure()
## Not all are included in .Rd due to timing restriction
## ##################################################################################

library(pseudoCure)

## A melanoma data example
data(e1684, package = "smcure")
e1684 <- na.omit(e1684)

## Define formula
fn <- ~ AGE + TRT * SEX; set.seed(123)

## Unpenalized PHMC model (~ 0.015s)
summary(fit1 <- pCure(fn, fn, FAILTIME, FAILCENS, e1684))

## Penalized PHMC model with tuning parameters selected by 10-fold cross validation
## User specifies the range of tuning parameters (~0.190s)
summary(update(fit1, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200))
## Auto selection of the range of tuning parameters (~2.344s)
summary(update(fit1, lambda1 = "auto", lambda2 = "auto"))
## Penalized PHMC model given tuning parameters (~0.044s)
summary(update(fit1, lambda1 = 0.006, lambda2 = 0.022))

## All of the above but with PHPH model
summary(fit2 <- update(fit1, model = 'p')) ## ~0.013s
summary(update(fit2, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200)) ## ~0.758s
summary(update(fit2, lambda1 = "auto", lambda2 = "auto")) ## ~6.053
summary(update(fit2, lambda1 = 0.006, lambda2 = 0.022)) ## ~0.019s

## A dental data example
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
Teeth1 <- do.call(rbind, lapply(split(Teeth, Teeth$id),
                                function(d) d[which.min(d$time),]))

fn <- ~ x52 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 
    x11 + x12 + x13 + x14 + x15 + x16 + 
    x20 + x21 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
    x30 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + 
    x40 + x42 + x43 + x44 + x45 + x46 + x48 + x49 + x50 + x51

## Unpenalized PHMC model (~ 7.792s)
summary(fit1 <- pCure(fn, fn, time, event, Teeth1))

## Penalized PHMC model with tuning parameters selected by 10-fold cross validation
## User specifies the range of tuning parameters (~3.817s)
summary(update(fit1, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200))
## Auto selection of the range of tuning parameters (~48.06s)
summary(update(fit1, lambda1 = "auto", lambda2 = "auto"))
## Penalized PHMC model given tuning parameters (~0.336s)
summary(update(fit1, lambda1 = 0.006, lambda2 = 0.022))

## All of the above but with PHPH model
summary(fit2 <- update(fit1, model = 'p')) ## ~7.696s
summary(update(fit2, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200)) ## ~
summary(update(fit2, lambda1 = "auto", lambda2 = "auto")) ## ~1.790s
summary(update(fit2, lambda1 = 0.006, lambda2 = 0.022)) ## ~


now <- Sys.time()
summary(update(fit2, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200)) ## ~8.882s
then <- Sys.time()
then - now

