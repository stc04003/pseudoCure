## A melanoma data example
data(e1684, package = "smcure")
e1684 <- na.omit(e1684)

## Define formula & fix random seed 
fn <- ~ AGE + TRT * SEX; set.seed(123)

## Unpenalized PHMC model (~ 0.046s)
summary(fit1 <- pCure(fn, fn, FAILTIME, FAILCENS, e1684))

## Penalized PHMC model with tuning parameters selected by 10-fold cross validation
## User specifies the range of tuning parameters
summary(update(fit1, lambda1 = 1:10 / 200, lambda2 = 1:10 / 200))
## Auto selection of the range of tuning parameters
summary(update(fit1, lambda1 = "auto", lambda2 = "auto"))

## Penalized PHMC model given tuning parameters
summary(update(fit1, lambda1 = 0.006, lambda2 = 0.022))
