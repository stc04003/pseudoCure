## ###########################################################################
## Load and prepare the Melanoma data 
## If latency component is involved, call minSi4(), otherwise, call minSi()
## If short-term component is involved, call minSi4(), otherwise, call minSi()
## ###########################################################################

library(survival)
library(microbenchmark)
library(Rcpp)
library(geepack)

sourceCpp("RcppCodes.cpp")

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
Teeth1 <- do.call(rbind, lapply(split(Teeth, Teeth$id), function(d) d[which.min(d$time),]))

## Test for cure
npcure::testmz(time, event, Teeth1)

n <- nrow(Teeth1)
tmax <- max(Teeth1$time[Teeth1$event > 0])
t0 <- quantile(Teeth1$time[Teeth1$event > 0], c(1:9 / 10, .95))


## ####################################
## Without variable selection
## ####################################

## Incidence component
KMs2 <- minSi4(Teeth$time, Teeth$event, c(t0, max(Teeth$time)))
S <- KMs2[1:length(t0), n + 1]
Si <- KMs2[1:length(t0), 1:n]
cure <- KMs2[length(t0) + 1, n + 1]
curei <- KMs2[length(t0) + 1, 1:n]

Teeth.inc <- Teeth1
Teeth.inc$id <- 1:n
Teeth.inc$curei <- n * (1 - cure) - (n - 1) * (1 - curei)

head(Teeth.inc[,1:10])
names(Teeth.inc)

fn <- ~ x52 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 
  x11 + x12 + x13 + x14 + x15 + x16 + 
  x20 + x21 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
  x30 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + 
  x40 + x42 + x43 + x44 + x45 + x46 + x48 + x49 + x50 + x51

fit.inc <- geese(as.formula(paste("curei ~", fn)[2]),
                 data = Teeth.inc, jack = TRUE, scale.fix = TRUE,
                 family = gaussian, mean.link = "logit")

summary(fit.inc)





## Latency component

Teeth.lat <- data.frame(
    id = rep(1:n, each = length(t0)),
    Si = n * (S - cure) / (1 - cure) - c((n - 1) * (Si - curei) / (1 - curei)),
    Teeth[rep(1:n, each = length(t0)),],
    t = kronecker(rep(1, n), diag(length(t0))))
Teeth.lat$Fi <- 1 - Teeth.lat$Si
rownames(Teeth.lat) <- NULL

head(Teeth.lat[,1:10])
names(Teeth.lat)

fn2 <- ~ 0 + t.1 + t.2 + t.3 + t.4 + t.5 + t.6 + t.7 + t.8 + t.9 + t.10 +
  x52 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 
  x11 + x12 + x13 + x14 + x15 + x16 + 
  x20 + x21 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + 
  x30 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + x39 + 
  x40 + x42 + x43 + x44 + x45 + x46 + x48 + x49 + x50 + x51

fit.lat <- geese(as.formula(paste("Fi ~", fn2)[2]), data = Teeth.lat, id = id, 
                jack = TRUE, scale.fix = TRUE, family = gaussian, mean.link = "cloglog")
summary(fit.lat)

