## ###########################################################################
## Load and prepare the Melanoma data 
## If latency component is involved, call minSi4(), otherwise, call minSi()
## If short-term component is involved, call minSi4(), otherwise, call minSi()
## ###########################################################################

library(survival)
library(microbenchmark)
library(Rcpp)
library(geepack)

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
## Incidence component
## ####################################

KMs2 <- minSi4(Teeth$time, Teeth$event, c(t0, max(Teeth$time)))
