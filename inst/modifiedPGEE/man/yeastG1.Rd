\name{yeastG1}
\alias{yeastG1}
\docType{data}
\title{Yeast cell-cycle gene expression data}
\description{A yeast cell-cycle gene expression data set collected in the CDC15 experiment of 
Spellman et al. (1998) where genome-wide mRNA levels of 6178 yeast 
open reading frames (ORFs) in a two cell-cycle period were measured at M/G1-G1-S-G2-M stages.
However, to better understand the phenomenon underlying cell-cycle process, it is important to identify transcription 
factors (TFs) that regulate the gene expression levels of cell cycle-regulated genes. In this study, 
we presented a subset of 283 cell-cycled-regularized genes observed over 4 time points at G1 stage and 
the standardized binding probabilities of a total of 96 TFs obtained from a mixture model approach of 
Wang et al. (2007) based on the ChIP data of Lee et al. (2002).}

\usage{data("yeastG1")}

\details{A data frame with 1132 observations (283 cell-cycled-regularized genes observed over 4 time points) 
with 99 variables (e.g., id, y, time, and 96 TFs).}

\references{
Lee, T.I., Rinaldi, N.J., Robert, F., Odom, D.T., Bar-Joseph, Z., Gerber, G.K., Hannett, N.M., 
Harbison, C.T., Thompson, C.M., Simon, I., et al. (2002). Transcriptional regulatory networks in
Saccharomyces cerevisiae. \emph{Science}, \bold{298}, 799--804.

Spellman, P.T., Sherlock, G., Zhang, M.Q., Iyer, V.R., Anders, K., Eisen, M.B., Brown, P.O., 
Botstein, D., and Futcher, B. (1998). Comprehensive identification of cell cycle regulated genes of the yeast
Saccharomyces cerevisiae by microarray hybridization. \emph{Molecular Biology of Cell}, \bold{9},
3273--3297.

Wang, L., Chen, G., and Li, H. (2007). Group SCAD regression analysis for microarray time course
gene expression data. \emph{Bioinformatics}, \bold{23}, 1486--1494.

Wang, L., Zhou, J., and Qu, A. (2012). Penalized generalized estimating equations 
for high-dimensional longitudinal data anaysis. \emph{Biometrics}, \bold{68}, 353--360.

}

\examples{
\dontrun{
library(PGEE)
# load data
data(yeastG1)
data <- yeastG1
# get the column names
colnames(data)[1:9]
# see some portion of yeast G1 data
head(data,5)[1:9]

# define the input arguments
formula <- "y ~.-id"
family <- gaussian(link = "identity")
lambda.vec <- seq(0.01,0.2,0.01)
# find the optimum lambda
cv <- CVfit(formula = formula, id = id, data = data, family = family, scale.fix = TRUE,
scale.value = 1, fold = 4, lambda.vec = lambda.vec, pindex = c(1,2), eps = 10^-6,
maxiter = 30, tol = 10^-6)
# print the results
print(cv)

# see the returned values by CVfit
names(cv)
# get the optimum lambda
cv$lam.opt

#fit the PGEE model
myfit1 <- PGEE(formula = formula, id = id, data = data, na.action = NULL,
family = family, corstr = "independence", Mv = NULL,
beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
scale.value = 1, lambda = cv$lam.opt, pindex = c(1,2), eps = 10^-6, 
maxiter = 30, tol = 10^-6, silent = TRUE)

# get the values returned by myfit object
names(myfit1)
# get the values returned by summary(myfit) object
names(summary(myfit1))
# see a portion of the results returned by coef(summary(myfit1))
head(coef(summary(myfit1)),7)

# see the variables which have non-zero coefficients
index1 <- which(abs(coef(summary(myfit1))[,"Estimate"]) > 10^-3)
names(abs(coef(summary(myfit1))[index1,"Estimate"]))

# see the PGEE summary statistics of these non-zero variables
coef(summary(myfit1))[index1,]

# fit the GEE model
myfit2 <- MGEE(formula = formula, id = id, data = data, na.action = NULL,
family = family, corstr = "independence", Mv = NULL,
beta_int = c(rep(0,dim(data)[2]-1)), R = NULL, scale.fix = TRUE,
scale.value = 1, maxiter = 30, tol = 10^-6, silent = TRUE)

# get the GEE summary statistics of the variables that turned out to be
# non-zero in PGEE analysis
coef(summary(myfit2))[index1,]

# see the significantly associated TFs in PGEE analysis
names(which(abs(coef(summary(myfit1))[index1,"Robust z"]) > 1.96))

# see the significantly associated TFs in GEE analysis
names(which(abs(coef(summary(myfit2))[,"Robust z"]) > 1.96))

}
}
