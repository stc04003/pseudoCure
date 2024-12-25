
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.2.0-6666ff.svg)](https://cran.r-project.org/)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--12--25-yellowgreen.svg)](/commits/master)

------------------------------------------------------------------------

## **pseudoCure**

------------------------------------------------------------------------

## pseudoCure: Analysis of survival data with cure fraction and variable selection: A pseudo-observations approach

The **pseudoCure** package implements a pseudo-observation approach for
survival data with a cure fraction. The modeling framework is based on
the Cox proportional hazards mixture cure model and the bounded
cumulative hazard model.

## Installation

Install and load the package from GitHub using

``` r
> devtools::install_github("stc04003/pseudoCure")
> library(pseudoCure)
> packageVersion("pseudoCure")
```

## Reference

Su, C.-L., Chiou, S., Lin, F.-C., and Platt, R. W. (2022) Analysis of
survival data with cure fraction and variable selection: A
pseudo-observations approach *Statistical Methods in Medical Research*,
**31**(11): 2037–2053.
