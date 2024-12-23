% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Teeth500}
\alias{Teeth500}
\title{Dental data for illustration}
\format{
A data frame containing the following variables:
\describe{
  \item{time}{tooth survival time subject to right censoring.}
  \item{event}{Tooth loss status: 1 = lost, 0 = not lost.}
  \item{molar}{Molar indicator; 1 = molar tooth, 0 = non-molar tooth.}
  \item{mobil}{Mobility score, on a scale from 0 to 5.}
  \item{bleed}{Bleeding on probing, expressed as a percentage.}
  \item{plaque}{Plaque score, expressed as a percentage.}
  \item{pocket}{Periodontal probing depth.}
  \item{cal}{Clinical Attachment Level.}
  \item{fgm}{Free Gingival Margin.}
  \item{filled}{Number of filled surfaces.}
  \item{decay_new}{New decayed surfaces.}
  \item{decay_recur}{Recurrent decayed surfaces.}
  \item{crown}{Crown indicator; 1 = tooth has a crown, 0 = no crown.}
  \item{endo}{Endodontic therapy indicator; 1 = endo therapy performed, 0 = no endo therapy.}
  \item{filled_tooth}{Filled tooth indicator; 1 = filled, 0 = not filled.}
  \item{decayed_tooth}{Decayed tooth indicator; 1 = decayed, 0 = not decayed.}
  \item{total_tooth}{Total number of teeth.}
  \item{gender}{Gender; 1 = male, 0 = female}
  \item{diabetes}{Diabetes indicator; 1 = diabetes, 0 = no diabetes.}
  \item{tobacco_ever}{Tobacco use indicator; 1 = had tobacco use, 0 = never had tobacco use.}
}

A data frame with 500 observations and 20 variables.
}
\usage{
data(Teeth500)
}
\description{
Data on the survival of teeth with many predictors
}
\details{
The data is a subset of the original dataset included in the \code{MST} package
under the name \code{Teeth}.
This subset contains the time to the first tooth loss due to periodontal reasons.
}
\references{
Calhoun, Peter and Su, Xiaogang and Nunn, Martha and Fan, Juanjuan (2018) Constructing Multivariate Survival Trees: The MST Package for R. \emph{Journal of Statistical Software}, \bold{83}(12).
}
