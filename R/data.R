#' Dental data for illustration
#'
#' Data on the survival of teeth with many predictors
#' 
#' @format A data frame containing the following variables:
#' \describe{
#'   \item{time}{tooth survival time subject to right censoring.}
#'   \item{event}{Tooth loss status: 1 = lost, 0 = not lost.}
#'   \item{molar}{Molar indicator; TRUE = molar tooth, FALSE = non-molar tooth.}
#'   \item{mobil}{Mobility score, on a scale from 0 to 5.}
#'   \item{bleed}{Bleeding on probing, expressed as a percentage.}
#'   \item{plaque}{Plaque score, expressed as a percentage.}
#'   \item{pocket}{Periodontal probing depth (tooth-level mean).}
#'   \item{cal}{Clinical Attachment Level (tooth-level mean).}
#'   \item{fgm}{Free Gingival Margin (tooth-level mean).}
#'   \item{filled}{Number of filled surfaces.}
#'   \item{decay_new}{New decayed surfaces.}
#'   \item{decay_recur}{Recurrent decayed surfaces.}
#'   \item{crown}{Crown indicator; TRUE = tooth has a crown, FALSE = no crown.}
#'   \item{endo}{Endodontic therapy indicator; TRUE = endo therapy performed, FALSE = no endo therapy.}
#'   \item{filled_tooth}{Filled tooth indicator; TRUE = filled, FALSE = not filled.}
#'   \item{decayed_tooth}{Decayed tooth indicator; TRUE = decayed, FALSE = not decayed.}
#'   \item{filled_pct}{Mean number of filled surfaces.}
#'   \item{decay_new_pct}{Mean number of new decayed surfaces.}
#'   \item{decay_recur_pct}{Mean number of recurrent decayed surfaces.}
#'   \item{filled_tooth_pct}{Percentage of filled teeth.}
#'   \item{decayed_tooth_pct}{Percentage of decayed teeth.}
#'   \item{missing_tooth_pct}{Percentage of missing teeth.}
#'   \item{total_tooth}{Total number of teeth.}
#'   \item{gender}{Gender; Female vs Male.}
#'   \item{diabetes}{Diabetes indicator; TRUE = diabetes, FALSE = no diabetes.}
#'   \item{tobacco_ever}{Tobacco use indicator; TRUE = had tobacco use, FALSE = never had tobacco use.}
#' }
#'
#' @details
#' The data is a subset of the original dataset included in the \code{MST} package
#' under the name \code{Teeth}.
#' This subset contains the time to the first tooth loss due to periodontal reasons.
#' 
#' @usage data(Teeth500)
#' @docType data
#' @name Teeth500
#' @rdname Teeth500
#' @format A data frame with 500 observations and 26 variables.
#'
#' @references
#' Calhoun, Peter and Su, Xiaogang and Nunn, Martha and Fan, Juanjuan (2018) Constructing Multivariate Survival Trees: The MST Package for R. \emph{Journal of Statistical Software}, \bold{83}(12).
#'
NULL



