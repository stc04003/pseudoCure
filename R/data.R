#' Dental data for illustration
#'
#' @description A data on the survival of teeth with many predictors
#' The variables are as follows:
#' \describe{
#'   \item{time}{tooth survival time subject to right censoring}
#'   \item{event}{Tooth loss status: 1 = lost, 0 = not lost}
#'   \item{molar}{Indicates whether the tooth is a molar (logical).}
#'   \item{mobil}{Mobility score (on a scale 0-5)}
#'   \item{bleed}{Bleeding on probing (percentage)}
#'   \item{plaque}{Plaque score (percentage)}
#'   \item{pocket}{Periodontal probing depth (tooth-level mean, in mm)}
#'   \item{cal}{Clinical attachment level (tooth-level mean, in mm)}
#'   \item{fgm}{Free gingival margin (tooth-level mean, in mm)}
#'   \item{filled}{Number of filled surfaces on the tooth}
#'   \item{decay_new}{Number of newly decayed surfaces}
#'   \item{decay_recur}{Number of recurrently decayed surfaces}
#'   \item{crown}{Indicates whether the tooth has a crown (factor)}
#'   \item{endo}{Indicates whether the tooth has undergone endodontic therapy (factor)}
#'   \item{filled_tooth}{Indicates whether the tooth is filled (factor)}
#'   \item{decayed_tooth}{Indicates whether the tooth is decayed (factor)}
#'   \item{filled_sum}{Total number of filled surfaces across all teeth}
#'   \item{decay_new_sum}{Total number of newly decayed surfaces across all teeth}
#'   \item{decay_recur_sum}{Total number of recurrently decayed surfaces across all teeth}
#'   \item{filled_tooth_sum}{Total number of filled teeth}
#'   \item{decayed_tooth_sum}{Total number of decayed teeth}
#'   \item{missing_tooth_sum}{Total number of missing teeth}
#'   \item{total_tooth}{Total number of teeth}
#'   \item{gender}{Gender of the patient (factor)}
#'   \item{diabetes}{Indicates whether the patient has diabetes mellitus (factor)}
#'   \item{tobacco_ever}{Indicates whether the patient has a history of tobacco use (factor)}
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
NULL



