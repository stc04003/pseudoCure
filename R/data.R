#' Dental data for illustration
#'
#' @description A data on the survival of teeth with many predictors
#' The variables are as follows:
#' \describe{
#'   \item{id}{subjects identification}
#'   \item{time}{tooth survival time subject to right censoring}
#'   \item{event}{Tooth loss status: 1 = lost, 0 = not lost}
#'   \item{mobil}{Mobility score (on a scale 0-5)}
#'   \item{bleed}{Bleeding on probing (percentage)}
#'   \item{plaque}{Plaque score (percentage)}
#'   \item{pocket_mean}{Periodontal probing depth (tooth-level mean, in mm)}
#'   \item{pocket_max}{Periodontal probing depth (tooth-level maximum, in mm)}
#'   \item{cal_mean}{Clinical attachment level (tooth-level mean, in mm)}
#'   \item{cal_max}{Clinical attachment level (tooth-level maximum, in mm)}
#'   \item{fgm_mean}{Free gingival margin (tooth-level mean, in mm)}
#'   \item{fgm_max}{Free gingival margin (tooth-level maximum, in mm)}
#'   \item{filled}{Number of filled surfaces on the tooth}
#'   \item{decay_new}{Number of newly decayed surfaces}
#'   \item{decay_recur}{Number of recurrently decayed surfaces}
#'   \item{dfs}{Decayed and filled surfaces (total count)}
#'   \item{crown}{Indicates whether the tooth has a crown (factor)}
#'   \item{endo}{Indicates whether the tooth has undergone endodontic therapy (factor)}
#'   \item{pontic}{Bridge pontic}
#'   \item{filled_tooth}{Indicates whether the tooth is filled (factor)}
#'   \item{decayed_tooth}{Indicates whether the tooth is decayed (factor)}
#'   \item{bleed_ave}{Mean bleeding on probing (percentage)}
#'   \item{plaque_ave}{Mean plaque index (percentage)}
#'   \item{pocket_mean_ave}{Mean periodontal probing depth (mean of tooth-level means, in mm)}
#'   \item{pocket_max_ave}{Mean periodontal probing depth (mean of tooth-level maximums, in mm)}
#'   \item{cal_mean_ave}{Mean clinical attachment level (mean of tooth-level means, in mm)}
#'   \item{cal_max_ave}{Mean clinical attachment level (mean of tooth-level maximums, in mm)}
#'   \item{fgm_mean_ave}{Mean free gingival margin (mean of tooth-level means, in mm)}
#'   \item{fgm_max_ave}{Mean free gingival margin (mean of tooth-level maximums, in mm)}
#'   \item{mg_ave}{Mean mucogingival defect}
#'   \item{filled_sum}{Total number of filled surfaces across all teeth}
#'   \item{filled_ave}{Mean number of filled surfaces per tooth}
#'   \item{decay_new_sum}{Total number of newly decayed surfaces across all teeth}
#'   \item{decay_new_ave}{Mean number of newly decayed surfaces per tooth}
#'   \item{decay_recur_sum}{Total number of recurrently decayed surfaces across all teeth}
#'   \item{decay_recur_ave}{Mean number of recurrently decayed surfaces per tooth}
#'   \item{dfs_sum}{Total number of decayed and filled surfaces across all teeth}
#'   \item{dfs_ave}{Mean number of decayed and filled surfaces per tooth}
#'   \item{filled_tooth_sum}{Total number of filled teeth}
#'   \item{filled_tooth_avg}{Percentage of filled teeth}
#'   \item{decayed_tooth_sum}{Total number of decayed teeth}
#'   \item{decayed_tooth_ave}{Percentage of decayed teeth}
#'   \item{missing_tooth_sum}{Total number of missing teeth}
#'   \item{missing_tooth_ave}{Percentage of missing teeth}
#'   \item{total_tooth}{Total number of teeth}
#'   \item{dft}{Number of Decayed and Filled Teeth}
#'   \item{gender}{Gender of the patient (factor)}
#'   \item{diabetes}{Indicates whether the patient has diabetes mellitus (factor)}
#'   \item{tobacco_ever}{Indicates whether the patient has a history of tobacco use (factor)}
#'   \item{molar}{Indicates whether the tooth is a molar (logical).}
#' }
#'
#' @details
#' The data is a subset of the original dataset included in the \code{MST} package
#' under the name \code{Teeth}.
#' This subset contains the time to the first tooth loss due to periodontal reasons.
#' 
#' @usage data(Teeth1)
#' @docType data
#' @name Teeth1
#' @rdname Teeth1
#' @format A data frame with 5336 rows and 51 variables.
#'
#' @references
#' Calhoun, Peter and Su, Xiaogang and Nunn, Martha and Fan, Juanjuan (2018) Constructing Multivariate Survival Trees: The MST Package for R. \emph{Journal of Statistical Software}, \bold{83}(12).
NULL



