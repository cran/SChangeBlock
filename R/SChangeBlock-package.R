## usethis namespace: start
#' @import Rcpp 
#' @import stats
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices rgb
#' @importFrom nortest lillie.test
#' @importFrom nortest ad.test
#' @useDynLib SChangeBlock
## usethis namespace: end
NULL

#' @keywords internal
#' 
#' @description
#' Provides methods to detect structural changes in time series or random fields (spatial data). Focus is on the detection of abrupt changes or trends in independent data, but the package also provides a function to de-correlate data with dependence. The functions are based on the test suggested in Schmidt (2024) <doi:10.3150/23-BEJ1686> and the work in Görz and Fried (2025+).
#' 
#' 
"_PACKAGE"
