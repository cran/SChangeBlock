#' Block test p-value
#'
#' Returns the p-value of a test statistic according to the block test for structural changes.
#' 
#' @param tn test statistic
#' @param fun Character string; one of "gmd" (default), "var", "jb", "grubbs", "ANOVA", "ad", "sw" (see \code{\link{block_stat}} for details).
#' 
#' @returns A numeric value between 0 and 1.
#' 
#' @export
#' @importFrom nortest ad.test
block_pValue <- function(tn, fun = "gmd")
{
  if(fun == "gmd" | fun == "var") return(1 - pnorm(tn))
  else if(fun == "jb") return(1 - pchisq(tn, 2))
  else if(fun == "grubbs") ## Grubbs: no real p-value
    return(as.numeric(tn < crit.grubbs(attr(tn, "n"), 0.05))) ### 0.05 ????
  else if(fun == "ANOVA") return(1 - pf(tn, attr(tn, "k") - 1, attr(tn, "N") - attr(tn, "k")))
  else if(fun == "ad") return(nortest::ad.test(tn)$p.value) 
  else if(fun == "sw") return(shapiro.test(tn)$p.value)
  
  stop("fun unknown")
}


#' Block test for structural changes
#' 
#' A test to detect whether an underlying time series or random field is stationary (hypothesis)
#' or if there is a location shift present in a region (alternative).
#'
#' @param x times series or random field to be tested. Either a numeric vector or a numeric matrix.
#' @param s parameter for the size of the blocks, 0.5 < s < 1, block length \eqn{l_n = [n^s]}. Default is \code{\link{sOpt}}\code{(n, 0.6)}.
#' @param fun Character string; one of "gmd" (default), "var", "jb", "grubbs", "ANOVA", "ad", "sw" (see \code{\link{block_stat}} for details).
#' @param varEstim variance estimator or variance estimation of the whole field or times series.
#'                 Either a function to estimate the variance with, or a numeric value.
#'
#' @return A list of the class "htest" containing the following components:
#'  \item{statistic}{value of the test statistic (numeric).}
#'  \item{p.value}{p-value (numeric).}
#'  \item{alternative}{alternative hypothesis (character string).}
#'  \item{method}{name of the performed test (character string).}
#'  \item{data.name}{name of the data (character string).}
#'
#' @seealso [block_stat], [block_pValue]
#' 
#' @references \enc{Görz}{Goerz}, S. and Fried, R. (2025). "Detecting changes in the mean of spatial random fields on a regular grid", \emph{arXiv preprint} \href{https://arxiv.org/abs/2512.11599}{arXiv:2512.11599}
#' @references Schmidt, S. K. (2024). Detecting changes in the trend function of heteroscedastic time series. Bernoulli, 30(4), 2598-2622.
#'
#' @examples
#' # time series with a shift 
#' x <- arima.sim(model = list(ar = 0.5), n = 100)
#' x[1:50] <- x[1:50] + 1
#' block_test(x, sOpt(100, 0.6))
#' 
#' # field without a shift and ordinary variance
#' X <- genField(c(50, 50))
#' block_test(X, sOpt(50, 0.6), "var")
#' 
#' # field with a shift and ordinary variance
#' X <- genField(c(50, 50), type = 2)
#' block_test(X, sOpt(50, 0.6), "var")
#' 
#' #' # GMD test statistic, scaling variance estimated by the mad
#' block_stat(X, 0.6, fun = "var", varEstim = mad)
#' 
#' @export
block_test <- function(x, s, fun = "gmd", varEstim = var)
{
  stat <- block_stat(x, s, fun, varEstim)
  names(stat) <- "S"
  pval <- block_pValue(stat, fun)
  
  Dataname <- deparse(substitute(x))
  
  erg <- list(alternative = "two-sided", method = paste(fun, "block test for structural changes"), 
              data.name = Dataname, statistic = stat, blocksize = attr(stat, "blocksize"),
              p.value = pval#, lrv = list(method = method, param = attr(stat, "param"), value = attr(stat, "sigma")), 
              )
  
  class(erg) <- "htest"
  return(erg)
}
