#' Mixture distribution
#' 
#' Generates a random sample from a mixture normal distribution.
#'
#' @param n sample size, numeric.
#' @param q probability that observation is drawn from the contamination distribution, numeric.
#' @param h mean of the contamination distribution, numeric.
#' @param sigma standard deviation of the contamination distribution, numeric.
#'
#' @details The resulting sample is drawn from the distribution
#' \deqn{(1 - q)\mathcal{N}(0, 1)\; + \; q \mathcal{N}(h, \sigma^2).}
#'
#' @return Numeric vector of length n containing the random sample.
#'
#' @examples
#' # random sample with 0.01 chance of contamination distribution with mean 10
#' rmix(100)
#' 
#' # random sample with 0.01 chance of contamination distribution with standard deviation 10
#' # IMPORTANT: h needs to be set to 0!
#' rmix(100, h = 0, sigma = 1)
#' 
#' @export
rmix <- function(n, q = 0.01, h = 10, sigma = 1)
{
  x <- rnorm(n)
  index <- sample(n, rbinom(1, n, q))
  x[index] <- sigma * x[index] + h
  return(x)
}



