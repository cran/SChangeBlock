############# HELPERS ##########################################################

# #' Object dimensions
# #'
# #' Function to return the dimensions of both a vector and other objects
# #'     (like matrices and data.frames)
# #'
# #'
# dim2 <- function(x)
# {
#   if(is.vector(x)) return(c(length(x), 0))
#   return(dim(x))
# }

# RsSizes <- function(n, lower = 0.5, upper = 1, step = 0.1)
# {
#   s <- seq(lower, upper, step)
#   ln <- round(n^s)
#   bn <- floor(n / ln)
# 
#   data.frame(n = n, s = s, ln = ln, bn = bn, "ln*bn" = ln * bn, diff = n - (ln * bn))
# }
# 
# RsOpt <- function(n, s = 0.6)
# {
#   sapply(n, function(n)
#   {
#     res <- RsSizes(n, 0.5, 1, 0.001)
#     res <- res[res$diff == 0, ]
#     res$s[which.min(abs(res$s - s))]
#   })
# }

mu <- function(x, ln, fun = "mean")
{
  if(!is.function(fun)) fun <- eval(parse(text = fun))
  if(is.vector(x)) x <- matrix(x)
  if(is.na(ln[2])) ln[2] <- 1
  n <- dim(x)
  
  res <- sapply(seq(ln[1] + 1, n[1] + 1, ln[1]), function(i)
  {
    sapply(seq(ln[2] + 1, n[2] + 1, ln[2]), function(j)
    {
      fun(x[(i - ln[1]):(i - 1), (j - ln[2]):(j - 1)])
    })
  })
  return(res)
}

#' Skewness and kurtosis
#' 
#' Compute the skewness and the kurtosis of the data vector \code{x}.
#' 
#' @param x numeric vector.
#' 
#' @returns A numeric value.
#' 
#' @examples 
#' skewness(rnorm(100))
#' skewness(rexp(100, 2))
#' 
#' kurtosis(rnorm(100))
#' kurtosis(rt(100, 5))
#' @rdname higherMoments
#'@export
skewness <- function(x)
{
  n <- length(x)
  sigma <- var(x) * (n - 1) / n
  return(mean((x - mean(x))^3) / sigma^(3/2))
}

#' @rdname higherMoments
#'
#'@export
kurtosis <- function(x)
{
  n <- length(x)
  sigma <- var(x) * (n - 1) / n
  return(mean((x - mean(x))^4) / sigma^2)
}

#' Grubbs outlier test
#' 
#' Computes the test statistic and the critical value of the outlier test according to Grubbs.
#' 
#' @param x numeric vector.
#' @param varEstim Variance estimation or estimation function. Either a numeric value or a function taking one argument.
#' 
#' @return numeric value.
#' 
#' @examples
#' x <- rnorm(100)
#' grubbs(x) > crit.grubbs(100, 0.05)
#' 
#' # add outlier
#' x[1] <- x[1] + 100
#' grubbs(x) > crit.grubbs(100, 0.05)
#' 
#' 
#' @rdname grubbs
#' @export
grubbs <- function(x, varEstim = var)
{
  if(is.function(varEstim)) varEstim <- varEstim(x)
  return(max(abs(x - mean(x))) / sqrt(varEstim))
}


#' Critical value of the Grubbs test
#' 
#' @param n sample size; positive numeric value.
#' @param alpha significance level; between 0 and 1.
#' 
#' @rdname grubbs
#'@export
crit.grubbs <- function(n, alpha = 0.05)
{
  (n - 1) / sqrt(n) * sqrt(qt(alpha / 2 / n, n-2)^2 /
                             (n - 2 + qt(alpha / 2 / n, n-2)^2))
}
