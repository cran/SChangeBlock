#'Bandwidth estimation 
#'
#'Calculate MSE-optimal bandwidths according to Andrews (1991).
#'
#'@param X numeric vector or matrix.
#'@param p1,p2 exponents for sample size n resp. estimated dependency, between 0 and 1.
#'@param lag lag to which the autocorrelations are to be estimated. Integer > 0 but smaller than the length resp. number of rows and columns of X.
#'
#'@details Bandwidth \eqn{\boldsymbol{b}^{(n,m)} = (b_1^{(n)}, b_2^{(m)})} is estimated via
#'         \deqn{b_i^{(k)} = \min\left(k-1, \max\left(1, k^{0.3} \left(\frac{2\rho_i}{1 - \rho_i^2}\right)^{0.3}\right)\right),} 
#'         where \eqn{\rho_1} and \eqn{\rho_2} are the mean row- resp. column-wise Spearman autocorrelations to lag 1.
#'
#'@returns A numeric vector containing one or two elements, depending on if a vector or matrix is supplied.
#'         In case of a matrix: the first value is the bandwidth for the row-wise 
#'         and the second one for the column-wise estimation.
#'         
#'@references Andrews, D. W. (1991). “Heteroskedasticity and autocorrelation consistent covariance 
#'            matrix estimation”. In: Econometrica: Journal of the Econometric Society, pp. 817–858.
#'
#'@examples 
#'X1 <- genField(c(50, 50), Theta = genTheta(1, 0.4))
#'bandwidth(X1, 0.3, 0.3)
#'
#'Theta <- matrix(c(0.08, 0.1, 0.08, 0.8, 1, 0.8, 0.08, 0.1, 0.08), ncol = 3)
#'X2 <- genField(c(50, 50), Theta = Theta)
#'bandwidth(X2, 1/3, 2/3)
#'
#' @export
bandwidth <- function(X, p1 = 0.3, p2 = 0.3, lag = 1)
{
  if(is.matrix(X))
  {
    n <- nrow(X)
    m <- ncol(X)
    
    # ln <- round(n^s)
    # lm <- round(m^s)

    if(lag < 1 | lag >= n | lag >= m) stop("Wrong lag supplied!")
    
    # Res1 <- numeric(floor(n / ln) * floor(m / lm))
    # Res2 <- Res1
    # 
    # index <- 1
    # 
    # sapply(seq(1, n - ln + 1, ln), function(i)
    # {
    #   t(sapply(seq(1, m - lm + 1, lm), function(j)
    #   {
    #     Res1[index] <<- mean(sapply(i:(i+ln-1), function(l) cor(X[l, j:(j+lm-2)], X[l, (j+1):(j+lm-1)], method = "spearman")))
    #     Res2[index] <<- mean(sapply(j:(j+lm-1), function(k) cor(X[i:(i+ln-2), k], X[(i+1):(i+ln-1), k], method = "spearman")))
    #     
    #     index <<- index + 1
    #   }))
    # })
    # 
    # rho1 <- abs(median(Res1))
    # rho2 <- abs(median(Res2))

    rho1 <- abs(mean(sapply(1:n, function(i) cor(X[i, 1:(m-lag)], X[i, -c(1:lag)], method = "spearman"))))
    rho2 <- abs(mean(sapply(1:m, function(i) cor(X[1:(n-lag), i], X[-c(1:lag), i], method = "spearman"))))

    param1 <- min(max(round(m^(p1) * ((2 * rho1) / (1 - rho1^2))^(p2)), 0), m-1)
    param2 <- min(max(round(n^(p1) * ((2 * rho2) / (1 - rho2^2))^(p2)), 0), n-1)
    
    if(is.na(param1)) param1 <- 1
    if(is.na(param2)) param2 <- 1
    
    return(c(param1, param2))
  } else
  {
    n <- length(X)
    if(lag < 1 | lag >= n) stop("Wrong lag supplied!")
    
    rho <- abs(cor(X[1:(n-lag)], X[(1+lag):n], method = "spearman"))
    return(min(max(round(n^p1 * (rho / (1 - rho))^p2), 1), n-1))
  }
}


#' De-correlation
#' 
#' De-correlates a random field or time series, so that the resulting values can be treated as independent.
#' 
#' @param X Random Field, numeric matrix, or time series
#' @param lags numeric vector containing two integer values: the bandwidths for the row- reps. column-wise autocovariance estimation.
#'             (Up to which lag should the autocovariances be estimated?) Lags must be smaller than the dimensions of X.
#' @param method 1L: square root of the matrix via [robcp::modifChol()], inversion via [solve()] \cr
#'               2L: square root and inversion via singular value decomposition \cr
#'               3L: square root via [expm::sqrtm()], inversion via [solve()]
#' @param separable if the autocovariance function is (assumed to be) separable in the two directions of X, 
#'                  those two autocovariances can be estimated separately and then combined (after square root and inversion) as a Kronecker product.
#' @param M numeric vector containing two integer values, only needed for \code{type = 1}, see [autocov()].
#' @param type 0: ordinary autocovariance estimation, 1: difference-based autocovariance estimation. See [autocov()].
#'                  
#' @details
#' The contents of \code{X} are ordered into a vector \eqn{x} column-wise. The autocovariance matrix \eqn{\Sigma} of \eqn{x} is estimated by [autocov()]. 
#' \eqn{\Sigma} is taken the square root of and being inverted using the functions specified in \code{method}. Then
#' \deqn{y = \Sigma^{-\frac{1}{2}} (x - \bar{x}).}
#' Then \eqn{y} is ordered back into a matrix \eqn{Y} with the same dimension as \eqn{X}.
#' 
#'                  
#' @returns De-correlated random field or time series; same data type and size as input \code{X}.
#'
#' @examples
#' x <- arima.sim(list(ar = 0.4), 200)
#' y <- decorr(x, 3)
#' \donttest{
#' oldpar <- par(mfrow = c(2, 2))
#' acf(x)
#' pacf(x)
#' acf(y)
#' pacf(y)
#' par(oldpar) }
#' 
#' X <- genField(c(20, 20), Theta = genTheta(1, 0.4))
#' Y <- decorr(X, c(2, 2))
#' 
#'
#' @importFrom robcp modifChol
#' @importFrom expm sqrtm
#' @export
decorr <- function(X, lags, method = 1L, separable = FALSE, M = 1, type = 0)
{
  # if(is.ts(X)) 
  # {
  #   tspX <- tsp(X)
  #   clsX <- NULL
  # } else
  # { 
  #   tspX <- NULL
  #   clsX <- class(X)
  # }
  
  attrX <- attributes(X)
  
  if(is.vector(X) | is.null(dim(X)))
  {
    lags[2] <- 0
    M[2] <- 0
    X <- as.matrix(X)
  } else 
  {
    if(length(lags) == 1) lags[2] <- lags[1]
    if(length(M) == 1) M[2] <- M[1]
  }
  
  x <- as.vector(X)
  X <- as.matrix(X)

  if(separable)
  {
    A1 <- autocov(X, c(1, lags[2]), direction = 1, M = M, type = type)
    A2 <- autocov(X, c(lags[1], 1), direction = 2, M = M, type = type)
    
    A.invsqrt <- kronecker(invsqrt(A1, method), invsqrt(A2, method))
  } else 
  {
    # browser()
    A <- autocov(X, lags, M = M, type = type)
    A.invsqrt <- invsqrt(A, method)
  }

  y <- t(A.invsqrt) %*% (x - mean(x))
  if(ncol(X) > 1) Y <- matrix(y, ncol = ncol(X)) else Y <- as.vector(y)
  
  # if(!is.null(tspX)) 
  # {
  #   Y <- ts(Y)
  #   tsp(Y) <- tspX
  # } else
  # {
  #   class(Y) <- clsX
  # }
  
  attributes(Y) <- attrX
  
  return(Y)
}

invsqrt <- function(M, method = 1L)
{
  if(method == 1L)
  {
    Mchol <- tryCatch(chol(M), error = function(e) e)
    if("error" %in% class(Mchol))
    {
      Mchol <- modifChol(M)
    }
    
    res <- solve(Mchol)
    
    return(res)
  } else if(method == 2L)
  {
    res <- svd(M)
    return(res$u %*% diag(1 / sqrt(res$d)) %*% t(res$v))
  } else if(method == 3L)
  {
    Mchol <- tryCatch(chol(M), error = function(e) e)
    if("error" %in% class(Mchol))
    {
      Mchol <- sqrtm(M)
    }
    return(solve(Mchol))
  }
  {
    stop("Wrong method")
  }
}
