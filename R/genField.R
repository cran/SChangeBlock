#' Generate random fields
#'
#'@description
#'Generate random fields on a regular grid. Dependency can be modeled to some extend.
#'
#'
#' @param n dimensions of the requested random field. Numeric vector.
#' @param distr function specifying the error distribution.
#' @param type change region type (integer or character). See \code{\link{changeRegion}}.
#' @param H height of the location shift (numeric value).
#' @param Theta matrix specifying the dependency between observations of a 2-dim random field. 
#'            Has to be a 2-dim. matrix where there is an odd number of entries in both dimensions. 
#'            Explanations are given in [genTheta].
#' @param q dependency parameter for the model order of the MA field (integer > 0). 
#' @param param dependency parameters for the model parameters of the MA field (numeric vector.) 
#'              Both \code{q} and \code{param} are ignored if \code{Theta} is supplied.
#'              A dependency matrix is generated using \code{\link{genTheta}}.
#' @param ... additional arguments for the generation of \code{\link{changeRegion}}.
#'
#' @details
#' The dependent random field is generated as follows: 
#' Denote with \eqn{(e_{ij})} the error matrix from which the random field \eqn{(Y_{ij})}
#' (under the hypothesis, i.e. without any location shift) is generated. 
#'
#' @returns
#' The function returns a matrix of dimension \code{n} that has the \code{\link{class}} "RandomField".
#'
#' @examples
#' genField(c(50, 50))
#' genField(c(50, 50), type = 3, Theta = genTheta(1, 0.2))
#' 
#' @seealso [changeRegion], [genTheta]
#'
#' @export
genField <- function(n, distr = rnorm, type = 0L, H = 100, 
                     Theta = NULL, q = NULL, param = NULL, ...) 
{
  ## is there supposed to be dependency? If yes, the error matrix has to be
  ## larger than the requested field:
  if(!is.null(Theta))
  {
    pn <- dim(Theta)
    if(any(pn %% 2 == 0)) stop("Dimensions of Theta are incorrect!")
    pn2 <- n + pn - 1
  } else if(!is.null(q) && !is.null(param))
  {
    pn2 <- n + 2 * q
  } else
  {
    pn2 <- n
  }
  
  # error matrix
  X <- matrix(distr(prod(pn2)), nrow = pn2[1])
  
  if(!is.null(Theta) || !(is.null(q) || is.null(param))) 
  {
    # multiply error matrix X with dependency matrix Theta:
    X <- dependency(X, Theta, q, param)
  } 
  
  ## add change region:
  if(type != 3L && type != 0L)
  {
    alt <- numeric(prod(n))
    altIndex <- changeRegion(n, type = type, ...)
    alt[altIndex] <- H / prod(n)^(0.5)# / (2 * length(altIndex))
    X <- X + alt
  } else if(type == 3L)
  {
    X <- X + changeRegion(n, type = 3L) * H
  }
  
  class(X) <- "RandomField"
  return(X)
}


#' Change Region
#' 
#' Generates the indices for different types of change regions.
#'
#' @param n dimensions of the random field. Numeric vector.
#' @param s parameter for the size of the blocks, 0.5 < s < 1, block length \eqn{l_n = [n^s]}.
#' @param type change region type (integer or character). See "Details".
#' @param middle,delta,distFun parameters for type 4L. See "Details".
#' 
#' @return 
#' \item{Types 1a, 1b, 1c, 2, 4, 5}{A vector of indices.}
#' \item{Type 3}{A numeric (n x n) matrix containing the heights of the shifts at the corresponding locations.}
#' 
#' @details
#' Change region types: 
#' * 1L or "1a": exactly one block is shifted
#' * "1b": size of the shift region is one block, but the shift region lies in two blocks
#' * "1c": size of the shift region is one block, but the shift region lies in four blocks
#' * 2L: exactly half of the data is shifted
#' * 3L: there is a steady increase from left to right
#' * 4L: 
#' * 5L: for demonstration purposes: the field is divided into 10 "columns". Every other column is shifted
#' * 6L: 
#' * 7L: a "circle" including everything within \code{distFun} \code{delta} from \code{middle} is shifted
#' 
#' @seealso [genField]
#' 
#' @examples
#' changeRegion(c(50, 50), 0.6, "1a")
#' changeRegion(c(50, 50), type = 2)
#' changeRegion(c(50, 50), type = 3L)
#' changeRegion(c(50, 50), type = 7L, middle = c(10, 10))
#' 
#' @export
changeRegion <- function(n, s, type = 1L, middle, delta = 0.15, distFun = dist)
{
  if(missing(s)) s <- sOpt(n, 0.6)
  ln <- round(n^s)

  if(is.na(n[2]))
  {
    n[2] <- 1
    ln[2] <- 1
  }

  alt <- numeric(prod(n))

  if(type == "1a" || type == 1L)
  {
    # one whole block
    index <- sapply(seq(1, ln[1] * n[2], n[1]), function(x) x:(x + ln[1] - 1))
  } else if(type == "1b")
  {
    # shift lies in two blocks
    index <- sapply(seq(1, ln[1] * n[2], n[1]),
                    function(x) (x + floor(ln[1] / 2)):(x + floor(1.5 * ln[1]) - 1))
  } else if(type == "1c")
  {
    # shift lies in four blocks
    index <- sapply(seq(floor(ln[1] / 2) * n[2], floor(1.5 * ln[1]) * n[2] - 1, n[1]),
                    function(x) (x + floor(ln[1] / 2) + 1):(x + floor(1.5 * ln[1])))
  } else if(type == 2L)
  {
    # half of the data
    index <- 1:(prod(n) / 2)
  } else if(type == 7L)
  {
    # "circle" including everything within distFun delta from middle

    if(missing(middle)) middle <- floor(n / 2)
    delta <- delta * distFun(rbind(c(1, 1), n))
    
    index <- do.call(rbind, sapply(1:n[1], function(x) 
    {
      res <- which(Vectorize(function(y) distFun(rbind(c(x, y), middle)))(1:n[2]) < delta)
      if(length(res) > 0) return(data.frame(x, res))
      return(NULL)
    }))
    index <- (index[, 2] - 1) * n[1] + index[, 1]
  } else if(type == 3L)
  {
    # steady increase from left to right
    return(matrix(rep(1 / sqrt(prod(n)) * (0:(n[2]-1)) / (n[2] - 1), each = n[1]), nrow = n[1]))
  } else if(type == 5L)
  {
    # field is divided into 10 "columns", every other column is shifted
    i <- floor(n[1] / 10)
    index <- c((n[2] * i + 1):(n[2] * 2 * i), (n[2] * 3 * i + 1):(n[2] * 4 * i), 
               (n[2] * 5 * i + 1):(n[2] * 6 * i), (n[2] * 7 * i + 1):(n[2] * 8 * i), 
               (n[2] * 9 * i + 1):(n[2] * 10 * i))
  } else if(type == 6L)
  {
    upper <- min(n)
    index <- unlist(sapply(1:upper, function(i) ((i - 1) * n[1] + 1):((i - 1) * n[1] + i)))
  } else if(type == 4L)
  {
    i1 <- floor(n / 4)# + as.integer((n) %% 4 > 0)
    i2 <- floor(n / 2)
    # i2 <- floor(n / 2) + as.integer((n) %% 4 > 0)
    
    index <- c(sapply(seq(i1[2] * n[1] + 1, i2[2] * n[1], n[1]), function(x) x:(x + i2[1] - 1)), 
               # rectangle
               sapply(seq(i1[1] + 1, i1[2] * n[1], n[1]), function(x) x:(x + (i2[1] - i1[1]) - 1)))
               # square
  }

  # alt[index] <- H / prod(n)^(0.43)
  return(as.integer(index))
}
