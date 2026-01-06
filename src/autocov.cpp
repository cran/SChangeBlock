#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// double kTH(double x)
// {
//   if(fabs(x) <= 1) return (1 + cos(M_PI * x)) / 2;
//   return 0;
// }

double kFlatTop(double x)
{
  if(fabs(x) <= 0.5) return 1; else if(fabs(x) <= 1) return 2 - 2 * fabs(x);
  return 0;
}

// [[Rcpp::export]]
double gamma(NumericMatrix X, int h1, int h2)
{
  int n = X.nrow();
  int m = X.ncol();
  
  if(h1 < 0 && h2 < 0) 
  {
    h1 = -h1;
    h2 = -h2;
  }
  
  int i, j;
  double res = 0; 
  
  for(i = std::max(0, -h1); i < std::min(n - h1, n); i++)
  {
    for(j = std::max(0, -h2); j < std::min(m - h2, m); j++)
    {
      res += (X(i, j) * X(i + h1, j + h2));
    }
  }
    
  return res / (n * m);
}

// [[Rcpp::export]]
double gammaDiff(NumericMatrix X, int h1, int h2)
{
  int n = X.nrow();
  int m = X.ncol();
  
  if(h1 < 0 && h2 < 0) 
  {
    h1 = -h1;
    h2 = -h2;
  }
  
  int i, j;
  double res = 0; 
  
  for(i = std::max(0, -h1); i < std::min(n - h1, n); i++)
  {
    for(j = std::max(0, -h2); j < std::min(m - h2, m); j++)
    {
      res += pow(X(i, j) - X(i + h1, j + h2), 2);
    }
  }
 
 return res / (2 * (n - h1) * (m - h2));
}


//' Autocovariance matrix
//' 
//' Estimates the autocovariance matrix for a given data matrix X. Via the parameter \code{direction}, it is possible to estimate only 
//' row- or columnwise autocovariance matrices, which is useful if the autocovariance function is separable.
//' 
//' @param X numeric matrix,
//' @param b numeric vector containing two integer values: the bandwidths for the row- resp. column-wise estimation.
//'        (Up to which lag should the autocovariances be estimated?) If \code{direction} > 0: only one integer must be supplied. 
//'        Bandwidths must be smaller than the dimensions of X.
//' @param M numeric vector containing two integer values, only needed for \code{type = 1}, see Details.
//' @param direction 0: all directions, 1: only row-wise autocovariances, 2: only column-wise autocovariances.
//' @param type 0: ordinary autocovariance estimation, 1: difference-based autocovariance estimation. See Details.
//' 
//' @returns A numeric matrix of size \eqn{N \times N}. If \code{direction = 0} then \eqn{N =} \code{prod(dim(X))}. 
//'          If \code{direction = 1} then \eqn{N =} \code{ncol(X)}, if \code{direction = 2} then \eqn{N =} \code{nrow(X)}.
//'          
//' @details In this function, the autocovariance matrix of \code{X} is interpreted as the autocovariance matrix of 
//'          \code{x = as.vector(X)}, i.e. where \code{X} is ordered into a vector column-wise. 
//'          If \code{type = 0}, the autocovariance to lags \eqn{h_1, h_2} is estimated using the regular estimator
//'          \deqn{\hat{\gamma}_{\text{reg}}(h_1, h_2) = 
//'           \frac{1}{(n-h_1)(m-h_2)}\sum_{i = 1}^{n-h_1} \sum_{j = 1}^{m-h_2} (Y_{i, j} - \bar{Y})(Y_{i+h_1, j+h_2} - \bar{Y}).}
//'          If \code{type = 1}, the autocovariance to lags \eqn{h_1, h_2} is estimated by a difference-based version, inspired by 
//'          the estimator of Tecuapetla-Gómez and Munk (2017) for time series:
//'          \deqn{\hat{\gamma}_{\text{diff}}(h_1, h_2) = \hat{\sigma}_{\text{diff}}^2 -
//'           \frac{1}{2(n - h_1)(m - h_2)}\sum_{i = 1}^{n-h_1} \sum_{j = 1}^{m-h_2} (Y_{i, j} - Y_{i + h_1, j + h_2})^2}
//'          with
//'          \deqn{\hat{\sigma}_{\text{diff}}^2 = \frac{1}{4} \left(\frac{1}{n(m - M_2)}\sum_{i = 1}^n \sum_{j = 1}^{m - M_2} 
//'          (Y_{i, j} - Y_{i, j+M_2})^2 + \frac{1}{(n - M_1)m}\sum_{i = 1}^{n-M_1} \sum_{j = 1}^{m} (Y_{i, j} - Y_{i+M_1, j})^2 \right),}
//'          where \eqn{M_1 =} \code{M[1]}, \eqn{M_2 = } \code{M[2]}.
//' 
//' @examples
//' X <- genField(c(20, 20))
//' autocov(X, c(4, 4))[1:10, 1:100]
//' 
//' # if separable:
//' Sigma1 <- autocov(X, 4, direction = 1)
//' Sigma2 <- autocov(X, 4, direction = 2)
//' kronecker(Sigma1, Sigma2)[1:10, 1:100]
//' 
//' @references
//' Tecuapetla‐Gómez, I., & Munk, A. (2017). Autocovariance estimation in regression with a discontinuous signal and m‐dependent errors: A difference‐based approach. Scandinavian Journal of Statistics, 44(2), 346-368.
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix autocov(NumericMatrix X, NumericVector b, IntegerVector M = IntegerVector::create(1, 1), int direction = 0, int type = 0)
{
 NumericMatrix Y = clone(X);
 Y = X - mean(X);

 int n = X.nrow();
 int m = X.ncol(); 
 
 IntegerVector nm = {n, m};
 // NumericVector S = sOpt(nm, s);
 
 int b_n, b_m;
 int M_n, M_m;
 
 if(direction == 1) 
 {
   n = 1;
   b_n = 1;
   b_m = b[0] + 1;
   M_n = 0; 
   M_m = M[0];
 } else if(direction == 2)
 {
   m = 1;
   b_n = b[0] + 1; 
   b_m = 1;
   M_n = M[0];
   M_m = 0;
 } else
 {
   b_n = b[0] + 1;
   b_m = b[1] + 1;
   M_n = M[0];
   M_m = M[1];
 }
 
 int N = n * m; 
 
 int i, j; 
 double thresh = R_NegInf;
 
 if(b_n > n + 1 || b_m > m + 1)
 {
   stop("Bandwidths cannot be larger than X.");
 }
 
 NumericMatrix A(N, N); 
 double cov; 
 int h1, h2;
 
 // main diagonal: 
 double var;
 if(type == 0) 
 {
   var = gamma(Y, 0, 0);
 } else
 {
   var = (gammaDiff(Y, M_n, 0) + gammaDiff(Y, 0, M_m)) / 2;
   // var = gammaDiff(Y, M_n, M_m); 
   // Rprintf("var = %f\n", var);
   // for(i = 0; i < M_n; i++)
   // {
   //   var += gammaDiff(Y, i, M_m);
   // }
   // Rprintf("var = %f\n", var);
   // for(j = 0; j < M_m; j++)
   // {
   //   var += gammaDiff(Y, M_n, j);
   // }
   // Rprintf("var = %f\n", var);
   // 
   // var = var / (M_n + M_m + 1);
   // Rprintf("var = %f\n", var);
 }

 for(i = 0; i < N; i++)
 {
   A(i, i) = var;
 }

 if(b_n > 1 || b_m > 1)
 {
   // x direction
   for(h1 = 1; h1 < b_n; h1++)
   {
     if(type == 0) cov = gamma(Y, h1, 0);
      else cov = (var - gammaDiff(Y, h1, 0));
     
     if(cov / var >= thresh) 
     {
       for(i = 0; i < N - h1; i++)
       {
         if(i % n < n - h1)
         {
           A(i + h1, i) = cov; 
           A(i, i + h1) = cov;
         }
       }
     }
   }
   
   // y direction
   for(h2 = 1; h2 < b_m; h2++)
   {
     if(type == 0) cov = gamma(Y, 0, h2);
      else cov = (var - gammaDiff(Y, 0, h2));
     
     if(cov / var >= thresh) 
     {
       for(i = 0; i < N - h2 * n; i++)
       {
         A(i + h2 * n, i) = cov; 
         A(i, i + h2 * n) = cov;
       }
     }
     // x and y; only x goes into the negative
     for(h1 = -b_n + 1; h1 < b_n; h1++)
     {
       if(h1 != 0)
       {
         if(type == 0) cov = gamma(Y, h1, h2);
          else cov = (var - gammaDiff(Y, h1, h2));
         
         if(cov / var >= thresh) 
         {
           for(j = std::max(0, -h1); j < N - h2 * n - std::max(0, h1); j++)
           {
             if((j - std::max(0, -h1)) % n < n - abs(h1))
             {
               A(j + h1 + h2 * n, j) = cov;
               A(j, j + h1 + h2 * n) = cov;
             }
           }
         }
       }
     }
   }
 }
 
 return A;
}

// double gammaWrapper(NumericMatrix X, int h1, int h2, NumericVector s)
// {
//   int n = X.nrow();
//   int m = X.ncol();
//   int i, j;
//   int index = 0; 
//   
//   int bn = round(pow(n, s[0]));
//   int bm = round(pow(m, s[1]));
//   
//   NumericVector Res( floor(n / bn) * floor(m / bm) );
//   
//   for(i = 0; i < n - bn + 1; i = i + bn)
//   {
//     for(j = 0; j < m - bm + 1; j = j + bm)
//     {
//       Res[index] = gamma(X(Range(i, i + bn - 1), Range(j, j + bm - 1)), h1, h2);
//       index++;
//     }
//   }
//   
//   return median(Res);
// }

// //' Long run variance
// //' 
// //' Estimates the long run variance of a 2-dimensional matrix \code{X} using kernel
// //' density estimation and the Tukey-Hanning kernel function.
// //' 
// //' @param X numeric matrix,
// //' @param b numeric vector containing exactly two values: the bandwidths for the row- reps. column-wise estimation.
// //' 
// //' @details ??????? Tukey-Hanning kernel ????????
// //' 
// //' @return A numeric value.
// //' 
// //' @examples
// //' X1 <- genField(c(50, 50), Phi = genPhi(1, 0.4))
// //' b <- bandwidth(X1, 1/3, 2/3)
// //' lrv(X1, b)
// //' 
// //' Phi <- matrix(c(0.08, 0.1, 0.08, 0.8, 1, 0.8, 0.08, 0.1, 0.08), ncol = 3)
// //' X2 <- genField(c(50, 50), Phi = Phi)
// //' b <- bandwidth(X2, 1/3, 2/3)
// //' lrv(X2, b)
// //' 
// //' @export
// // [[Rcpp::export]]
// double lrv(NumericMatrix X, NumericVector b = 1)
// {
//   NumericMatrix Y = clone(X);
//   Y = X - mean(X);
//   
//   int b_n, b_m; 
//   
//   b_n = b[0]; 
//   if(b.length() == 2)
//   {
//     b_m = b[1];
//   } else
//   {
//     stop("l has to contain exactly two values!");
//   }
//   
//   double rows = 0, cols = 0, both = 0;
//   int h1, h2;
//   
//   if(b_n > 1 || b_m > 1)
//   {
//     for(h1 = 1; h1 < b_n; h1++)
//     {
//       rows += gamma(Y, h1, 0) * kTH(h1 / b_n);
//       
//       for(h2 = 1; h2 < b_m; h2++)
//       {
//         both += gamma(Y, h1, h2) * kTH(h1 / b_n) * kTH(h2 / b_n);
//       }
//     }
//     for(h2 = 1; h2 < b_m; h2++)
//     {
//       cols += gamma(Y, 0, h2) * kTH(h2 / b_m);
//     }
//   }
//   
//   double v = gamma(X, 0, 0);
//   double res = v + 2 * rows + 2 * cols + 4 * both;
//   if(res <= 0) return v;
//   return res;
// }
