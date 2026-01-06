#include <Rcpp.h>
using namespace Rcpp;

// //' @export
// // [[Rcpp::export]]
// NumericMatrix mtrx(RObject X)
// {
//   NumericMatrix x;
//   if(is<NumericVector>(X))
//   {
//     if(Rf_isMatrix(X))
//     {
//       x = (as<NumericMatrix>(X));
//     } else
//     {
//       NumericVector x2 = (as<NumericVector>(X));
//       int n = x2.size();
//       x2.attr("dim") = Dimension(n, 1);
//       x = as<NumericMatrix>(x2);
//     }      
//   } else
//   {
//     stop("Neither numeric vector nor numeric matrix!");
//   }
//   
//   
//   return x; 
// }

//' Mu
//' 
//' This function returns a vector of the block means for a given random field X.
//' 
//' @param x Numeric vector or matrix.
//' @param group strictly positive integer vector or matrix indicating the group (or block) of the corresponding observation in X.
//'        Overwrites `l` if specified.
//' @param l block length. Integer vector of length 1 or 2, depending on the number of dimensions of X, with strictly positive entries.
//' 
//' @return A numeric vector of length \code{floor(n[1] / l[1]) * floor(n[2] / l[2])}.
//' 
//' @examples 
//' X <- genField(c(50, 100), H = 100, type = 2)
//' M <- Mu(X, l = c(10, 20))
//' 
//' plot(X)
//' image(matrix(M, ncol = 5))
//' 
//' @export
// [[Rcpp::export]]
NumericVector Mu(RObject x, Rcpp::Nullable<Rcpp::RObject> group = R_NilValue, 
                 Rcpp::Nullable<Rcpp::IntegerVector> l = R_NilValue)
{
  NumericMatrix X;
  IntegerMatrix Group;
  NumericVector MuVec;
  int l_n, l_m;
  
  int n, m;
  
  if(l.isNotNull())
  {
    IntegerVector L(l);
    
    l_n = L[0];
    if(L.length() == 1)
    {
      l_m = l_n;
    } else
    {
      l_m = L[1];
    }
    
    if(l_n <= 0 || l_m <= 0) stop("Only positive values allowed for l");
  }
  
  if(is<NumericVector>(x))
  {
    if(Rf_isMatrix(x))
    {
      X = (as<NumericMatrix>(x));
      n = X.nrow();
      m = X.ncol();
      
      if(group.isNotNull()) Group = (as<IntegerMatrix>(group)); 
    } else
    {
      RObject y = clone(x);
      NumericVector x2 = (as<NumericVector>(y));
      n = x2.size();
      m = 1;
      x2.attr("dim") = Dimension(n, 1);
      X = as<NumericMatrix>(x2);
      
      if(group.isNotNull())
      {
        Rcpp::Nullable<Rcpp::RObject> g = clone(group);
        IntegerVector g2 = (as<IntegerVector>(g));
        g2.attr("dim") = Dimension(n, 1);
        Group = as<IntegerMatrix>(g2);
      } else if(l.isNotNull())
      {
        l_m = 1;
      } 
    }
  }
  
  // start with taking the means
  // first: by specified group
  if(group.isNotNull())
  {
    if(n != Group.nrow() || m != Group.ncol()) stop("'group' does not have the same dimensions as X.");
    
    int i, j;  
    
    double ngroups = max(Group);
    
    MuVec = NumericVector(ngroups);
    IntegerVector ln(ngroups);
    
    for(i = 0; i < ngroups; i++)
    {
      MuVec[i] = 0;
      ln[i] = 0;
    }
    
    for(i = 0; i < n; i++)
    {
      for(j = 0; j < m; j++)
      {
        MuVec[Group(i, j) - 1] += X(i, j);
        ln[Group(i, j) - 1]++;
      }
    }
    
    for(i = 0; i < ngroups; i++)
    {
      MuVec[i] /= ln[i];
    }
    // then: by specifying the block length(s) for same-sized blocks
  } else if(l.isNotNull())
  {
    int b_n = floor(X.nrow() / l_n);
    int b_m = floor(X.ncol() / l_m);
    
    MuVec = NumericVector(b_n * b_m);
    int i, j, i2, j2, blocklength;
    
    blocklength = l_n * l_m;
    
    for(i = 0; i < b_n; i++)
    {
      for(j = 0; j < b_m; j++)
      {
        MuVec(i * b_m + j) = 0;
        
        for(i2 = i * l_n; i2 < (i + 1) * l_n; i2++)
        {
          for(j2 = j * l_m; j2 < (j + 1) * l_m; j2++)
          {
            MuVec(i * b_m + j) += X(i2, j2);
          }
        }
        
        MuVec(i * b_m + j) /= blocklength;
      }
    }
  } else
  {
    stop("Either 'group' or 'l' needs to be specified");
  }
  
  return MuVec;
}
 // NumericVector Mu(RObject x, IntegerVector l)//, Nullable<IntegerVector> e = R_NilValue)
 // {
 //   int l_n, l_m; //ex, ey;
 //   
 //   l_n = l[0]; 
 //   if(l_n <= 0) stop("Only positive values allowed for l");
 //   
 //   NumericMatrix X;
 //   if(is<NumericVector>(x))
 //   {
 //     if(Rf_isMatrix(x))
 //     {
 //       X = (as<NumericMatrix>(x));
 //       if(l.length() == 1)
 //       {
 //         l_m = l_n;
 //       } else
 //       {
 //         l_m = l[1];
 //         if(l_m <= 0) stop("Only positive values allowed for l");
 //       }
 //     } else
 //     {
 //       RObject y = clone(x);
 //       NumericVector x2 = (as<NumericVector>(y));
 //       int n = x2.size();
 //       x2.attr("dim") = Dimension(n, 1);
 //       X = as<NumericMatrix>(x2);
 //       l_m = 1;
 //     }      
 //   } else
 //   {
 //     stop("Neither numeric vector nor numeric matrix!");
 //   }
 //   
 //   // if(e.isNull())
 //   // {
 //   //   ex = ey = 0;
 //   // } else
 //   // {
 //   //   IntegerVector e_(e);
 //   //   ex = e_(0);
 //   //   
 //   //   if(e_.length() == 1)
 //   //   {
 //   //     ey = 0; 
 //   //   } else
 //   //   {
 //   //     ey = e_(1);
 //   //   }
 //   // }
 //   
 //   int b_n = floor(X.nrow() / l_n);
 //   int b_m = floor(X.ncol() / l_m);
 //   
 //   NumericVector MuVec(b_n * b_m);
 //   int i, j, i2, j2, blocklength;
 //   
 //   blocklength = l_n * l_m;
 //   
 //   // MuVec(0) = 0;
 //   // count = 0;
 //   // for(i2 = ex; i2 < l_n; i2++)
 //   // {
 //   //   for(j2 = ey; j2 < l_m; j2++)
 //   //   {
 //   //     MuVec(0) += X(i2, j2);
 //   //     count++;
 //   //   }
 //   // }
 //   // MuVec(0) /= count;
 //   
 //   for(i = 0; i < b_n; i++)
 //   {
 //     for(j = 0; j < b_m; j++)
 //     {
 //       MuVec(i * b_m + j) = 0;
 //       
 //       for(i2 = i * l_n; i2 < (i + 1) * l_n; i2++)
 //       {
 //         for(j2 = j * l_m; j2 < (j + 1) * l_m; j2++)
 //         {
 //           MuVec(i * b_m + j) += X(i2, j2);
 //         }
 //       }
 //       
 //       MuVec(i * b_m + j) /= blocklength;
 //     }
 //   }
 //   
 //   return MuVec;
 // }
