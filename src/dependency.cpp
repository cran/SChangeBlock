#include <Rcpp.h>
using namespace Rcpp;

//' Dependency Matrix Theta
//' 
//' This function generates a symmetric dependency matrix \eqn{\Theta} of a 
//' specific type of spatial MA(q) model.
//' 
//' @param q model order (integer).
//' @param param MA parameter (numeric value between 0 and 1).
//' @param structure Character string, either "MA" or "AR" indicating the structure of the dependency matrix. Details below.
//' 
//' @return A matrix of size (2q + 1) x (2q + 1).
//' 
//' @details Symmetric spatial MA(q) model (or an approximation to a spatial AR(1) model) for 2-dim. random fields:
//'          \deqn{Y_{ij} = \sum_{k = -q}^q \sum_{l = -q}^q \theta_{kl} \varepsilon_{kl}.}
//'          \eqn{(\theta_{kl}) = \Theta}. \cr \cr
//'          For "MA": \deqn{\theta_{kl} = \code{param}^{|k - q - 1| + |l - q - 1|}.}
//'          For "AR": \deqn{\theta_{kl} = \tilde{\theta}_{kl} / \sqrt{\sum_{|k| \leq q} \sum_{|l| \leq q} \tilde{\theta}^2_{kl}} 
//'          \quad \text{ with } \quad \tilde{\theta}_{kl} = \code{param}^{\sqrt{k^2 + l^2}}.}
//'          
//' @examples
//' genTheta(1, 0.2, "MA")
//' 
//' genTheta(40, 0.2, "AR")
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix genTheta(int q, NumericVector param, String structure = "MA")
{
 if(q == 0) 
 {
   NumericMatrix P(1, 1);
   P(0, 0) = 1;
   return P;
 }
 
 for(int i = 0; i < param.length(); i++)
 {
   if(param[i] < 0 || param[i] > 1) stop("all values of param have to be between 0 and 1");
 }
 
 int p = 2 * q + 1;
 int i, j;
 NumericMatrix Theta(p, p);
 
 if(structure == "MA")
 {
   if(param.length() == 2 * q)
   {
     for(i = 0; i < p; i++)
     {
       for(j = 0; j < p; j++)
       {
         if(i != q || j != q) Theta(i, j) = param(abs(i - q) + abs(j - q) - 1);
       }
     }
     Theta(q, q) = 1;
   } else if(param.length() == 1)
   {
     for(i = 0; i < p; i++)
     {
       for(j = 0; j < p; j++)
       {
         Theta(i, j) = pow(param(0), abs(i - q) + abs(j - q));
       }
     }
   } else
   {
     stop("param has the wrong length! Must be either 1 or 2*q.");
   }
 } else if(structure == "AR")
 {
   double sumTheta = 0;
   
   if(param.length() == 1)
   {
     for(i = 0; i < p; i++)
     {
       for(j = 0; j < p; j++)
       {
         Theta(i, j) = pow(param(0), sqrt(pow(i - q, 2) + pow(j - q, 2)));
         sumTheta += pow(Theta(i, j), 2);
       }
     }
     sumTheta = sqrt(sumTheta);
   
     for(i = 0; i < p; i++)
     {
       for(j = 0; j < p; j++)
       {
         Theta(i, j) /= sumTheta;
       }
     }
   } else
   {
     stop("param has the wrong length! Must be 1.");
   }
 }

 return Theta;
}


// [[Rcpp::export]]
NumericMatrix dependency(NumericMatrix E, Nullable<NumericMatrix> Theta_ = R_NilValue, 
                         Nullable<IntegerVector> q_ = R_NilValue, Nullable<NumericVector> param_ = R_NilValue)
{
  int pn, pm;
  NumericMatrix Theta; 
  
  if(Theta_.isNull())
  {
    int q = (IntegerVector (q_))(0);
    pn = pm = 2 * q + 1;
    Theta = genTheta(q, (NumericVector) param_);
  } else
  {
    Theta = NumericMatrix (Theta_);
    pn = Theta.nrow();
    pm = Theta.ncol();
  }
  
  // Theta((pn - 1) / 2, (pm - 1) / 2) = 1;
  
  int n2q = E.nrow(); 
  int m2q = E.ncol();
  int n = n2q - pn + 1; 
  int m = m2q - pm + 1;
  
  NumericMatrix X(n, m);
  
  int i, j, i2, j2;
  
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < m; j++)
    {
      X(i, j) = 0;
      for(i2 = 0; i2 < pn; i2++)
      {
        for(j2 = 0; j2 < pm; j2++)
        {
          X(i, j) += Theta(i2, j2) * E[(i + i2) * m2q + j + j2];
        }
      }
    }
  }
  
  return X;
}
