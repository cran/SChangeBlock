#include <Rcpp.h>
using namespace Rcpp;

//' Gini's mean difference
//' 
//' Calculates Gini's mean difference of a given vector \code{x}.
//' 
//' @param x numeric vector.
//' 
//' @return A numeric value.
//' 
//' @examples
//' x <- rnorm(100)
//' GMD(x)
//' 
//' @export
// [[Rcpp::export]]
double GMD(NumericVector x)
{
  int i, j;
  int n = x.size();
  
  double sum = 0;
  
  for(j = 1; j < n; j++)
  {
    for(i = 0; i < j; i++)
    {
      sum += std::fabs(x[i] - x[j]);
    }
  }
  sum = sum * 2 / (n * (n - 1));
  return sum;
}
