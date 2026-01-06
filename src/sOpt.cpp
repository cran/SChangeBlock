#include <Rcpp.h>
using namespace Rcpp;


//' @param n integer value.
//' @param lower,upper lower and upper search border, between 0 and 1.
//' @param step size of the step for the search, between 0 and 1.
//' 
//' 
//' @return \code{sSizes} returns a data frame containing
//' \item{n}{the given sample size}
//' \item{s}{the exponent in question}
//' \item{ln}{the resulting block length}
//' \item{bn}{the corresponding number of block}
//' \item{ln.bn}{block length times number of blocks}
//' \item{diff}{difference between the given sample size and the number of observations covered by the blocks}
//' 
//' @examples
//' sSizes(50)
//' sSizes(50, 0.6, 0.8, 0.01)
//' 
//' @rdname sOpt
//' @export
// [[Rcpp::export]]
DataFrame sSizes(int n, double lower = 0.5, double upper = 1, double step = 0.1)
 {
   if(lower > upper) stop("invalid search interval");
   int N = floor((upper - lower) / step) + 1;
   NumericVector s(N);
   NumericVector ln(N);
   NumericVector bn(N);
   NumericVector lnbn(N);
   NumericVector diff(N);
   
   s(0) = lower; 
   ln(0) = round(pow(n, s(0)));
   bn(0) = floor(n / ln(0));
   lnbn(0) = ln(0) * bn(0);
   diff(0) = n - lnbn(0);
   
   for(int i = 1; i < N; i++)
   {
     s(i) = s(i-1) + step;
     ln(i) = round(pow(n, s(i)));
     bn(i) = floor(n / ln(i));
     lnbn(i) = ln(i) * bn(i);
     diff(i) = n - lnbn(i);
   }
   
   DataFrame df = DataFrame::create(Named("n") = n, 
                                    Named("s") = s,
                                    Named("ln") = ln, 
                                    Named("bn") = bn, 
                                    Named("ln*bn") = lnbn, 
                                    Named("diff") = diff);
   return df;
}

//' Optimal parameter s
//' 
//' Calculates the best parameter \eqn{\tilde{s}} for a given approximation s, such that \eqn{n \; % \; \[n^{s}\] = 0}.
//' 
//' @param n Sample size(s), numeric (vector).
//' @param s Desired exponent, \eqn{0.5 \leq s \leq 1}.
//' 
//' @return \code{sOpt} returns a numeric vector of the optimal exponent(s).
//' 
//' @examples 
//' sOpt(50, 0.6)
//' sOpt(100, 0.6)
//' 
//' @export
// [[Rcpp::export]]
NumericVector sOpt(IntegerVector n, double s = 0.6)
{
   int N = n.length();
   NumericVector out(N);
   DataFrame df; 
   NumericVector exp;
   NumericVector diff;
   int j, Ndiff;
   double minDiff;
   
   for(int i = 0; i < N; i++)
   {
     df = sSizes(n(i), 0.5, 1, 0.001);
     exp = df["s"];
     diff = df["diff"];
     Ndiff = diff.length();
     
     j = 0;
     out(i) = exp(0);
     minDiff = R_PosInf;
     
     while(j < Ndiff && (diff(j) != 0 || fabs(s - exp(j)) < minDiff))
     {
       if(diff(j) == 0) 
       {
         minDiff = fabs(s - exp(j));
         out(i) = exp(j);
       }
       
       j++;
     }
   }
   
   return out;
}
