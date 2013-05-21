\name{anova.op}
\alias{anova.op}
\title{
   ANOVA for objects of class \code{"op"}
}
\description{
   ANOVA table calculation for objects of class \code{"op"}
}
\usage{
\method{anova}{op}(mod1, mod2, y)
}
\arguments{
   \item{mod1}{null model operator matrix of class \code{"op"}.}
   \item{mod2}{alternative model operator matrix of class \code{"op"}.}
   \item{y}{data that the models are applied to}
}
\details{
   Calculates ANOVA as described in the reference below.
}
\value{
   An object of class \code{"anova"}.
}
\references{
W. S. Cleveland and S. J. Devlin. Locally weighted regression: An approach to regression analysis by local fitting. \emph{Journal of the American Statistical Association}, 83(403):596--610, 1988.
}
\author{Ryan Hafen}
\note{
   This function treats the model with the larger residual sum of squares as the null model and carries out the test accordingly.  With loess anova, the calculations typically are approximate since exact calculations would involve matrix products, but in this implementation, the calculations are exact since \code{operator} already operates at this level of time complexity.
}
\seealso{
\code{\link{loessOp}}, \code{\link{loessOp}}.
}
\examples{
# a simulated example with periodic component
rf <- function(x) { 
  n <- length(x) 
  sin(x * 2 * pi/200) + rnorm(n, sd = 0.5) + rep(c(0.5, 0.25, 
  0, -0.25, -0.5, -0.25, 0), ceiling(n/7))[1:n] 
}
n <- 200
x <- c(1:n) 
set.seed(8765) 
ysim <- rf(x)

# null model: no seasonality
lop <- loessOp(n, span=105, degree=2)

# alternative model: seasonality
sop <- stlOp(n, n.p=7, t.window=105, t.degree=2, 
   s.window="periodic") 

anova(sop$fit, lop, ysim)
# result:
# [1] Null model: loess with degree=2, span=105
# [2]  Alt model: stl with n.p=7
#                   seasonal smoothing: periodic
#                   trend smoothing: deg=2, span=105
#                   l.degree=2, l.window=7, inner=2
# 
# Analysis of Variance: numerator df 6, denominator df 187.09
# 
#        ENP    RSS F-value    Pr(>F)    
# [1,]  6.05 68.884                      
# [2,] 12.05 52.751  9.5281 3.903e-09 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
