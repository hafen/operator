\name{loessOp}
\alias{loessOp}
\alias{loessOp.default}
\alias{loessOp.loess}
\title{Loess operator matrix calculation}
\description{
Calculates the hat matrix, L, used to obtain fitted values for given loess smoothing parameters, for an equally-spaced design \code{1, \dots, n}.  Also allows extrapolation beyond the design points.
}
\usage{
loessOp(n, span, degree = 1, blend = 0, at = 1:n, 
   stats = TRUE)
\method{loessOp}{default}(n, span, degree=1, blend=0, at=1:n, stats=TRUE)
\method{loessOp}{loess}(x, at=1:x$n, stats=TRUE, blend=0)
}
\arguments{
  \item{n}{number of equally-spaced design points, assigned values \code{1, \dots, n}.}
  \item{x}{object of class \code{"loess"}.}
  \item{span}{number of design points to be used in the local neighborhood -- must be odd.}
  \item{degree}{degree of local polynomial (currently can be either 0, 1, or 2).}
  \item{blend}{the amount of blending to degree 0 smoothing at the endpoints.}
  \item{at}{which rows of the operator matrix to calculate.}
  \item{stats}{whether or not to calculate auxiliary statistics.}
}
\details{
   If all that is desired is the loess estimate, it is more efficient and flexible to use the built-in \code{loess()}.  The main purpose of this function is its use in stlOp or subsequent smoothings after using \code{stlOp()}.
}
\value{
   A list of class \code{"op"}.
  \item{O}{the operator matrix of dimension \code{length(at) x n}.}
  \item{at}{at as specified in \code{loessOp()}.}
  \item{span}{span as specified in \code{loessOp()}.}
  \item{deg}{deg as specified in \code{loessOp()}.}
  \item{var}{the squared l2 norm of each row in O -- used in variance calculations.}
  \item{stats}{only if \code{stats=TRUE}.  This is a list with the regulator matrix \code{lambda}=(I-O)'(I-O), the equivalent number of parameters \code{enp} (tr O'O), \code{delta1} (tr lambda), \code{delta2} (tr lambda^2), and \code{trace.hat} (tr O).}
}
\references{
W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{Ryan Hafen}
\note{
This requires o(n^2) storage and computation.  Rows of the operator matrix corresponding to interior design points will be identical, so keep this in mind.
}
\seealso{\code{\link{stlOp}}, \code{\link{plotVar}}, 
\code{\link{plotOp}}}
\examples{
n <- 200
L <- loessOp(n, span=91, degree=2, at=c(-9:(n+10)))
plotVar(L)
plotOp(L)

## get fitted values for some data
x <- c(1:n)
# generate some data
y <- sin(x*2*pi/n) + rnorm(x, sd=0.5) 
# get the fitted values
yhat <- predict(L, y)
# another way: yhat <- L$O %*% y

plot(x, y, xlim=range(L$at))
lines(L$at, yhat, col="red")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
