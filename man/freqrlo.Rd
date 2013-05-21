\name{freqrlo}
\alias{freqrlo}
\alias{freqrlo.op}
\alias{freqrlo.loess}
\alias{freqrlo.default}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Frequency response calculation for loess (or general) smoothing matrix
}
\description{
Returns the frequency response function for given loess smoothing parameters or for a specified operator matrix.   Can return values for asymmetric fits.
}
\usage{
\method{freqrlo}{default}(span, degree, M = 400, at = "symmetric")
\method{freqrlo}{op}(op, M = 400, at = 1)
\method{freqrlo}{loess}(x, M = 400, at = "symmetric")
}
\arguments{
   \item{span}{span of loess smooth (must be odd integer).}
   \item{degree}{degree of loess smooth (must be 0, 1, or 2).}
   \item{M}{density of the grid of fourier frequencies, which range from 0 to 0.5.}
   \item{at}{at which design point the frequency response is to be calculated.  This is interpreted as number of points from the end, with 1 being the endpoint.  If an operator matrix is supplied, it is the row of the matrix to use.  For symmetric neighborhoods, the value \code{"symmetric"} can be specified.}
   \item{op}{object of class \code{"op"}.  The operator matrix in this object is used to calculate the frequency response.}
   \item{x}{object of class \code{"loess"}.  Parameters from this objects are passed on to 
   \code{freqrlo.default}.}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
   A data frame of class \code{"freqr"} consisting of frequencies \code{fk} and frequency response values \code{frfv}.
}
\references{
W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

R. Shumway and D. Stoffer. \emph{Time series analysis and its applications.} Springer Verlag, 2000.
}
\author{Ryan Hafen}
% \note{
% %%  ~~further notes~~
% }

\seealso{
\code{\link{freqr}}, \code{\link{plot.freqr}}, \code{\link{loess}}, \code{\link{loessOp}}
}
\examples{
# symmetric loess smooth
a <- freqrlo(11, 2)
plot(a, strip=FALSE, axis=axis.period)

# endpoint loess smooth
b <- freqrlo(11, 2, at=1)
plot(b, strip=FALSE, axis=axis.period)

# accessing from loess
ll <- loess(rnorm(100) ~ c(1:100), degree=1, span=15/100)
ll.fr <- freqrlo(ll)
plot(ll.fr)

# STL overall operator
co2.stl <- stl2(co2, t=as.vector(time(co2)), n.p=12, 
   s.window=35, s.degree=1, t.degree=1, 
   s.blend=0.5, t.blend=0.5)

co2op <- stlOp(co2.stl)
co2fr1 <- freqrlo(co2op$fit, at=1, M=800)
plot(co2fr1, aspect=0.25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
