\name{freqr}
\alias{freqr}
\alias{freqr.stl}
\alias{freqr.stl2}
\alias{freqr.stlop}
\alias{freqr.default}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Frequency response calculation for STL decomposition
}
\description{
Calculates the frequency response function for the seasonal, trend, and high-pass filters of an STL decomposition at a specified grid of fourier frequencies for the symmetric smoothing case.  This can be used to visualize the power transfer function or phase shift function to evaluate the spectral properties of an STL model.
}
\usage{
\method{freqr}{default}(n.p, s.window, s.degree = 0, 
   t.window = NULL, t.degree = 1, l.window = nextodd(n.p), 
   l.degree = t.degree, periodic = FALSE, M = 400)
\method{freqr}{stl}(x, ...)
\method{freqr}{stl2}(x, ...)
\method{freqr}{stlop}(x, ...)
}
\arguments{
  \item{n.p, 
  s.window, 
  s.degree, 
  t.window, 
  t.degree, 
  l.window, 
  l.degree, 
  periodic}{STL parameters (see \code{\link{stlOp}}).}
  \item{x}{object of class \code{"stl"}, \code{"stl2"}, or \code{"stlop"}.  Parameters from these objects are passed on to \code{freqr.default}.}
  \item{M}{density of the grid of fourier frequencies, which range from 0 to 0.5.}
  \item{...}{}
}
% \details{
% }
\value{
   A data frame of class \code{"freqr"} consisting of frequencies \code{fk} and frequency response values \code{frfv}.
}
\references{
R. B. Cleveland, W. S. Cleveland, J. E. McRae, and I. Terpenning (1990) STL: A Seasonal-Trend Decomposition Procedure Based on Loess. \emph{Journal of Official Statistics}, 6, 3-73.

R. Shumway and D. Stoffer. \emph{Time series analysis and its applications.} Springer Verlag, 2000.
}
\author{
Ryan Hafen
}
\note{
To investigate the spectral properties of an STL model at an asymmetric smoothing location, one can use the \code{\link{freqrlo}} function with a \code{"stlop"} object to obtain the frequency response function at any smoothing point.
}
\seealso{
\code{\link{stlOp}}, \code{\link{freqrlo}}, \code{\link{plot.freqr}}
}
\examples{
co2.stl <- stl2(co2, t=as.vector(time(co2)), n.p=12, 
   s.window=35, s.degree=1, t.degree=1, 
   s.blend=0.5, t.blend=0.5)

co2op <- stlOp(co2.stl)

# symmetric filter
co2fr <- freqr(co2op, M=1000)
# seasonal, trend, and overall filter at endpoint
co2fr1 <- freqrlo(co2op$fit, at=1, M=800)
co2fr1.s <- freqrlo(co2op$seas, at=1)
co2fr1.t <- freqrlo(co2op$trend, at=1)
# overall filter at 10 points from end and middle
co2fr10 <- freqrlo(co2op$fit, at=10, M=800)
co2frm <- freqrlo(co2op$fit, at=234, M=800)

plot(co2fr, aspect=0.25)
plot(co2fr1, aspect=0.25)
plot(co2fr1.s, aspect=0.25)
plot(co2fr1.t, aspect=0.25)
plot(co2frm, aspect=0.25)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
