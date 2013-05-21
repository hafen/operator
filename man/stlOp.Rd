\name{stlOp}
\alias{stlOp}
\alias{stlOp.stl}
\alias{stlOp.stl2}
\alias{stlOp.default}
\title{Calculate operator matrices for the STL decomposition}
\description{
Calculate operator matrices for the STL decomposition for a univariate equally-spaced design time series.
}
\usage{
\method{stlOp}{default}(n, n.p, s.window, s.degree = 1, 
   t.window = NULL, t.degree = 1, fc.window = NULL, fc.degree = NULL,
   l.window = nextodd(n.p), l.degree = t.degree, critfreq
   = 0.05, inner = 2, n.ahead = 0, s.blend = 0, t.blend =
   0, l.blend = t.blend, fc.blend = NULL, fc.name = NULL,
   arma = NULL, stats = TRUE)
\method{stlOp}{stl}(x, n.ahead=0, s.blend=0, t.blend=0, 
   l.blend=t.blend, critfreq=0.05, arma=NULL, stats=TRUE)
\method{stlOp}{stl2}(x, n.ahead=0, arma=NULL, stats=TRUE)
}
\arguments{
   \item{n}{number of observations in the time series.}
   \item{x}{object of class \code{"stl"} or \code{"stl2"}.}
   \item{n.p}{the periodicity of the seasonal component.}
   
 	\item{s.window}{either the character string \code{"periodic"} or the span (in lags) of the loess window for seasonal extraction, which should be odd.  This has no default.}
 	\item{s.degree}{degree of locally-fitted polynomial in seasonal extraction.  Should be 0, 1, or 2.}
 	\item{t.window}{the span (in lags) of the loess window for trend extraction, which should be odd.  If \code{NULL}, the default, \code{nextodd(ceiling((1.5*period) / (1-(1.5/s.window))))}, is taken.}
 	\item{t.degree}{degree of locally-fitted polynomial in trend extraction.  Should be 0, 1, or 2.}
 	\item{l.window}{the span (in lags) of the loess window of the low-pass filter used for each subseries.  Defaults to the smallest odd integer greater than or equal to \code{frequency(x)} which is recommended since it prevents competition between the trend and seasonal components.  If not an odd integer its given value is increased to the next odd one.}
 	\item{l.degree}{degree of locally-fitted polynomial for the subseries low-pass filter.  Should be 0, 1, or 2.}
   \item{critfreq}{the critical frequency to use for automatic calculation of smoothing windows for the trend and high-pass filter.}
   \item{n.ahead}{number of time units into the future to compute the operator matrices for.  This will be used for prediction.}
   \item{fc.window}{vector of lengths of windows for loess smoothings for other trend frequency components after the original STL decomposition has been obtained.  The smoothing is applied to the data with the STL seasonal component removed.  A frequency component is computed by a loess fit with the window length equal to the first element of fc.window, the component is removed, another component is computed with the window length equal to the second element of fc.window, and so forth. In most cases, the values of the argument should be decreasing, that is, the frequency bands of the fitted components should increase. The robustness weights from original STL are used as weights in the loess fitting if specified.}
   \item{fc.degree}{vector of degrees of locally-fitted polynomial in the loess smoothings for the frequency components specified in fc.window. Values of 0, 1 and 2 are allowed. If the length of fc.degree is less than that of fc.window, the former is expanded to the length of the latter using rep; thus, giving the value 1 specifies a degree of 1 for all components.}
   \item{fc.name}{vector of names of the post-trend smoothing operations specified by \code{fc.window} and \code{fc.degree} (optional).}
   \item{inner}{integer; the number of \sQuote{inner} (backfitting) iterations; usually very few (2) iterations suffice.}
   \item{s.blend, t.blend, l.blend, fc.blend}{vectors of proportion of blending to degree 0 polynomials at the endpoints of the series.}
 	\item{arma}{object of class \code{"Arima"}, which specifies the arma model to use on the remainder component.}
   \item{stats}{whether or not to calculate auxiliary statistics for the overall STL fit operator matrix.}
}
\details{
   The STL operator matrix is the matrix that yields the fitted values (seasonal + trend) through a linear operation with the observed time series.  It is obtained through a series of linear filters.  If post-trend smoothing is specified through \code{fc.window, \dots}, then the overall operator matrix will be the seasonal operator plus the operator for each post-trend smoothing.  If ARMA modeling is specified, this will also be factored in to the calculation of the overall operator matrix.
}
\value{
   a list of class \code{"stlop"}.
   \item{fit}{the overall operator matrix, a list of class \code{"op"}, which is calculated as described in the details section.  If \code{stats=TRUE}, then auxiliary statistics will also be computed and will be available as an element of this list.}
   \item{seas}{the operator matrix for just the seasonal component, an object of class \code{"op"}.}
  \item{trend}{the operator matrix for just the trend component, a list of class \code{"op"}.}
  \item{fc}{a list of operator matrices corresponding to post-trend frequency components, if specified.}
  \item{at}{the values at which the operator matrices were calculated, which will either be \code{1, \dots, n} or \code{1, \dots, n + n.ahead}.}
  \item{pars}{parameters used in the fitting.}
}
\references{
R. B. Cleveland, W. S. Cleveland, J.E. McRae, and I. Terpenning (1990) STL: A Seasonal-Trend Decomposition Procedure Based on Loess. \emph{Journal of Official Statistics}, 6, 3-73.

W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{Ryan Hafen}
\note{
   This function does a lot of \code{n}x\code{n} matrix multiplication.  Be aware of this when choosing \code{n}.  While the loess operator matrices are calculated in C, the matrix multiplication happens in R.  This is because of the speed of BLAS, which is especially good if R has been compiled with threaded BLAS.
}
\seealso{\code{\link{stl}}, \code{\link{plotVar}}, 
\code{\link{loessOp}}, \code{\link{predict.stlop}}}
\examples{
sop <- stlOp(50, n.p=7, s.window="periodic")
plotVar(sop)
plotVar(sop$seas)

# stl blending the trend
t.blend <- list(deg=0, span=11, blend.p=0.5, n.b=20)
sop2 <- stlOp(50, n.p=7, s.window="periodic", t.blend=t.blend)
plotVar(sop2)
plotOp(sop$fit)

# predicting ahead
sop <- stlOp(100, n.p=7, s.window=51, s.degree=1, n.ahead=7)
plotVar(sop)

# now stl + further loess smoothing, predicting ahead 7 days
# first get day-of-week component, then smooth with span 1001
# followed by smooth of span 91
n <- 200
rop <- stlOp(n, 7, s.window="periodic", t.window=39, n.ahead=7, 
   fc.window=c(1001, 91), fc.degree=c(1, 2), fc.blend=c(0, 0.5))
plotVar(rop)

# do same thing the hard way (without specifying fc.degree, etc.)
# just to illustrate handling elements of stlop objects
tmp <- stlOp(n, 7, s.window="periodic", t.window=39, n.ahead=7)
dow <- tmp$seas$O
fc1 <- loessOp(n, span=1001, degree=1, at=1:207)$O
fc2 <- loessOp(n, span=91, degree=2, at=1:207, blend=0.5)$O
rop2 <- dow + fc1 - fc1 %*% dow[1:n,] + fc2 - fc2 %*% dow[1:n,] - 
   fc2 %*% fc1[1:n,] + fc2 %*% fc1[1:n,] %*% dow[1:n,]

}
