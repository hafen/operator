\name{predict.op}
\alias{predict.op}
\alias{predict.stlop}
\title{Obtain predicted values for a given operator matrix and data vector}
\description{
   Returns fitted values, optionally with standard errors.
}
\usage{
\method{predict}{op}(op, y, se = FALSE, newdata = op$at, interval = c("none", 
   "confidence", "prediction"), level = 0.95)
\method{predict}{stlop}(stlop, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{op}{object of class \code{"op"} or \code{"opblend"}.}
  \item{stlop}{object of class \code{"stlop"}.}
  \item{y}{data vector with same length as the number of columns of the operator matrix in \code{op} or the fit operator matrix in \code{stlop}.}
  \item{se}{should standard errors be computed?}
  \item{newdata}{at which values should fits be computed?  Valid values are 1 to \code{n.ahead}.}
  \item{interval}{should confidence or prediction intervals be computed?}
  \item{level}{level for the confidence or prediction intervals.}
  \item{...}{}
}
\details{
	If \code{se=TRUE} or \code{interval} is \code{"confidence"} or \code{"prediction"}, and the argument \code{op} does not have a \code{stats} element, the auxiliary statistics will be computed.  Also, if 
	\code{interval="prediction"} and \code{newdata} was not specified, it will by default be all values beyond the length of the original series, up to \code{n.ahead}.
}
\value{If \code{se = FALSE} and \code{interval="none"}, a vector giving the prediction for each point in the design space. If \code{se = TRUE} or \code{interval} is \code{"confidence"} or \code{"prediction"}, a list containing a data frame \code{data} with components
	\item{x}{the time values.}
	\item{y}{the observed series values.}
   % another data frame \code{fits} with components
	\item{at}{the time values at which the fit was computed.}
   \item{fit}{the fitted/predicted values.}
   \item{se.fit}{an estimated standard error for each predicted value.}
	\item{lower}{lower confidence/prediction limit, if requested.}
	\item{upper}{upper confidence/prediction limit, if requested and additional elements.}
   \item{residual.scale}{the estimated scale of the residuals used in computing the standard errors.}
   \item{df}{an estimate of the effective degrees of freedom used in estimating the residual scale, intended for use with t-based confidence intervals.}
}
\references{
R. B. Cleveland, W. S. Cleveland, J. E. McRae, and I. Terpenning (1990) STL: A Seasonal-Trend Decomposition Procedure Based on Loess. \emph{Journal of Official Statistics}, 6, 3-73.

W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{Ryan Hafen}
\note{
   All of this can be done for loess fitting easily using \code{loess()} and \code{predict.loess()}.  More important is \code{predict.stlop()} for obtaining predicted values and standard errors for the STL decomposition.
}
\seealso{\code{\link{stlOp}}, \code{\link{loessOp}}}
\examples{
co2op <- stlOp(468, n.p = 12, l.window = 13, t.window = 19, 
   s.window = 35, s.degree = 1, n.ahead = 36)

co2pi <- predict(co2op, co2, interval="prediction") 

plot(co2pi, start=400, xlab="Time", ylab="CO2Concentration(ppm)") 

# a simulated example
rf <- function(x) { 
  n <- length(x) 
  sin(x * 2 * pi/200) + rnorm(n, sd = 0.5) + rep(c(0.5, 0.25, 
  0, -0.25, -0.5, -0.25, 0), ceiling(n/7))[1:n] 
}
n <- 200 
x <- 1:n 
set.seed(8765) 
ysim <- rf(x)

sop <- stlOp(200, n.p=7, t.window=105, t.degree=2, s.window="periodic") 
ysimCI <- predict(sop, ysim, interval="confidence", level=0.99) 
plot(ysimCI, fcol="gray", CIalpha=0.75) 

}