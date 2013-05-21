\name{cp}
\alias{cp}
\alias{cp.op}
\alias{cp.stlop}
\alias{cp.loess}
\title{
Calculate Mallows Cp Statistic for Linear Operator
}
\description{
Calculates Mallows Cp statistic for linear operators such as loess and the STL time series modeling method.
}
\usage{
\method{cp}{op}(object, y, sigmasq = 1)
\method{cp}{stlop}(object, y, sigmasq = 1)
\method{cp}{loess}(object, y, sigmasq = 1)
}
\arguments{
  \item{object}{objects of class \code{"op"} or \code{"stlop"}, or \code{"loess"} which provide information from which the Cp statistic is calculated.  If an ARMA model was included in the \code{"stlop"} object, this will be used.}
  \item{y}{the data that the operator is to be applied to.}
  \item{sigmasq}{an unbiased estimate of the residual variance.}
}
\details{
The Cp statistic is calculated as
\deqn{\hat M = RSS/\hat \sigma^2 - tr(\Lambda) + \nu}{M = RSS/sigmahat^2 - tr(Lambda) + nu}
where \eqn{RSS}{RSS} is the residual sum of squares, \eqn{\hat \sigma^2}{sigmahat^2} is the estimated residual variance of an unbiased model, \eqn{\Lambda}{Lambda} is the regulator matrix \eqn{(I - L)'(I-L)}, and \eqn{\nu}{nu} is \eqn{tr(L)}.

If an unbiased estimate of the residual variance is not known, the function can be called with default \code{sigmasq = 1} and then Cp can be reconstructed using the output of the function.
}
\value{
   A data frame consisting of the smoother degrees of freedom \code{df}, which is \eqn{\nu}{nu} in the above equation, the Cp statistic \code{cp}, the estimated residual standard deviation \eqn{\hat \sigma}{sigmahat} as \code{sigmahat}, the trace of the regulator matrix, \eqn{tr(\Lambda)}{tr(Lambda)} as \code{delta1}, and the residual sum of squares.
}
\references{
W. S. Cleveland and S. J. Devlin. Locally weighted regression: An approach to regression analysis by local fitting. \emph{Journal of the American Statistical Association}, 83(403):596--610, 1988.
}
\author{Ryan Hafen}
% \note{
% %%  ~~further notes~~
% }
% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
n <- 500
x <- c(1:n)
set.seed(2456)
y <- sin(x/100) + rnorm(n, sd=0.5)

cpres <- NULL
for(sp in seq(251, 351, by=5)) {
   ll <- loess(y ~ x, degree=2, span=sp/n)
   cpres <- rbind(cpres,
      cp(ll, sigmasq=0.5^2)
   )
}

xyplot(cp ~ df, data=cpres)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
