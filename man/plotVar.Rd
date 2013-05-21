\name{plotVar}
\alias{plotVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of the variance factor of the loess or stl fits at each design 
point}
\description{
Trellis plot of the values of the variance factor of the loess or stl fits 
at each design point.
}
\usage{
plotVar(op, ylab = "variance", xlab = "x", 
   panel = var.panel, auto.key = 
   list(lines = TRUE, points = FALSE, columns = 3), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{op}{object of class "op", "stlop", or "opblend"}
  \item{ylab, xlab, panel, auto.key, ...}{trellis plotting parameters}
}
\details{
   If \code{op} is of class "op", then the values of \code{op\$var} are 
   plotted.  If \code{op} is of class "opblend", the variance of the two 
   fits and the blended variance are overplotted.  If \code{op} is of 
   class "stlop", the variance for the fitted values, trend, and seasonal
   components are overplotted.
}
\value{
   an object of class "trellis"
}
\references{
R. B. Cleveland, W. S. Cleveland, J.E. McRae, and I. Terpenning (1990) 
STL: A Seasonal-Trend Decomposition Procedure Based on Loess. 
\emph{Journal of Official Statistics}, 6, 3-73.

W.S. Cleveland, E. Grosse and W.M. Shyu (1992) Local regression
models. Chapter 8 of \emph{Statistical Models in S} eds J.M. Chambers
and T.J. Hastie, Wadsworth & Brooks/Cole.
}
\author{Ryan Hafen}
\seealso{\code{\link{stl2}}}
\examples{

n <- 200
L1 <- loessOp(n, span=91, degree=2, at=c(-9:(n+10)))
L2 <- loessOp(n, span=45, degree=0, at=c(-9:(n+10)))
Lb <- opBlend(L1, L2, n.b=50, blend.p=0.5)

p <- plotVar(L1)
p
plotVar(L2, ylim=p$y.limits)
plotVar(Lb)

sop <- stlOp(50, n.p=7, s.window="periodic", s.degree=1, n.ahead=7)
plotVar(sop)
plotVar(sop$S)

}
