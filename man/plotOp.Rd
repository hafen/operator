\name{plotOp}
\alias{plotOp}
\title{A "heat map" plot of the operator matrix}
\description{
Plots the row and column values of the operator matrix as a ``heat map''.
}
\usage{
plotOp(op, nbreaks = 30, panel = op.panel, aspect = "xy", 
   xlab = "Column", ylab = "Row", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{op}{object of class \code{"op"}, \code{"opblend"}, or \code{"stlop"}.}
  \item{nbreaks}{number of breaks in the color graduation.}
  \item{panel, aspect, xlab, ylab, ...}{override the default arguments to the lattice call.}
}
\details{
   The color white correspond to values around 0.
}
\value{
   an object of class \code{"trellis"}.
}
\references{
R. B. Cleveland, W. S. Cleveland, J.E. McRae, and I. Terpenning (1990) STL: A Seasonal-Trend Decomposition Procedure Based on Loess. \emph{Journal of Official Statistics}, 6, 3-73.

W. S. Cleveland, E. Grosse and W. M. Shyu (1992) Local regression models. Chapter 8 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
}
\author{Ryan Hafen}
\seealso{\code{\link{stl2}}}
\examples{
sop <- stlOp(75, n.p=7, s.window="periodic")
plotOp(sop$fit)
}
