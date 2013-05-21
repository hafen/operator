\name{plot.freqr}
\alias{plot.freqr}
\title{
   Plot power transfer function of objects of class \code{"freqr"}
}
\description{
   Plots the power transfer function of objects of class \code{"freqr"}.
}
\usage{
\method{plot}{freqr}(x, critfreq = NA, 
   xlab = "Frequency", ylab = "Transfer Function", 
   layout = c(1, length(unique(x$which))), 
   type = c("g", "l"), as.table = TRUE, 
   panel = freqrPanel, n.p = NULL, ...)
}
\arguments{
  \item{x}{objects of class \code{"freqr"}.}
  \item{critfreq}{if specified, calculates the critical frequencies for each component plotted that yield a transfer function of \code{critfreq}.  Then a vertical line is plotted where this occurs.}
  \item{n.p}{if specifying \code{critfreq}, must also specify \code{n.p} (see \code{\link{stlOp}}).}
  \item{xlab, ylab, layout, type, as.table, panel, \dots}{arguments passed on to \code{xyplot()}.}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
   An object of class \code{"trellis"}.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{Ryan Hafen}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{freqr}}, \code{\link{freqrlo}}
}
% \examples{
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
