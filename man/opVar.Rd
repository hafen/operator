\name{opVar}
\alias{opVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute sum of squared values in each row of operator matrix}
\description{
   Compute sum of squared values in each row of operator matrix
}
\usage{
opVar(O)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{O}{a matrix}.
}
\value{
   vector of length \code{nrow(O)}.
}
\author{Ryan Hafen}
\seealso{\code{\link{stlOp}}, \code{\link{loessOp}}, \code{\link{plotVar}}}
