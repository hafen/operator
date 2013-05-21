

blendprop <- function(y, degree, span, span2=nextodd((span - 1)/2), at=1) {
   
   require(operator)

	nextodd <- function(x) {
		x <- round(x)
      x2 <- ifelse(x%%2==0, x+1, x)
		as.integer(x2)
	}

   n <- length(y)
   span3 <- (span - 1)/2
   if(at > span3) at <- span3
      
   ll2 <- loessOp(n=span, span=span, degree=degree)$O[span - at + 1,]
   ll0 <- loessOp(n=span, span=ifelse(degree==1, span, span2), degree=0)$O[span - at + 1,]
   llS <- loessOp(n=span, span=span, degree=degree)$O[(span-1)/2 + 1,]
   
   ind <- (span + 1):(n - span)
   a <- rep(0, length(ind))
   b <- rep(0, length(ind))
   c <- rep(0, length(ind))
   
   for(i in seq_along(ind)) {
      a[i] <- sum(y[(ind[i] - span + at):(ind[i] + at - 1)] * ll0)
      b[i] <- sum(y[(ind[i] - span + at):(ind[i] + at - 1)] * ll2)
      c[i] <- sum(y[(ind[i] - span3):(ind[i] + span3)] * llS)
   }
      
   d <- a - b
   e <- b - c

   p <- -sum(d*e)/sum(d^2)
   
   if(p > 1) p <- 1
   if(p < 0) p <- 0

   p
}