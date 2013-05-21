loessOp <- function(n, span, degree=1, blend=0, at=1:n, stats=TRUE) UseMethod("loessOp")

loessOp.default <- function(n, span, degree=1, blend=0, at=1:n, stats=TRUE) {
   at <- sort(at)

   if(!all(c(1:n) %in% at)) # can only do stats if at has all design points
      stats <- FALSE

   if(span %% 2 != 1)
      stop("span must be odd!")

   nextodd <- function(x) {
       x <- round(x)
       x <- ifelse(x%%2==0, x+1, x)
       as.integer(x)
   }
   
   if(! degree %in% c(0, 1, 2))
      stop("degree must be 0, 1, or 2!")
   
   m <- length(at)
   out <- .C("loess_op", as.integer(n), as.integer(degree), as.integer(span), as.integer(at), as.integer(m), res=double(m*n), PACKAGE="operator")

   O <- matrix(data=out$res, nrow=m, byrow=TRUE)
   # attr(O, "span") <- span
   # attr(O, "deg") <- degree
   
   var <- NULL
   if(all(diff(at) == 1)) {
      var <- opVar(O)
   }
   
   ret <- list(
         O = O,
         at = at,
         var=var,
         pars=list(
            deg=degree,
            span=span
      )
   )
   class(ret) <- "op" 

   # blending
	if(blend > 0 && blend <= 1 && degree >= 1) {
      if(degree == 2)
         sp0 <- nextodd((span+1)/2)
      if(degree == 1)
         sp0 <- span

      O2 <- loessOp(n, span=sp0, degree=0, at=at, stats=FALSE)
      ret <- opBlend(ret, O2, blend)
   }

   if(stats) {
      stats <- getstats(ret)
      ret <- c(ret, list(stats=stats))
   }
   
   ret$call <- match.call()
   class(ret) <- "op" 
   return(ret)
}


loessOp.loess <- function(x, at=1:x$n, stats=TRUE, blend=0) {
   sp <- round(x$par$span * x$n)
   if (sp %% 2 == 0) 
       sp <- sp + 1
   sp <- as.integer(sp)
   
   loessOp.default(x$n, span=sp, degree=x$par$deg, blend=blend, at=at, stats=stats)
}

