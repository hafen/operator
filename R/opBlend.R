opBlend <- function(from, to, blend.p, n.b=as.integer(from$pars$span/2)) {

   O1 <- from$O
   O2 <- to$O

   if(!all(from$at==to$at))
      stop("both operator matrices must be evaluated at the same design points")

   if(any(diff(from$at) > 1))
      stop("operator matrix must be evaluated at each design point")

   n <- ncol(O1)
   n2 <- nrow(O1)

   beg <- min(from$at)
   end <- max(from$at)

   n.after <- end - n
   n.before <- 1 - beg

   w1 <- matrix(nrow=n.b, ncol=n, data=rep(seq(1, 1 - blend.p, length=n.b), n))
   # if there is predicting ahead, keep the same blending proportion
   w1 <- rbind(w1, matrix(nrow=n.after, ncol=n, data=1 - blend.p))

   w11 <- matrix(nrow=n.b, ncol=n, data=rep(rev(seq(1, 1 - blend.p, length=n.b)), n))
   w11 <- rbind(matrix(nrow=n.before, ncol=n, data=1 - blend.p), w11)

   res <- O1
   right.index <- max(1, (n2 - n.b - n.after + 1)):n2
   left.index <- 1:min(n2, (n.b + n.before))

   if(length(right.index) > n2/2) {
      m <- median(from$at)
      right.index <- right.index[from$at[right.index] <= m]
      left.index <- left.index[from$at[left.index] > m]
      w11 <- w11[1:length(left.index),]
      w1 <- w1[(nrow(w1) - length(right.index) + 1):nrow(w1),]
   }

   w2 <- 1 - w1
   w22 <- 1 - w11

   res[right.index,] <- O1[right.index,] * w1 + O2[right.index,] * w2
   res[left.index,] <- O1[left.index,] * w11 + O2[left.index,] * w22

   ret <- list(
         O = res,
         at = from$at,
         blend.p = blend.p,
         n.b = n.b,
         var = opVar(res),
         from.var = from$var,
         to.var = to$var,
         from.span = from$pars$span,
         from.deg= from$pars$deg,
         to.span = to$pars$span,
         to.deg = to$pars$deg
      )
   class(ret) <- c("opblend", "op")
   return(ret)
}

