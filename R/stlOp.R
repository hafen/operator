stlOp <- function(n, n.p, s.window, s.degree=1, t.window=NULL, t.degree=1, fc.window=NULL, fc.degree=NULL, l.window=nextodd(n.p), l.degree=t.degree, critfreq=0.05, inner=2, n.ahead=0, s.blend=0, t.blend=0, l.blend=t.blend, fc.blend=NULL, fc.name=NULL, arma=NULL, stats=TRUE) UseMethod("stlOp")

stlOp.default <- function(n, n.p, s.window, s.degree=1, t.window=NULL, t.degree=1, fc.window=NULL, fc.degree=NULL, l.window=nextodd(n.p), l.degree=t.degree, critfreq=0.05, inner=2, n.ahead=0, s.blend=0, t.blend=0, l.blend=t.blend, fc.blend=NULL, fc.name=NULL, arma=NULL, stats=TRUE) {

   n.p <- as.integer(n.p)

   if(n.ahead < 0)
      stop("n.ahead must be a non-negative integer")
   
   if(!is.null(arma)) {
      if(!any(class(arma) == "Arima"))
         stop("argument 'arma' must be of class 'Arima' (fitted with arima() in stats package)")
   }
   
   periodic <- FALSE

   nextodd <- function(x) {
       x <- round(x)
       x <- ifelse(x%%2==0, x+1, x)
       as.integer(x)
   }
   
   # deg.check <- function(deg) {
   #     degname <- deparse(substitute(deg))
   #     deg <- as.integer(deg)
   #     if (deg < 0 || deg > 2)
   #         stop(degname, " must be 0, 1, or 2")
   #     deg
   # }
   
   wincheck <- function(x) {
      x <- nextodd(x)
	   if(any(x <= 0)) stop("Window lengths must be positive.")
      x
	}
	
	degcheck <- function(x) {
	   if(! all(x==0 | x==1 | x==2)) stop("Smoothing degree must be 0, 1, or 2")
	}
	
   # new approximation for t.degree
   get.t.window <- function(t.dg, s.dg, n.s, n.p, omega) {
      if(t.dg == 0) t.dg <- 1
      if(s.dg == 0) s.dg <- 1
      
      coefs_a <- data.frame(a = c(0.000103350651767650, 3.81086166990428e-06
      ), b = c(-0.000216653946625270, 0.000708495976681902))
      coefs_b <- data.frame(a = c(1.42686036792937, 2.24089552678906
      ), b = c(-3.1503819836694, -3.30435316073732), c = c(5.07481807116087, 
      5.08099438760489))
      coefs_c <- data.frame(a = c(1.66534145060448, 2.33114333880815
      ), b = c(-3.87719398039131, -1.8314816166323), c = c(6.46952900183769, 
      1.85431548427732))
      
      # estimate critical frequency for seasonal
      betac0 <- coefs_a$a[s.dg] + coefs_a$b[s.dg] * omega
      betac1 <- coefs_b$a[s.dg] + coefs_b$b[s.dg] * omega + coefs_b$c[s.dg] * omega^2
      betac2 <- coefs_c$a[s.dg] + coefs_c$b[s.dg] * omega + coefs_c$c[s.dg] * omega^2
      f_c <- (1 - (betac0 + betac1 / n.s + betac2 / n.s^2)) / n.p
      
      # choose 
      betat0 <- coefs_a$a[t.dg] + coefs_a$b[t.dg] * omega
      betat1 <- coefs_b$a[t.dg] + coefs_b$b[t.dg] * omega + coefs_b$c[t.dg] * omega^2
      betat2 <- coefs_c$a[t.dg] + coefs_c$b[t.dg] * omega + coefs_c$c[t.dg] * omega^2

      betat00 <- betat0 - f_c

      n.t <- nextodd((-betat1 - sqrt(betat1^2 - 4*betat00*betat2)) / (2 * betat00))
      
      n.t
   }
   
   # blend.check <- function(blend) {}

	if(is.null(l.window)) {
	   l.window <- nextodd(n.p)
	} else {
	   l.window <- wincheck(l.window)
	}
   
   if (is.character(s.window)) {
       if (is.na(pmatch(s.window, "periodic"))) 
           stop("unknown string value for s.window")
       else {
           periodic <- TRUE
           s.window <- 10 * n + 1
           s.degree <- 0
       }
   }

	degcheck(s.degree); degcheck(t.degree); degcheck(l.degree)
   # s.degree <- deg.check(s.degree)
   # t.degree <- deg.check(t.degree)
   # l.degree <- deg.check(l.degree)

	if (is.null(t.window)) {	   
      # t.window <- nextodd(ceiling(1.5 * n.p/(1 - 1.5/s.window)))
		t.window <- get.t.window(t.degree, s.degree, s.window, n.p, critfreq)
	} else {
	   t.window <- wincheck(t.window)
	}

	n2 <- n + n.ahead
	at <- c(1:n2)
   aa <- ceiling(n/n.p)
   aa2 <- ceiling(n2/n.p)
	

   I <- matrix(nrow=n, ncol=n, data=0)
   I[1:n, 1:n] <- diag(nrow=n, ncol=n)


	cs <- rep(1:n.p, aa)[1:n] # cycle-subseries indices
   C <- matrix(nrow=n.p*(aa2+2), ncol=n, data=0)

   if(periodic) {
		for(i in 1:n.p) {
			csi <- which(cs==i)
			C[((0:(aa2+2-1))*(n.p) + i), csi] <- matrix(data=rep(1/length(csi), length(csi)), nrow=aa2+2, ncol=length(csi), byrow=TRUE)
		}
   } else {
		for(i in 1:n.p) {
			csi <- which(cs==i)
			
			tmp <- loessOp(length(csi), span=s.window, degree=s.degree, at=0:(aa2 + 1), stats=FALSE)

         # blend seasonal smoothing
         if(s.blend > 0 && s.blend <= 1 && s.degree >= 1) {
            if(s.degree == 2)
               s.sp0 <- nextodd((s.window+1)/2)
            if(s.degree == 1)
               s.sp0 <- s.window

            tmp2 <- loessOp(length(csi), span=s.sp0, degree=0, at=0:(aa2 + 1), stats=FALSE)
            tmp <- opBlend(tmp, tmp2, s.blend)
         }

			C[((0:(aa2+2-1))*(n.p) + i), csi] <- tmp$O
		}
   }

	C <- C[1:(n2 + 2*n.p),]

   P <- matrix(nrow=n2, ncol=n2 + 2*n.p, data=0)
   P[, (n.p+1):(n2 + n.p)] <- diag(1, n2, n2)

   # PC <- P %*% C

   T <- loessOp(n, span=t.window, degree=t.degree, at=at, stats=FALSE)
   
   if(t.blend > 0 && t.blend <= 1 && t.degree >= 1) {
      if(t.degree == 2)
         t.sp0 <- nextodd((t.window+1)/2)
      if(t.degree == 1)
         t.sp0 <- t.window

      T2 <- loessOp(n, span=t.sp0, degree=0, at=at, stats=FALSE)
      T <- opBlend(T, T2, t.blend)
   }

   T <- T$O

   # get the values for each row of the moving average matrix
   b1 <- c(0:(n.p-1))*3
   b1[1] <- 1
   b <- c(b1, b1[n.p]+1, rev(b1))/(n.p^2*3)

   l.b <- length(b)

   MA <- matrix(nrow=n2, ncol=n2 + 2*n.p, data=0)
   for(i in 1:n2) {
      MA[i, i:(i+l.b-1)] <- b
   }

   # loess smoothing of the moving averages operator matrix
   L1 <- loessOp(n2, span=l.window, degree=l.degree, at=at, stats=FALSE)

   if(l.blend > 0 && l.blend <= 1 && l.degree >= 1) {
      if(l.degree == 2)
         l.sp0 <- nextodd((l.window+1)/2)
      if(l.degree == 1)
         l.sp0 <- l.window

      L2 <- loessOp(n2, span=l.sp0, degree=0, at=at, stats=FALSE)
      L1 <- opBlend(L1, L2, l.blend)$O
   } else {
      L1 <- L1$O
   }

   # now finally we have L
   L <- L1 %*% MA

   S <- (P - L) %*% C

   # first iteration...
   S.new <- S
   T.new <- T %*% (I - S.new[1:n,])

   # subsequent iterations...
   if(inner > 1) {
      for(i in 2:inner) {
      	S.new <- S %*% (I - T.new[1:n,])
      	T.new <- T %*% (I - S.new[1:n,])
      }      
   }
   
   # now post-trend smoothing, if specified...
   
   fc <- NULL
   fc.number <- 0
   fc.res <- NULL
   if(!is.null(fc.window)) {
      fc.number <- length(fc.window)
      fc <- list()
      length(fc) <- fc.number
      
      fc.window <- wincheck(fc.window)
      if(is.null(fc.degree)) fc.degree <- 1
      if(length(fc.degree) < fc.number)
         fc.degree <- c(fc.degree, rep(fc.degree[length(fc.degree)], fc.number - length(fc.degree)))
      fc.cumulative <- matrix(nrow=n2, ncol=n, data=0) # keep sum of all previous fc smoothings
      fc.cumulative2 <- S.new # keep sum of all fc operators

      degcheck(fc.degree)

      if(is.null(fc.name))
         fc.name <- paste("fc.", fc.window, sep="")

      if(is.null(fc.blend))
         fc.blend <- rep(t.blend, fc.number)
         
      if(length(fc.blend) < fc.number)
         fc.blend <- c(fc.blend, rep(0, fc.number - length(fc.blend)))
         
      fc.res <- list()
      for(ii in 1:fc.number) {      

         FCtmp <- loessOp(n, span=fc.window[ii], degree=fc.degree[ii], at=at, stats=FALSE)

         if(fc.blend[ii] > 0 && fc.blend[ii] <= 1 && fc.degree[ii] >= 1) {
            if(fc.degree[ii] == 2)
               fc.sp0 <- nextodd((fc.window[ii]+1)/2)
            if(fc.degree[ii] == 1)
               fc.sp0 <- fc.window[ii]

            FC0 <- loessOp(n, span=fc.sp0, degree=0, at=at, stats=FALSE)
            FCtmp <- opBlend(FCtmp, FC0, fc.blend[ii])
         }
         
         fc[[ii]] <- FCtmp$O %*% (I - S.new[1:n,] - fc.cumulative[1:n,])
         
         fc.cumulative <- fc.cumulative + FCtmp$O
         fc.cumulative2 <- fc.cumulative2 + fc[[ii]]
      }
   
      if(any(is.null(fc.name)) || length(fc.name) < fc.number)
         fc.name <- paste("trend", c(1:fc.number))

      names(fc) <- fc.name
      O <- fc.cumulative2
   } else {
      O <- S.new + T.new
   }
   
   O <- list(O=O,
             at=at,
             var=opVar(S.new + T.new),
             type="operator",
             call=match.call(),
             pars=data.frame(
                s.degree=s.degree,
                t.degree=t.degree,
                l.degree=l.degree,
                s.window=s.window,
                t.window=t.window,
                l.window=l.window,
                inner=inner,
                n.p=n.p, periodic=periodic)
   )
   class(O) <- "op"
   
   # ARMA to remove autocorrelation
   if(!is.null(arma)) {
      tmp <- armaOp(arma)
      if(n.ahead > 0) {
         tmp <- predict(tmp, n.ahead=n.ahead)
      } else {
         tmp <- tmp$O
      }
      O2 <- tmp + O$O - tmp %*% O$O[1:ncol(tmp),]      
      
      O2 <- list(O=O2,
                at=at,
                var=opVar(O2),
                type="operator with ARMA",
                pars=data.frame(
                   s.degree=s.degree,
                   t.degree=t.degree,
                   l.degree=l.degree,
                   s.window=s.window,
                   t.window=t.window,
                   l.window=l.window,
                   inner=inner,
                   n.p=n.p, 
                   periodic=periodic),
                armamod=arma
      )
      
      class(O2) <- "op"
      
      if(stats) {
         stat <- getstats(O2)
         O2 <- c(O2, list(stats=stat))
      }
      
      class(O2) <- "op"
   }
   
   S <- list(O=S.new,
             at = at,
             var=opVar(S.new),
             type="seasonal"
   )
   class(S) <- "op"
   T <- list(O=T.new,
             at = at,
             var=opVar(T.new),
             type="trend"
   )
   class(T) <- "op"

   if(stats) {
      stat <- getstats(O)
      O <- c(O, list(stats=stat))
      class(O) <- "op" 
   }

   ret <- list(fit=O,
               seas=S,
               trend=T,
               fc=fc,
               at=at,
               pars=list(
                  s.degree=s.degree,
                  t.degree=t.degree,
                  l.degree=l.degree,
                  fc.degree=fc.degree,
                  s.window=s.window,
                  t.window=t.window,
                  l.window=l.window,
                  fc.window=fc.window,
                  s.blend=s.blend,
                  t.blend=t.blend,
                  l.blend=l.blend,
                  fc.blend=fc.blend,
                  critfreq=critfreq,
                  inner=inner,
                  fc.number=fc.number,
                  fc.name=fc.name,
                  n.p=n.p, periodic=periodic
               )
   )
   if(!is.null(arma)) {
      ret <- c(ret, list(armafit=O2))
   }
   
   ret$call <- match.call()
   class(ret) <- "stlop"
   
   return(ret)
}

# methods for getting parameters from stl and stl2
stlOp.stl2 <- function(x, n.ahead=0, arma=NULL, stats=TRUE) {
   if(!(x$pars$outer == 1))
      warning("outer > 1, operator currently does not have robustness implemented")
   
   stlOp.default(n=length(x$data$raw), 
                 n.p=x$pars$n.p, 
                 s.window=x$pars$win$s.window, 
                 s.degree=x$pars$deg$s.degree, 
                 t.window=x$pars$win$t.window, 
                 t.degree=x$pars$deg$t.degree, 
                 fc.degree=x$pars$fc$fc.degree,
                 fc.window=x$pars$fc$fc.window,
                 l.window=x$pars$win$l.window, 
                 l.degree=x$pars$deg$l.degree, 
                 inner=x$pars$inner, 
                 n.ahead=n.ahead, 
                 s.blend=x$pars$blend$s.blend, 
                 t.blend=x$pars$blend$t.blend, 
                 l.blend=x$pars$blend$l.blend, 
                 fc.blend=x$pars$fc$fc.blend,
                 fc.name=x$pars$fc$fc.name,
                 critfreq=x$pars$critfreq,
                 arma=arma,
                 stats=stats)
}

stlOp.stl <- function(x, n.ahead=0, s.blend=0, t.blend=0, l.blend=t.blend, critfreq=0.05, arma=NULL, stats=TRUE) {
   if(!(x$outer == 0))
      warning("outer > 0, operator currently does not have robustness implemented")
   
   stlOp.default(n=length(x$time.series[,1]), 
                 n.p=frequency(x$time.series[,1]), 
                 s.window=as.integer(x$win[1]), 
                 s.degree=as.integer(x$deg[1]), 
                 t.window=as.integer(x$win[2]), 
                 t.degree=as.integer(x$deg[2]), 
                 l.window=as.integer(x$win[3]), 
                 l.degree=as.integer(x$deg[3]), 
                 inner=x$inner, 
                 n.ahead=n.ahead, 
                 s.blend=s.blend, 
                 t.blend=t.blend, 
                 l.blend=l.blend, 
                 critfreq=critfreq,
                 arma=arma,
                 stats=stats)
}

remainder.stlop <- function(stlop, y) {
   as.vector(y - stlop$fit$O %*% y)
}

