# get the loess filter coefficients
filterCoef <- function(win, deg) {
   win2 <- (win - 1)/2
   d <- loessOp(2*win, span=win, degree=deg, at=win2+1, stats=FALSE)$O[1, 1:win]
   d
}

freqr <- function(...) UseMethod("freqr")

freqr.default <- function(n.p, s.window, s.degree=0, t.window=NULL, t.degree=1, l.window=nextodd(n.p), l.degree=t.degree, periodic=FALSE, M=400) {

   N <- 2*M
   fk <- c(0:M)/N

   nextodd <- function(x) {
      x <- round(x)
      x <- ifelse(x%%2==0, x+1, x)
      as.integer(x)
   }

   deg.check <- function(deg) {
      degname <- deparse(substitute(deg))
      deg <- as.integer(deg)
      if (deg < 0 || deg > 2)
         stop(degname, " must be 0, 1, or 2")
      deg
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

	# check parameters
   s.degree <- deg.check(s.degree)
   t.degree <- deg.check(t.degree)
   l.degree <- deg.check(l.degree)

   if (is.null(t.window)) 
       t.window <- nextodd(ceiling(1.5 * n.p/(1 - 1.5/s.window)))

   s.window <- nextodd(s.window)
   t.window <- nextodd(t.window)
   l.window <- nextodd(l.window)

   t.window2 <- (t.window - 1)/2
   s.window2 <- (s.window - 1)/2
   l.window2 <- (l.window - 1)/2

   dt <- filterCoef(t.window, t.degree)

   jj <- -t.window2:t.window2
   tt <- apply(cos(2*pi*outer(fk, jj)), 1, function(x) sum(dt*x))
   dtt <- data.frame(fk=fk, frfv=tt, which="Trend Filter")

   if(!periodic) {
      ds <- filterCoef(s.window, s.degree)

      kk <- -s.window2:s.window2
      cc <- apply(cos(2*pi*outer(fk*n.p, kk)), 1, function(x) sum(ds*x))
      dcc <- data.frame(fk=fk, frfv=cc, which="Cycle-Subseries Filter")

   } else { # if periodic, spikes at (1/n.p and harmonics)

      dfk <- 1/N
      ind <- (1/n.p)/dfk
      mind <- floor((M + 1)/ind)
      ind <- c(1, round(ind*c(1:mind)+ 1))

      cc <- rep(0, M + 1)
      cc[ind] <- 1
      # TODO: add warning here it it looks like grid isn't fine enough to accurately show location of spikes

      dcc <- data.frame(fk=fk, frfv=cc, which="Cycle-Subseries Filter")      
   }

   dl <- filterCoef(l.window, l.degree)

   ## L is a convolution of two moving averages of length 3, 
   # a moving average of length 12
   # and loess with l.window, l.degree
   dl <- convolve(
      convolve(
         convolve(rep(1/3, 3), rev(rep(1/n.p, n.p)), type="open"),
         rev(rep(1/n.p, n.p)), type="open"
      ),
      rev(dl), type="open"
   )

   l.window2 <- (length(dl) - 1)/2
   ii <- -l.window2:l.window2
   ll <- apply(cos(2*pi*outer(fk, ii)), 1, function(x) sum(dl*x))
   dhh <- data.frame(fk=fk, frfv=1-ll, which="High Pass Filter")

   dss <- data.frame(fk=fk, frfv=(1 - ll)*cc, which="Seasonal Filter")

   all <- rbind(dtt, dcc, dhh, dss)

   class(all) <- c("freqr", "data.frame")
   all
}

freqr.stl2 <- function(x, ...) {
   freqr.default(n.p=x$pars$n.p, 
     s.window=x$pars$win$s.window, 
     s.degree=x$pars$deg$s.degree, 
     t.window=x$pars$win$t.window, 
     t.degree=x$pars$deg$t.degree, 
     l.window=x$pars$win$l.window, 
     l.degree=x$pars$deg$l.degree, 
     periodic=x$pars$periodic, ...)
}

freqr.stl <- function(x, ...) {
   periodic <- FALSE
   
   # strange way to check if seasonal is periodic
   if(length(x$time.series[,1]) * 10 + 1 == x$win[1] && x$deg[1] == 0) {
      periodic <- TRUE
   }

   freqr.default(n.p=frequency(x$time.series[,1]), 
     s.window=as.integer(x$win[1]), 
     s.degree=as.integer(x$deg[1]), 
     t.window=as.integer(x$win[2]), 
     t.degree=as.integer(x$deg[2]), 
     l.window=as.integer(x$win[3]), 
     l.degree=as.integer(x$deg[3]), 
     periodic=periodic, ...)
}

freqr.stlop <- function(x, ...) {
   freqr.default(n.p=x$pars$n.p,
      s.window=x$pars$s.window,
      s.degree=x$pars$s.degree,
      t.window=x$pars$t.window,
      t.degree=x$pars$t.degree,
      l.window=x$pars$l.window,
      l.degree=x$pars$l.degree, 
      periodic=x$pars$periodic, ...)
}

# TODO: better way to get n.p to getcrit? - maybe store it as an attribute to freqr?
plot.freqr <- function(x, critfreq=NA, xlab="Frequency", ylab="Transfer Function", layout=c(1, length(unique(x$which))), type=c("g", "l"), as.table=TRUE, panel=freqrPanel, n.p=NULL, ...) {
   xyplot(abs(frfv)^2 ~ fk | which, data=x,
      xlab=xlab,
      ylab=ylab,
      layout=layout,
      type=type,
      as.table=as.table,
      panel=panel,
      critfreq=critfreq,
      dat=x,
      n.p=n.p,
      ...
   )
}

# TODO: make it so lim only gets used if it is an stl decomposition
freqrPanel <- function(x, y, critfreq=0.05, dat, n.p, ...) {
   vlines <- NULL

   # add in vertical lines if lim is specified
   if(is.numeric(critfreq)) {
      if(is.null(n.p)) {
         warning("must specify n.p")
      } else {
         lms <- getcrit(dat, lim=critfreq, n.p=n.p)         
      
         pn <- panel.number()
         if(pn == 1) {
            vlines <- c(lms$ft_low)
         } else if(pn == 2) {
            vlines <- c(lms$fc_low, lms$fc_upp)
         } else if(pn == 3) {
            vlines <- c(lms$fh_low, lms$fh_upp)
         }
      }
   }
   panel.xyplot(x, y, ...)
   if(!is.null(vlines))
      panel.abline(v=vlines, lty=2, col="black")
}

getcrit <- function(a, lim=0.05, n.p) {
   if(!("freqr" %in% class(a)))
      stop("object is not of class 'freqr'")

	ff <- subset(a, which=="Trend Filter")$fk
	tt <- subset(a, which=="Trend Filter")$frfv^2
	cc <- subset(a, which=="Cycle-Subseries Filter")$frfv^2
	hh <- subset(a, which=="High Pass Filter")$frfv^2

	tlow <- min(which(tt < lim))
	ft_low <- approxfun(c(tt[tlow-1], tt[tlow]), c(ff[tlow-1], ff[tlow]))(lim)

	clow <- min(which(cc < lim))
	fc_low <- approxfun(c(cc[clow-1], cc[clow]), c(ff[clow-1], ff[clow]))(lim)
	fc_upp <- 1/n.p - fc_low

	hh2 <- hh[1:which.min(ifelse(diff(hh)<0, 0, 1))]
	hlow <- max(which(hh2 < lim))
	hmid <- max(which(hh2 < 0.5))
	hupp <- max(which(hh2 < 1 - lim))

	fh_low <- approxfun(c(hh[hlow], hh[hlow+1]), c(ff[hlow], ff[hlow+1]))(lim)
	fh_mid <- approxfun(c(hh[hmid], hh[hmid+1]), c(ff[hmid], ff[hmid+1]))(0.5)
	fh_upp <- approxfun(c(hh[hupp], hh[hupp+1]), c(ff[hupp], ff[hupp+1]))(1-lim)

	data.frame(ft_low=ft_low, fc_low=fc_low, fc_upp=fc_upp, fh_low=fh_low, fh_mid=fh_mid, fh_upp=fh_upp)
}


# now just for loess
freqrlo <- function(...) UseMethod("freqrlo")

freqrlo.default <- function(span, degree, M=400, at="symmetric") {
	fk <- c(0:M)/(2*M)

   span2 <- (span - 1)/2

   if(at == "symmetric")
      at <- span2 + 1

   if(at > span2 + 1) {
      at <- span2 + 1
      warning("specified value of 'at' taken to mean 'symmetric'")
   }

   if(at < span2) {
      ind <- c(1:span)
   } else {
      ind <- (at - span2):(at + span2)
   }

   dt <- loessOp(2*span, span=span, degree=degree, at=at, stats=FALSE)$O[1, ind]

   # dt <- filterCoef(span, degree)
   jj <- -span2:span2
   fr <- apply(exp(-2i*pi*outer(fk, jj)), 1, function(x) sum(dt*x))
	
	dfr <- data.frame(fk=fk, frfv=fr, which="Loess Filter")
   class(dfr) <- c("freqr", "data.frame")
   dfr
}

freqrlo.op <- function(op, M=400, at=1) {
	fk <- c(0:M)/(2*M)
   
   n <- ncol(op$O)
   dt <- op$O[at,]
   jj <- 1:n - at
   fr <- apply(exp(-2i*pi*outer(fk, jj)), 1, function(x) sum(dt*x))

   dfr <- data.frame(fk=fk, frfv=fr, which=paste("Filter at ", at, sep=""))
   class(dfr) <- c("freqr", "data.frame")
   dfr
}

freqrlo.loess <- function(x, M=400, at="symmetric") {
   nextodd <- function(x) {
      x <- round(x)
      x <- ifelse(x%%2==0, x+1, x)
      as.integer(x)
   }
   
   sp <- as.integer(x$pars$span * x$n)
   sp <- nextodd(sp)
   deg <- x$pars$degree
   freqrlo.default(sp, deg, M, at)
}

# ## A beginning attempt to add alternative axis labels (period scale) on the top
axis.period <- function(side,...) {
   if(side=="top") {
      # row <- lattice.getStatus("current.focus.row")
      # column <- lattice.getStatus("current.focus.column")
      # panel.layout <- trellis.currentLayout("panel")

      if(packet.number()==1)
      panel.axis(side=side, at=1/c(12, 7, 4, 3, 2), 
         labels=paste("1/", c(12, 7, 4, 3, 2), sep=""), outside=TRUE, rot=0
      ) 
   }
   else
      axis.default(side=side,...) 
} 

# a <- freqrlo(11, 2)
# plot(a, strip=F, axis=axis.period)





