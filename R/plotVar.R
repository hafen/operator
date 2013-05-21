plotVar <- function(op, ylab="variance", xlab="x", panel=var.panel, auto.key=list(lines=TRUE, points=FALSE, columns=3), ...) {

	require(lattice)

   var.panel <- function(x, y, ...) {
      # panel.polygon(c(op$at[1]-1000, op$n.b, op$n.b, op$at[1]-1000), c(0, 0, 100, 100), border=NA, fill="lightgray")
      panel.grid(h=-1, v=FALSE)
      if(op$at[1]<1)
         panel.abline(v=1, lty=2, col="black")
      if(tail(op$at, 1) > n)
         panel.abline(v=n, lty=2, col="black")
      panel.xyplot(x, y, ...)
   }

   if(class(op) == "opblend") {
      n <- ncol(op$O)
      
      d <- with(op, make.groups(
         from.var = from.var,
         to.var = to.var,
         blend.var = var)
      )
      d$x <- rep(op$at, 3)
      
      lab <- paste(round(op$blend.p*100), "% blend from deg=", op$from.deg, " span=", op$from.span, " to deg=", op$to.deg, " span=", op$to.span, sep="")
      p <- xyplot(data ~ x | lab, groups=which, data=d, type="l", 
         panel=panel,
         xlab=xlab,
         ylab=ylab,
         auto.key=auto.key,
         ...
      )
   } else if(class(op) == "stlop") {
      n <- ncol(op$fit$O)

      if(op$pars$fc.number > 0) {
         d <- with(op, make.groups(
            "Total variance" = op$fit$var,
            "Seasonal variance" = op$seas$var
         ))
         for(i in 1:op$pars$fc.number) {
            d <- rbind(d, data.frame(data=opVar(op$fc[[i]]), which=op$pars$fc.name[i]))
         }
         d$x <- rep(op$seas$at, 2 + op$pars$fc.number)
      } else {
         d <- with(op, make.groups(
            "Total variance" = op$fit$var,
            "Seasonal variance" = op$seas$var,
            "Trend variance" = op$trend$var)
         )
         d$x <- rep(op$seas$at, 3)
      }
      
      p <- xyplot(data ~ x, groups=which, data=d, type="l", 
         panel=panel,
         xlab=xlab,
         ylab=ylab,
         auto.key=auto.key,
         ...
      )
   } else {
      n <- ncol(op$O)
      lab <- ""
      if(class(op) == "op")
         lab <- paste("deg=", op$pars$deg, " span=", op$pars$span, sep="")
         
      # sometimes all are equal ("periodic"), but not at several decimal places...
      # if(all.equal(op$var, rep(op$var[1], length(op$var))))
      #    yy <- rep(op$var[1], length(op$var))
      # else
         yy <- op$var

      p <- xyplot(yy ~ op$at | lab, panel=panel, type="l", xlab=xlab, ylab=ylab, ...)
   }
   p
}

