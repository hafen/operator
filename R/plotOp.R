############################################################################
### levelplot of the operator matrix
############################################################################

plotOp <- function(op, nbreaks=30, panel=op.panel, aspect="xy", xlab="Column", ylab="Row", ...) {
   M <- op$O
   n <- ncol(M)
   require(lattice)
   m.dim <- dim(M)

   op.panel <- function(x, y, z, ...) {
      panel.levelplot(x, y, z, ...)
      if(op$at[1]<1)
         panel.abline(h=1, lty=2, col="black")
      if(tail(op$at, 1) > n)
         panel.abline(h=n, lty=2, col="black")
   }
   
   if(is.null(m.dim))
      stop("Must supply a matrix!")

   if(length(m.dim) != 2)
      stop("Matrix must have dimension 2")

   # x is rows
   # y is columns
   gg  <- expand.grid(x=op$at, y=1:m.dim[2])
   # as.vector lists the matrix by row
   gg$z <- as.vector(M)
   gg.breaks <- do.breaks(range(gg$z), nbreaks)
   # shift the breaks
   gg.d <- gg.breaks[2] - gg.breaks[1]
   gg.breaks <- gg.breaks - gg.d/2
   gg.breaks <- c(gg.breaks, gg.breaks[nbreaks+1] + gg.d)
   # find where 0 occurs (want 0 to be white)
   gg.zero <- which.min(abs(gg.breaks))
   gg.zero <- gg.zero + ifelse(gg.breaks[gg.zero] < 0, 0.5, -0.5)
   # set colors so that 0 should be white
   redwhite.pal <- colorRampPalette(c("red", "white"), space = "rgb")
   whiteblue.pal <- colorRampPalette(c("white", "blue", "orange"), space = "rgb")
   col.regions <- c(redwhite.pal((gg.zero-2)*30), whiteblue.pal((length(gg.breaks) + 2 - gg.zero)*30))

   p <- levelplot(z ~ y*x, data=gg, xlab=xlab, ylab=ylab, at=gg.breaks, col.regions=col.regions,
      panel=panel,
      aspect=aspect,
      ...
   )
   p
}

