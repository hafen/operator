
predict.op <- function(op, y, se=FALSE, newdata=op$at, interval = c("none", "confidence", "prediction"), level = 0.95) {
   n <- ncol(op$O)
   if(length(y) != n)
      stop("must supply data (y) with same length as design 1, ..., ", n)

   interval <- match.arg(interval)

   # default behavior is to only predict ahead if
   # prediction intervals are specified
   if(interval=="prediction" & missing(newdata) & max(op$at > n)) {
      newdata <- op$at[op$at > n]
   }

   where <- op$at %in% newdata
   
   possible <- newdata %in% op$at
   if(any(!possible)) {
      warning("Some of newdata is not possible to obtain prediction for.  Make sure that when you obtain the operator matrix, you use the n.ahead option.")
   }
   
   fit <- op$O %*% y
   
   if(se || interval != "none") {
      if(length(op$stats) == 0) {
         stats <- getstats(op)
      } else {
         stats <- op$stats
      }
      
      start <- which(op$at == 1)
      ind <- start:(n + start - 1)
      
      residual.scale <- sqrt(sum((fit[ind] - y)^2)/stats$delta1)
      se.fit <- sqrt(op$var) * residual.scale
      se.pred <- sqrt(1 + op$var) * residual.scale
      df <- stats$delta1^2/stats$delta2

      fits <- data.frame(at=op$at[where], fit=as.vector(fit[where]), se.fit=se.fit[where])

      if (interval != "none") {
         if(level > 0 || level < 1) {
            aa <- (1 - level)/2

            if(interval == "confidence") {
               upper <- fit - qt(aa, df) * se.fit
               lower <- fit + qt(aa, df) * se.fit
            } else { # prediction intervals

               if (any(newdata < n))
                  warning("Predictions on already observed data refer to _future_ responses, which does not make sense for time series \n")

               upper <- fit - qt(aa, df) * se.pred
               lower <- fit + qt(aa, df) * se.pred
            }
            fits$lower <- lower[where]
            fits$upper <- upper[where]
         } else {
            warning("confidence level must be between 0 and 1, not calculated")
         }         
      }
      
      yy <- rep(NA, length(op$at))
      yy[ind] <- y
      data <- data.frame(x=op$at, y=yy, fit=fit)
      
      ret <- list(data=data, fits=fits, residual.scale=residual.scale, df=df, interval=interval, level=level)
      class(ret) <- "opPred"
      
      return(ret)
   } else {
      return(as.vector(fit[where]))
   }
}

getstats <- function(op) {
   if(!any(class(op) %in% c("op", "opblend")))
      stop("statistics calculation needs object of class 'op' or 'opblend'")
   
   n <- ncol(op$O)
   start <- which(op$at == 1)
   ind <- start:(n + start - 1)

   enp <- sum(diag(t(op$O[ind, 1:n]) %*% op$O[ind, 1:n]))
   Lbar <- diag(rep(1, n)) -  op$O[ind, 1:n]
   Lbar2 <- t(Lbar) %*% Lbar
   Lbar3 <- Lbar2 %*% Lbar2
   delta1 <- sum(diag(Lbar2))
   delta2 <- sum(diag(Lbar3))
   df <- delta1^2 / delta2

   list(enp=enp, delta1=delta1, delta2=delta2, trace.hat=sum(diag(op$O[ind, 1:n])), lambda=Lbar2)
}

predict.stlop <- function(stlop, y, ...) {
   if("armafit" %in% names(stlop)) {
      predict.op(stlop$armafit, y=y, ...)
   } else {
      predict.op(stlop$fit, y=y, ...)
   }
}

# remainder.stlop <- function(stlop, y) {
#    if("armafit" %in% names(stlop)) {
#       y - predict.op(stlop$armafit, y=y)
#    } else {
#       y - predict.op(stlop$fit, y=y)
#    }
# }


plot.opPred <- function(pred, fcol="black", pcol="darkgray", start=1, CIalpha=0.3, panel=mypanel, intcol="#00009C", ...) {
   stopifnot(class(pred) == "opPred")
   stopifnot(pred$interval %in% c("confidence", "prediction"))
   
   startInd <- pred$fits$at >= start
   
	mypanel <- function(x, y, CIxbox, CIybox, fit, fitx, ...) {
      panel.xyplot(x, y, ...)
      panel.polygon(CIxbox, CIybox, border=NA, col=intcol, alpha=CIalpha)
      panel.lines(fitx, fit, col=fcol)
   }

   p <- xyplot(y + fit ~ x, data=pred$data,
		panel=mypanel,
      col=c(pcol, fcol),
      type=c("p", "l"), distribute.type=TRUE,
      CIxbox = c(pred$fits$at[startInd], rev(pred$fits$at[startInd])),
      CIybox = c(pred$fits$lower[startInd], rev(pred$fits$upper[startInd])),
      fit = pred$fits$fit[startInd],
      fitx = pred$fits$at[startInd],
      prepanel=function(x, y) list(ylim=c(
         max(c(y, pred$fits$upper[startInd]), na.rm=TRUE), 
         min(c(y, pred$fits$lower[startInd]), na.rm=TRUE))
      ),
      subset = x >= start,
      ...
   )
   p
}

cp <- function(object, y, sigmasq=1) UseMethod("cp")

# cp.default <- function(...) {
# }

cp.op <- function(object, y, sigmasq=1) {

   if(length(object$stats) == 0) {
      stats <- getstats(object)
   } else {
      stats <- object$stats
   }

   resid <- y - predict(object, y, newdata=1:ncol(object$O))
   sigmahat <- sqrt(sum(resid^2)/stats$delta1)
   cp1 <- sum(resid^2)
   cp2 <- - stats$delta1 + stats$enp

   cp <- cp1/sigmasq + cp2

   data.frame(df=stats$enp, cp=cp, sigmahat=sigmahat, d1=stats$delta1, RSS=cp1, span=L$pars$span, deg=L$pars$deg)
}

cp.loess <- function(object, y, sigmasq=1) {
   resid <- resid(object)
   sigmahat <- sqrt(sum(resid^2)/object$one.delta)
   cp1 <- sum(resid^2)
   cp2 <- - object$one.delta + object$enp

   cp <- cp1/sigmasq + cp2

   data.frame(df=object$enp, cp=cp, sigmahat=sigmahat, d1=object$one.delta, RSS=cp1, span=object$pars$span * object$n, deg=object$pars$degree)
}

cp.stlop <- function(object, y, sigmasq=1) {
   if("armafit" %in% names(object)) {
      cp.op(object$armafit, y, sigmasq)
   } else {
      cp.op(object$fit, y, sigmasq)
   }
}

callsummary <- function(x) {
   ret <- "no information"
   if("call" %in% names(x)) {
      mod <- substr(as.character(x$call)[1], 1, 3)
      
      if(mod == "stl") {
         ret <- paste(
"stl with n.p=", x$pars$n.p, "
                  seasonal smoothing: ", ifelse(x$pars$periodic, "periodic", paste("deg=", x$pars$s.degree, ", span=", x$pars$s.window, sep="")), "
                  trend smoothing: deg=", x$pars$t.degree, ", span=", x$pars$t.window, "
                  l.degree=", x$pars$l.degree, ", l.window=", x$pars$l.window, ", inner=", x$pars$inner, sep="")         
      } else if(mod == "loe") {
         ret <- paste("loess with degree=", x$par$deg,", span=", x$pars$span, sep="")
      }
   }
   ret
}

anova.op <- function(mod1, mod2, y) {

   stopifnot(class(mod1) %in% c("op", "opblend"))
   stopifnot(class(mod2) %in% c("op", "opblend"))
   stopifnot(ncol(mod1$O) == ncol(mod2$O))
   stopifnot(length(y) == ncol(mod1$O))
   
   if(length(mod1$stats) == 0) {
      stats1 <- getstats(mod1)
   } else {
      stats1 <- mod1$stats
   }

   if(length(mod2$stats) == 0) {
      stats2 <- getstats(mod2)
   } else {
      stats2 <- mod2$stats
   }

   # yhat1 <- predict(mod1, y)
   # yhat2 <- predict(mod2, y)
   # resid1 <- yhat1 - y
   # resid2 <- yhat2 - y
   # 
   # RSS1 <- sum(resid1^2)
   # RSS2 <- sum(resid2^2)
   y <- as.vector(y)
   RSS1 <- (t(y) %*% stats1$lambda) %*% y
   RSS2 <- (t(y) %*% stats2$lambda) %*% y

   # null model is the one with the larger RSS
   if(RSS1 < RSS2) {             # mod2 is null
      num <- RSS2 - RSS1
      denom <- RSS1
      lambdadiff <- stats2$lambda - stats1$lambda
      delta1 <- stats1$delta1
      delta2 <- stats1$delta2
      enp <- c(stats2$enp, stats1$enp); rss <- c(RSS2, RSS1)
      descr <- paste("[1] Null model: ", callsummary(mod2), "\n", "[2]  Alt model: ", callsummary(mod1), sep="")
   } else {                      # mod1 is null
      num <- RSS1 - RSS2
      denom <- RSS2
      lambdadiff <- stats1$lambda - stats2$lambda
      delta1 <- stats2$delta1
      delta2 <- stats2$delta2      
      enp <- c(stats1$enp, stats2$enp); rss <- c(RSS1, RSS2)
      descr <- paste("[1] Null model: ", callsummary(mod1), "\n", "[2]  Alt model: ", callsummary(mod2), sep="")
   }

   nu1 <- sum(diag(lambdadiff))
   nu2 <- sum(diag(lambdadiff %*% lambdadiff))

   Fhat <- c(NA, (num/nu1)/(denom/delta1))

   dfnum <- nu1^2/nu2
   dfden <- delta1^2/delta2
   p <- pf(Fhat, dfnum, dfden, lower.tail=FALSE)

   # ret <- list(p=p, F=Fhat, df1=nu1^2/nu2, df2=delta1^2/delta2, RSSn=RSSn, RSSa=RSSa)
      
   ret <- data.frame(ENP = round(enp, 2), RSS = rss, `F-value` = Fhat, 
       `Pr(>F)` = p, check.names = FALSE)
   attr(ret, "heading") <- paste(descr, "\n\n", "Analysis of Variance: numerator df ", format(round(dfnum, 2)), ", denominator df ", format(round(dfden, 2)), "\n", sep = "")
   class(ret) <- c("anova", "data.frame")
   ret
}


predict.armaOp <- function(armaop, n.ahead) {
	aop <- armaop$O
	armamod <- armaop$mod

	n <- length(resid(armamod))
	phi <- armamod$model$phi
	theta <- -armamod$model$theta

	if(length(theta)==1 && theta[1] == 0) theta <- numeric(0)

	np <- length(phi)
	nq <- length(theta)

	aop <- rbind(aop, matrix(nrow=n.ahead, ncol=n, data=0))

	for(i in 1:n.ahead) {
		ll <- rep(0, n)

		# AR part
		if(np > 0) {
			for(j in 1:np) {
				ind <- n + i - (j - 1) - 1
				if(ind <= n) {
					z <- rep(0, n)
					z[ind] <- phi[j]
					ll <- ll + z
				} else {
					ll <- ll + aop[ind,] * phi[j]
				}
			}
		}

		# MA part
		if(nq > 0 && i <= nq) {
			for(j in 1:nq) {
				ind <- n + i - (j - 1) - 1
				if(ind <= n) {
					z <- lll <- rep(0, n)
					z[ind] <- 1
					lll[1:ind] <- aop[ind, 1:ind]
					ll <- ll + (z - lll) * (-theta[j])
				}
			}
		}
		aop[n + i,] <- ll
	}
	aop
}



# closely mimics anova.loess
# but calculation of delta2 is not exact...
# anova.op <- function(mod1, mod2, y) {
# 
#    stopifnot(class(mod1) %in% c("op", "opblend"))
#    stopifnot(class(mod2) %in% c("op", "opblend"))
#    stopifnot(ncol(mod1$O) == ncol(mod2$O))
#    stopifnot(length(y) == ncol(mod1$O))
#    
#    objects <- list(mod1, mod2)
#    
#    nmodels <- 2
#    models <- lapply(objects, callsummary)
#    
#    delta1 <- sapply(objects, function(x) x$stats$delta1)
#    delta2 <- sapply(objects, function(x) x$stats$delta2)
#    enp <- sapply(objects, function(x) x$stats$enp)
#    max.enp <- order(enp)[nmodels]
# 
#    modnames <- c("Null model", "  Alt model")[order(enp)]
# 
#    descr <- paste(modnames, ": ", models, sep = "", collapse = "\n")
#    # descr <- paste("Model ", format(1:nmodels), ": ", models, 
#        # sep = "", collapse = "\n")
# 
#    rss <- sapply(objects, function(x) sum((y - predict(x, y))^2))
#    s <- sqrt(rss/delta1)
# 
#    d1diff <- abs(diff(delta1))
#    dfnum <- c(d1diff^2/abs(diff(delta2)))
#    dfden <- (delta1^2/delta2)[max.enp]
#    Fvalue <- c(NA, (abs(diff(rss))/d1diff)/s[max.enp]^2)
# 
#    pr <- pf(Fvalue, dfnum, dfden, lower.tail = FALSE)
#    ans <- data.frame(ENP = round(enp, 2), RSS = rss, `F-value` = Fvalue, 
#        `Pr(>F)` = pr, check.names = FALSE)
#    attr(ans, "heading") <- paste(descr, "\n\n", "Analysis of Variance: numerator df ", format(round(dfnum, 2)), ", denominator df ", 
#        format(round(dfden, 2)), "\n", sep = "")
#    class(ans) <- c("anova", "data.frame")
#    ans
# }
