armaOp <- function(armamod) {
	n <- length(resid(armamod))
	
	phi <- armamod$model$phi
	theta <- -armamod$model$theta

	p <- length(phi)
	q <- length(theta)

	pq <- max(p, q)

	if(p < q) {
		phi <- c(phi, rep(0, q - p))
	} else if(q < p) {
		theta <- c(theta, rep(0, p - q))
	}

	pii <- rep(0, n - 1)

	pii[1] <- theta[1] - phi[1]

	for(j in 2:(n-1)) {
		for(i in 1:min(pq, j - 1)) {
			pii[j] <- pii[j] + theta[i]*pii[j - i]
		}
		if(j <= pq)
			pii[j] <- pii[j] + theta[j] - phi[j]
	}
	pii <- c(0, -pii)

	O <- matrix(nrow=n, ncol=n, data=0)
	
	for(i in 2:n) {
		O[i, 1:i] <- rev(pii[1:i])
	}
	ret <- list(
	 O = O,
	 mod = armamod   
	)
   class(ret) <- "armaOp"
   ret
}


# # new way, where theta_i=0 corresponding to e_1, ..., e_p
# 
# if(p > 0) {
#  for(i in 1:p) {
#     O[i,i] <- 1
#  }     
# }
# 
# for(i in (p+1):n) {
# 
#  # now need to recalculate pi each time
#  pii <- rep(0, i - 1)
# 
#  pii[1] <- phi[1] - theta[1]
#    tmp <- phi[1]
# 
#  if(i > 2) {
#     for(j in 2:(i-1)) {
#          # end <- ifelse(i - j < p, TRUE, FALSE)
# 
#        if(j <= p) {
#           for(k in 1:min(j, q)) {
#              pii[j] <- pii[j-k] * theta[k] + phi[j] 
#           }              
#        } else {
#           for(k in 1:min(j, q)) {
#              pii[j] <- pii[j] + end * theta[k] * pii[j-k]
#           }
#        }
# 
#     }     
#  }
#  pii <- c(0, pii)
# 
# x[3] - (phi[1] * x[2] + phi[2] * x[1])
# 
# x[4] - ((phi[1] - theta[1])*x[3] + (phi[2] + theta[1]*phi[1])*x[2] + theta[1]*phi[2]*x[1])
# 
# phi[2] + theta[1]*phi[1] - theta[1]^2
