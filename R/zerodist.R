"zerodist" <- function(x, y, z, zero = 0.0) {
	zero = zero*zero # work with squares
	# calculates matrix with pairwise distances for 
	# coordinate vectors x and y:
	D <- outer(x, x, "-")^2 
	diag(D) <- 1
	if (!any(D <= zero))
		return(numeric(0))
	if (!missing(y))
		D <- D + outer(y, y, "-")^2
	diag(D) <- 1
	if (!any(D <= zero))
		return(numeric(0))
	if (!missing(z))
		D <- D + outer(z, z, "-")^2
	diag(D) <- 1
	n <- length(x)
	index <- 1:(n*n)
	z <- index[as.vector(D) <= zero]
    ret <- cbind(((z - 1) %/% n) + 1, ifelse(z %% n == 0, n, z %% n))
	if (zero > 0) {
		ret = cbind(ret, sqrt(as.vector(D)[z]))
    	ret = matrix(ret[ret[, 1] < ret[, 2], ], ncol = 3)
		colnames(ret) = c("left", "right", "dist")
	} else {
    	ret = matrix(ret[ret[, 1] < ret[, 2], ], ncol = 2)
		colnames(ret) = c("left", "right")
	}
	ret
}
