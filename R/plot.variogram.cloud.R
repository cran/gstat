"plot.variogram.cloud" <-
function (x, identify = FALSE, xlim = c(0, max(x$dist)), 
	ylim = c(0, max(x$gamma)), xlab = "distance", ylab = "semivariance", ...)
{
    if (identify) {
        plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
			ylab = ylab, ...)
        head = floor(x$np/2^16) + 1
        tail = floor(x$np%%2^16) + 1
        labs = paste(head, tail, sep = ",")
        sel = identify(x$dist, x$gamma, labs)
        return(cbind(head, tail)[sel, ])
    } else {
		x$np = rep(1, length(x$gamma))
		plot.variogram(x, xlim = xlim, ylim = ylim, xlab = xlab, 
			ylab = ylab, ...)
	}
}
