"plot.variogram.cloud" <-
function (x, identify = FALSE, digitize = FALSE, 
	xlim = c(0, max(x$dist)), ylim = c(0, max(x$gamma)), 
	xlab = "distance", ylab = "semivariance", ...) 
{
    if (identify || digitize) {
        plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, ...)
        head = floor(x$np/2^16) + 1
        tail = floor(x$np%%2^16) + 1
		if (identify) {
        	labs = paste(head, tail, sep = ",")
        	sel = identify(x$dist, x$gamma, labs)
		} else {
			poly = locator(n = 512, type = "l")
			sel = point.in.polygon(x$dist, x$gamma, poly$x, poly$y)
		}
        ret = data.frame(cbind(head, tail)[sel, ])
        class(ret) = c("point.pairs", "data.frame")
        return(ret)
	} else {
        x$np = rep(1, length(x$gamma))
        plot.variogram(x, xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, ...)
    }
}
