"plot.variogramCloud" <-
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
			print("mouse-left identifies, mouse-right stops")
        	labs = paste(head, tail, sep = ",")
        	sel = identify(x$dist, x$gamma, labs)
		} else {
			print("mouse-left digitizes, mouse-right closes polygon")
			poly = locator(n = 512, type = "l")
			if (!is.null(poly))
				sel = point.in.polygon(x$dist, x$gamma, poly$x, poly$y)
			else stop("digitized selection is empty")
		}
        ret = data.frame(cbind(head, tail)[sel, ])
        class(ret) = c("pointPairs", "data.frame")
        return(ret)
	} else {
        x$np = rep(1, length(x$gamma))
        plot.gstatVariogram(x, xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, ...)
    }
}
