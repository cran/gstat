"plot.variogramCloud" <-
function (x, identify = FALSE, digitize = FALSE, 
	xlim = c(0, max(x$dist)), ylim = c(0, max(x$gamma)), 
	xlab = "distance", ylab = "semivariance", keep = FALSE, ...) 
{
    if (identify || digitize) {
        plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, ...)
        head = floor(x$np/2^16) + 1
        tail = floor(x$np%%2^16) + 1
		if (identify) {
			print("mouse-left identifies, mouse-right stops")
        	labs = paste(head, tail, sep = ",")
        	sel = identify(x$dist, x$gamma, labs, pos = keep)
		} else {
			print("mouse-left digitizes, mouse-right closes polygon")
			poly = locator(n = 512, type = "l")
			if (!is.null(poly))
				sel = point.in.polygon(x$dist, x$gamma, poly$x, poly$y)
			else stop("digitized selection is empty")
		}
		ret = data.frame(cbind(head, tail)[sel, ])
		class(ret) = c("pointPairs", "data.frame")
        if (keep) {
			if (identify) {
				attr(x, "sel") = sel
				attr(x, "text") = labs[sel$ind]
			} else  # digitize
				attr(x, "poly") = poly
			attr(x, "ppairs") = ret
			return(x)
		} else 
        	return(ret)
	} else {
		sel = attr(x, "sel")
		lab = attr(x, "text")
		poly = attr(x, "poly")
		if (!is.null(sel) && !is.null(lab)) {
        	plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            	ylab = ylab, ...)
			text(x$dist[sel$ind], x$gamma[sel$ind], labels=lab, pos= sel$pos)
		} else if (!is.null(poly)) {
        	plot(x$dist, x$gamma, xlim = xlim, ylim = ylim, xlab = xlab, 
            	ylab = ylab, ...)
			lines(poly$x, poly$y)
		} else {
        	x$np = rep(1, length(x$gamma))
        	plot.gstatVariogram(x, xlim = xlim, ylim = ylim, xlab = xlab, 
        	    ylab = ylab, ...)
		}
    }
}
