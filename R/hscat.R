hscat = function(formula, data, breaks, ...) {
	stopifnot(!missing(breaks))
	x = variogram(formula, data, cloud = TRUE, ...)
	.BigInt = attr(x, ".BigInt")
	x$left = x$np%%.BigInt + 1
	x$right = x$np%/%.BigInt + 1
	x$class = cut(x$dist, breaks = breaks)
	y = model.frame(formula, data)[[1]]
	x$xx = y[x$left]
	x$yy = y[x$right]
	lab = as.character(formula)[2]
	panel = function(x,y,subscripts, ...) {
		xr = c(min(x),max(x))
		llines(xr, xr)
		lpoints(x,y,...)
		ltext(min(x), max(y), paste("r =", signif(cor(x,y),3)), adj=c(0,0.5))
	}
	xyplot(xx~yy|class, x, panel = panel,
		main = "lagged scatterplots", xlab = lab, ylab = lab, ...)
}

