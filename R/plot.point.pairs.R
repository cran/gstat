"plot.point.pairs" <-
function(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate", 
ylab = "y coordinate", ...) {
	xyplot(ycol ~ xcol, aspect = mapasp(x = xcol, y = ycol), 
		panel = panel.point.pairs, xlab = xlab, ylab = ylab, pairs = x,
		...)
}
