"plot.point.pairs" <-
function(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate", 
	ylab = "y coordinate", col.line = 2, ...) {
	xyplot(ycol ~ xcol, aspect = mapasp(x = xcol, y = ycol), 
		panel = panel.point.pairs, xlab = xlab, ylab = ylab, pairs = x,
		col.line = col.line, ...)
}
