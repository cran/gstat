"plot.pointPairs" <-
function(x, data, xcol = data$x, ycol = data$y, xlab = "x coordinate", 
	ylab = "y coordinate", col.line = 2, line.pch = 0, ...) {
	xyplot(ycol ~ xcol, aspect = mapasp(x = xcol, y = ycol), 
		panel = panel.pointPairs, xlab = xlab, ylab = ylab, pairs = x,
		col.line = col.line, line.pch = line.pch, ...)
}
