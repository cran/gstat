"select.spatial" <-
function(x = data$x, y = data$y, data, pch = "+", n = 512) {
	plot(x, y, pch = pch, asp = 1)
	pol = locator(n = n, type = "o")
	which(point.in.polygon(x, y, pol$x, pol$y))
}
