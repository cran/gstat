"select.spatial" <-
function(x = data$x, y = data$y, data, pch = "+") {
	require(MASS)
	eqscplot(x, y, pch = pch)
	pol = locator(n = 512, type = "o")
	sel = 1:length(x)
	sel[point.in.polygon(x, y, pol$x, pol$y)]
}
