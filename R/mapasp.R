mapasp <- function(data, x = data$x, y = data$y) {
	# calculates aspect ratio for levelplot of geographic data,
	# using proportial units (compare eqscplot)
	diff(range(y))/diff(range(x))
}
