mapasp <- function(data, x = data$x, y = data$y, get.number = FALSE) {
	# calculates aspect ratio for levelplot of geographic data,
	# using proportial units (compare eqscplot)
	if (version$major >= 2 && !get.number)
		return("iso")
	else
	diff(range(y))/diff(range(x))
}
