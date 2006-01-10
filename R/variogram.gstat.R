"variogram.gstat" = function (object, ...) {
	if (!inherits(object, "gstat"))
		stop("first argument should be of class gstat")
	y = list()
	locations = list()
	X = list()
	beta = list()
	grid = list()
	for (i in seq(along = object$data)) {
		d = object$data[[i]]
		raw = gstat.formula(d$formula, d$data)
		y[[i]] = raw$y
		locations[[i]] = raw$locations
		X[[i]] = raw$X
		beta[[i]] = raw$beta
		grid[[i]] = raw$grid
		if (d$degree != 0)
			stop("degree != 0: residual variograms wrt coord trend using degree not supported")
	}
	names(y) = names(locations) = names(X) = names(object$data)
	# call variogram.default() next:
	variogram(y, locations, X, trend.beta = beta, grid = grid, g = object, ...)
}
