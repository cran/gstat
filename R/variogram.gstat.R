"variogram.gstat" <-
function (object, ...) {
	if (!inherits(object, "gstat"))
		stop("first argument should be of class gstat")
	y = list()
	locations = list()
	X = list()
	beta = list()
	for (i in 1:length(object$data)) {
		d = object$data[[i]]
		raw = gstat.formula(d$formula, d$locations, d$data)
		y[[i]] = raw$y
		locations[[i]] = raw$locations
		X[[i]] = raw$X
		beta[[i]] = raw$beta
	}
	names(y) = names(locations) = names(X) = names(object$data)
	variogram.default(y, locations, X, trend.beta = beta, ...)
}
