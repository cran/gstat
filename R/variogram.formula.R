"variogram.formula" <-
function (object, locations = try.coordinates(data), data, ...) 
{
	# gstat.formula takes care of the case where locations contains
	# both data and coordinates --- see there.
	## ret = gstat.formula(object, locations, data)
	## variogram(object = ret$y, locations = ret$locations, X = ret$X, ...)
	g = gstat(formula = object, locations = locations, data = data)
	variogram(g, ...)
}
