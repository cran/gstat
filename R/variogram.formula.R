"variogram.formula" <-
function (object, locations, data, ...) 
{
#	if (missing(locations) && inherits(data, "spatial.data.frame"))
#		locations = sp.formula(data)
	ret = gstat.formula(object, locations, data)
	variogram(object = ret$y, locations = ret$locations, X = ret$X, ...)
}
