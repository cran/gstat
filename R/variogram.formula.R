"variogram.formula" <-
function (object, locations, data, ...) 
{
	ret = gstat.formula(object, locations, data)
	variogram(object = ret$y, locations = ret$locations, X = ret$X, ...)
}
