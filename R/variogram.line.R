"variogram.line" <-
function(object, maxdist, n=200, min=1.0e-6*maxdist, ids=c(0,0), 
	dir = c(1,0,0), ...)
{
	if (missing(object))
		stop("model is missing");
	if (!inherits(object, "variogram.model"))
		stop("model should be of mode variogram.model (use function vgm)")
	if (missing(maxdist))
		stop("maxdist is missing");
	if (length(dir) != 3)
		stop("dir should be numeric vector of length 3")
	pars = c(min,maxdist,n,dir)
	load.variogram.model(object, ids) # loads object into gstat 
	ret = .Call("gstat_variogram_values", as.integer(ids),
		as.numeric(pars), PACKAGE = "gstat")
	.Call("gstat_exit", 0, PACKAGE = "gstat");
	data.frame(dist=ret[[1]], gamma=ret[[2]])
}
