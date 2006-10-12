# $Id: variogramLine.q,v 1.3 2006-10-31 14:36:55 edzer Exp $

"variogramLine" <-
function(object, maxdist, n=200, min=1.0e-6 * maxdist, dir = c(1,0,0), 
	covariance = FALSE,	..., debug.level = 0)
{
	if (missing(object))
		stop("model is missing");
	if (!inherits(object, "variogramModel"))
		stop("model should be of mode variogramModel (use function vgm)")
	if (missing(maxdist))
		stop("maxdist is missing");
	if (length(dir) != 3)
		stop("dir should be numeric vector of length 3")
	.Call("gstat_init", as.integer(debug.level))
	pars = c(min,maxdist,n,dir)
	load.variogram.model(object, c(0,0)) # loads object into gstat 
	ret = .Call("gstat_variogram_values", as.integer(c(0,0)),
		as.numeric(pars), as.integer(covariance))
	.Call("gstat_exit", 0);
	data.frame(dist=ret[[1]], gamma=ret[[2]])
}

"variogram.line" <- function(..., deprecate = FALSE) {
	if (deprecate)
		cat("variogram.line is DEPRECATED, please use variogramLine instead\n")
	variogramLine(...)
}
