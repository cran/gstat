# $Id: krige.cv.q,v 1.15 2006-05-15 16:26:54 edzer Exp $

if (!isGeneric("krige.cv"))
	setGeneric("krige.cv", function(formula, locations, ...)
		standardGeneric("krige.cv"))

krige.cv.locations = function (formula, locations, data = sys.frame(sys.parent()), 
	model = NULL, ..., beta = NULL, nmax = Inf, nmin = 0, maxdist = Inf, 
	nfold = nrow(data), verbose = FALSE) {

	gstat.cv(gstat(g = NULL, id = "var1", formula = formula, locations = 
		locations, data = data, model = model, beta = beta, nmax = nmax, 
		nmin = nmin, maxdist = maxdist, ...), nfold = nfold, verbose = verbose)
}
setMethod("krige.cv", c("formula", "formula"), krige.cv.locations)

krige.cv.spatial = function (formula, locations, model = NULL, ..., beta = NULL,
	nmax = Inf, nmin = 0, maxdist = Inf, nfold = nrow(locations), verbose = FALSE) {

	# data = locations 
	gstat.cv(gstat(g = NULL, id = "var1", formula = formula,
		data = locations, model = model, beta =
		beta, nmax = nmax, nmin = nmin, maxdist = maxdist,
		...), nfold = nfold, verbose = verbose)
}
setMethod("krige.cv", c("formula", "Spatial"), krige.cv.spatial)
