if (!isGeneric("krige"))
	setGeneric("krige", function(formula, locations, ...)
		standardGeneric("krige"))

"krige.locations" <-
function (formula, locations, data = sys.frame(sys.parent()), 
	newdata, model = NULL, ..., beta = NULL, nmax = Inf, nmin = 0, 
	maxdist = Inf, block = numeric(0), nsim = 0, indicators = FALSE, 
	na.action = na.pass)
{
    g = gstat(formula = formula, locations = locations, data = data, 
		model = model, beta = beta, nmax = nmax, nmin = nmin, 
		maxdist = maxdist, ...)
    predict.gstat(g, newdata = newdata, block = block, nsim = nsim,
		indicators = indicators, na.action = na.action)
}
setMethod("krige", c("formula", "formula"), krige.locations)


krige.spatial <- function(formula, locations, newdata, model = NULL, ..., 
	beta = NULL, nmax = Inf, nmin = 0, maxdist = Inf, block = numeric(0), 
	nsim = 0, indicators = FALSE, na.action = na.pass)
{
	# locations = coordinates(arg2)
    g = gstat(formula = formula, # locations = locations, 
		data = locations, 
		model = model, beta = beta, nmax = nmax, nmin = nmin, 
		maxdist = maxdist, ...)
    predict.gstat(g, newdata = newdata, block = block, nsim = nsim,
		indicators = indicators, na.action = na.action)
}
setMethod("krige", c("formula", "Spatial"), krige.spatial)
setMethod("krige", c("formula", "NULL"), krige.spatial)

if (!isGeneric("idw"))
	setGeneric("idw", function(formula, locations, ...)
		standardGeneric("idw"))

idw.locations <-
function (formula, locations, data = sys.frame(sys.parent()), 
		newdata, nmax = Inf, nmin = 0, maxdist = Inf, block = numeric(0), 
		na.action = na.pass, idp = 2.0) {
	krige(formula, locations, data, newdata, nmax = nmax, nmin = nmin,
		maxdist = maxdist, block = block, na.action = na.action,
		set = list(idp = idp))
}
setMethod("idw", c("formula", "formula"), idw.locations)

idw.spatial <-
function (formula, locations, 
		newdata, nmax = Inf, nmin = 0, maxdist = Inf, block = numeric(0), 
		na.action = na.pass, idp = 2.0) {
	krige(formula, locations, newdata, nmax = nmax, nmin = nmin,
		maxdist = maxdist, block = block, na.action = na.action,
		set = list(idp = idp))
}
setMethod("idw", c("formula", "Spatial"), idw.spatial)
