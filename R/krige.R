"krige" <-
function (formula, locations, data = sys.frame(sys.parent()), 
	newdata, model = NULL, beta = NULL, nmax = Inf, nmin = 0, 
	maxdist = Inf, block = numeric(0), nsim = 0, indicators = FALSE, 
	na.action = na.pass, ...)
{
#	if (missing(locations) && inherits(data, "spatial.data.frame"))
#		locations = sp.formula(data)
    g = gstat(formula = formula, locations = locations, model = model,
		data = data, beta = beta, nmax = nmax, nmin = nmin, 
		maxdist = maxdist, ...)
    predict.gstat(g, newdata = newdata, block = block, nsim = nsim,
		indicators = indicators, na.action = na.action)
}
