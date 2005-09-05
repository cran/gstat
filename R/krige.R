"krige" <-
function (formula, locations = try.coordinates(data), data = sys.frame(sys.parent()), 
	newdata, model = NULL, ..., beta = NULL, nmax = Inf, nmin = 0, 
	maxdist = Inf, block = numeric(0), nsim = 0, indicators = FALSE, 
	na.action = na.pass)
{
	if (has.coordinates(locations)) { # shift arguments:
		if (!is(data, "Spatial")) # another shift:
			stop("if data derives from Spatial, so should newdata")
		if (!missing(newdata) && is(newdata, "variogramModel"))
			model = newdata
		newdata = data
		data = locations
		locations = coordinates(data)
	}
    g = gstat(formula = formula, locations = locations, model = model,
		data = data, beta = beta, nmax = nmax, nmin = nmin, 
		maxdist = maxdist, ...)
    predict.gstat(g, newdata = newdata, block = block, nsim = nsim,
		indicators = indicators, na.action = na.action)
}

idw <-
function (formula, locations = try.coordinates(data), data = sys.frame(sys.parent()), 
		newdata, nmax = Inf, nmin = 0, maxdist = Inf, block = numeric(0), 
		na.action = na.pass, idp = 2.0) {
	krige(formula, locations, data, nmax = nmax, nmin = nmin,
		maxdist = maxdist, block = block, na.action = na.action,
		set = list(idp = 2))
}
