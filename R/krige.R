"krige" <-
function (formula, locations = try.coordinates(data), data = sys.frame(sys.parent()), 
	newdata, model = NULL, beta = NULL, nmax = Inf, nmin = 0, 
	maxdist = Inf, block = numeric(0), nsim = 0, indicators = FALSE, 
	na.action = na.pass, ...)
{
	if (has.coordinates(locations)) { # shift arguments:
		if (has.coordinates(data)) # another shift:
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
