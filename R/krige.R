"krige" <-
function (formula, locations, data, newdata, model = NULL, beta = NULL, 
	nmax = Inf, maxdist = Inf, block = numeric(0), nsim = 0, 
	indicators = FALSE, ...) 
{
    g = gstat(formula = formula, locations = locations, model = model, 
        data = data, beta = beta, nmax = nmax, maxdist = maxdist, ...)
    predict.gstat(g, newdata = newdata, block = block, nsim = nsim,
    	indicators = indicators)
}
