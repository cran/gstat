"krige" <-
function (formula, locations, data, newdata, model = NULL, beta = NULL, 
	nmax = Inf, block = numeric(0), nsim = 0, ...) 
{
    g = gstat(formula = formula, locations = locations, model = model, 
        data = data, beta = beta, nmax = nmax, ...)
    predict.gstat(g, newdata = newdata, block = block, nsim = nsim)
}
