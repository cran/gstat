"krige.cv" <-
function (formula, locations, data, model = NULL, beta = NULL, nmax = Inf, 
	maxdist = Inf, nfold = nrow(data), verbose = TRUE, ...)
{
	nc = 2 + length(attr(terms(~x+y),"term.labels"))
	ret = data.frame(matrix(NA, nrow(data), nc))
	if (nfold < nrow(data))
		fold = sample(nfold, nrow(data), replace = TRUE)
	else
		fold = 1:nrow(data)
	for (i in sort(unique(fold))) {
		sel = (1:nrow(data))[fold == i]
    	g = gstat(formula = formula, locations = locations, model = model, 
			data = data[-sel,], beta = beta, nmax = nmax, 
			maxdist = maxdist, ...)
    	x = predict.gstat(g, newdata = data[sel,])
    	ret[sel,] = x
		if(verbose)
			print(paste("fold", i))
	}
	names(ret) = names(x)
	observed = gstat.formula(formula, locations, data)$y
	residual = observed - ret["var1.pred"]
	zscore = residual / sqrt(ret["var1.var"])
	ret = data.frame(ret, observed = observed, residual = residual, 
		zscore = zscore, fold = fold)
	names(ret) = c(names(x), "observed", "residual", "zscore", "fold")
	ret
}
