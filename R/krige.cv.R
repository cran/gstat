"krige.cv" <- function (formula, locations, data = sys.frame(sys.parent()), 
	model = NULL, ..., beta = NULL, nmax = Inf, nmin = 0, maxdist = Inf, 
	nfold = nrow(data), verbose = FALSE) {

	gstat.cv(gstat(g = NULL, id = "var1", formula = formula, locations = 
		locations, data = data, model = model, beta = beta, nmax = Inf, 
		nmin = 0, maxdist = maxdist, ...), nfold = nfold, verbose = verbose)
}
		
#{
#	if (inherits(data, "gstatVariogram"))
#		model = data
#	if (inherits(data, "Spatial"))
#		locations = data
#	if (has.coordinates(locations)) {
#		data = locations
#		locations = coordinates(data)
#		nc = 2 + dim(locations)[2]
#		spdf = TRUE
#	} else {
#		nc = 2 + length(attr(terms(locations), "term.labels"))
#		spdf = FALSE
#	}
#	ret = data.frame(matrix(NA, nrow(data), nc))
#	if (nfold < nrow(data))
#		fold = sample(nfold, nrow(data), replace = TRUE)
#	else
#		fold = 1:nrow(data)
#	for (i in sort(unique(fold))) {
#		sel = which(fold == i)
#    	g = gstat(formula = formula, locations = locations[sel, ], model = model, 
#			data = data[-sel, ], beta = beta, nmax = nmax, nmin = nmin,
#			maxdist = maxdist, ...)
#    	x = predict.gstat(g, newdata = data[sel, ])
#    	ret[sel, ] = x
#		if (verbose)
#			print(paste("fold", i))
#	}
#	names(ret) = names(x)
#	observed = gstat.formula(formula, locations, data)$y
#	residual = observed - ret["var1.pred"]
#	zscore = residual / sqrt(ret["var1.var"])
#	ret = data.frame(ret, observed = observed, residual = residual, 
#		zscore = zscore, fold = fold)
#	names(ret) = c(names(x), "observed", "residual", "zscore", "fold")
#	if (spdf)
#		coordinates(ret) = dimnames(locations)[[2]]
#	ret
#}
