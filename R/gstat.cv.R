"gstat.cv" <-
function (object, nfold = nrow(object$data[[1]]$data), remove.all = FALSE, 
	verbose = FALSE, all.residuals = FALSE, ...) 
{
	if (!inherits(object, "gstat")) 
		stop("first argument should be of class gstat")
	var1 = object$data[[1]]
	data = var1$data
	formula = var1$formula
	locations = var1$locations
	if (all.residuals)
		nc = length(object$data)
	else
		nc = 2 + length(attr(terms(locations), "term.labels"))
	ret = data.frame(matrix(NA, nrow(data), nc))
	if (missing(nfold)) 
		nfold = nrow(data)
	if (nfold < nrow(data)) 
		fold = sample(nfold, nrow(data), replace = TRUE)
	else fold = 1:nrow(data)

	if (all.residuals || (remove.all && length(object$data) > 1)) {
		all.data = list()
		for (v in 1:length(object$data))
			all.data[[v]] = object$data[[v]]$data
	}

	for (i in sort(unique(fold))) {
		sel = which(fold == i)
		object$data[[1]]$data = data[-sel, ]
		if (remove.all && length(object$data) > 1) {
			for (v in 2:length(object$data)) {
				varv = object$data[[v]]
				varv$data = all.data[[v]]
				atv = gstat.formula(varv$formula, varv$locations, 
				  varv$data)$locations
				at1 = gstat.formula(formula, locations, data[sel, 
				  ])$locations
				all = rbind(atv, at1)
				if (length(attr(terms(~x+y), "term.labels")) == 2) # 2-D
					zd = zerodist(all[, 1], all[, 2])
				else 
					zd = zerodist(all[, 1], all[, 2], all[, 3])
				skip = zd[, 1]
				object$data[[v]]$data = varv$data[-skip, ]
			}
		}
		x = predict.gstat(object, newdata = data[sel, ], ...)
		if (verbose) 
			print(paste("fold", i))
		if (all.residuals) {
			for (i in 1:length(object$data)) {
				var.i = object$data[[i]]
				data.i = all.data[[i]]
				formula.i = var.i$formula
				locations.i = var.i$locations
				observed = gstat.formula(formula.i, locations.i, data.i)$y[sel]
				pred.name = paste(names(object$data)[i], "pred", sep = ".")
				residual = as.numeric(observed - x[pred.name])
				ret[sel, i] = residual
			}
		} else 
			ret[sel, 1:nc] = x[, 1:nc]
	}

	if (! all.residuals) {
		names(ret) = names(x)[1:nc]
		observed = gstat.formula(formula, locations, data)$y
		pred.name = paste(names(object$data)[1], "pred", sep = ".")
		residual = observed - ret[pred.name]
		var.name = paste(names(object$data)[1], "var", sep = ".")
		zscore = residual/sqrt(ret[var.name])
		ret = data.frame(ret, observed = observed, residual = residual, 
			zscore = zscore, fold = fold)
		names(ret) = c(names(x)[1:nc], "observed", "residual", "zscore", "fold")
	} else {
		names(ret) = names(object$data)
	}
	ret
}
