"gstat.formula.predict" <-
function (formula, locations, newdata, na.action) 
{
	if (has.coordinates(newdata)) {
		locs = coordinates(newdata)
		newdata = as.data.frame(newdata)
	} else {
		# resolve locations:
		terms.l = terms(locations)
    	attr(terms.l, "intercept") = 0
    	mf.locs = model.frame(terms.l, newdata, na.action = na.action)
    	locs = model.matrix(terms.l, mf.locs)
	}

	# resolve formula:
	terms.f = delete.response(terms(formula))
    mf.f = model.frame(terms.f, newdata, na.action = na.action)
    X = model.matrix(terms.f, mf.f)

	if (inherits(locations, "formula") && NROW(locs) != NROW(X)) { 
		# NA's were filtered in one, but not the other:
		mf.locs = model.frame(terms.l, newdata, na.action = na.pass)
    	mf.f =    model.frame(terms.f, newdata, na.action = na.pass)
		valid.pattern = !(apply(cbind(mf.f, mf.locs), 1,
			             function(x) any(is.na(x))))
		X    = model.matrix(terms.f, mf.f   [valid.pattern, , drop = FALSE])
		locs = model.matrix(terms.l, mf.locs[valid.pattern, , drop = FALSE])
		if (NROW(locs) != NROW(X))
			stop("NROW(locs) != NROW(X): this should not occur")
	}
    list(locations = as.matrix(locs), X = as.matrix(X))
}
