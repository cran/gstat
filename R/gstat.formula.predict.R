"gstat.formula.predict" <-
function (formula, newdata, na.action) 
{
	if (is(newdata, "SpatialPolygons")) {
		locs = t(sapply(slot(newdata, "polygons"), function(x) slot(x, "labpt")))
		colnames(locs) = c("x", "y")
		if (is(newdata, "SpatialPolygonsDataFrame"))
			newdata = as.data.frame(newdata)
		else
			newdata = data.frame(a = rep(1, nrow(locs)))
	} else {
		locs = coordinates(newdata)
		newdata = as.data.frame(newdata)
	} 

	# resolve formula:
	terms.f = delete.response(terms(formula))
    mf.f = model.frame(terms.f, newdata, na.action = na.action)
    X = model.matrix(terms.f, mf.f)

	if (NROW(locs) != NROW(X)) { 
		# NA's were filtered in X, but not in coords:
    	mf.f =    model.frame(terms.f, newdata, na.action = na.pass)
		valid.pattern = !(apply(mf.f, 1, function(x) any(is.na(x))))
		X    = model.matrix(terms.f, mf.f   [valid.pattern, , drop = FALSE])
		locs = locs[valid.pattern, ]
		if (NROW(locs) != NROW(X))
			stop("NROW(locs) != NROW(X): this should not occur")
	}
    list(locations = as.matrix(locs), X = as.matrix(X))
}
