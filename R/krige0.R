dst2 = function(s, s0) { # assumes 2D:
	d = outer(s[,1], s0[,1], "-") ^ 2
	if (ncol(s) > 1) {
		for (i in 2:ncol(s))
			d = d + outer(s[,i], s0[,i], "-") ^ 2
	}
	d
}

extractFormula = function(formula, data, newdata) {
	# extract y and X from data
    m = model.frame(terms(formula), as(data, "data.frame"))
    y = model.extract(m, "response")
    if (length(y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)

	# extract x0 from newdata
    terms.f = delete.response(terms(formula))
    mf.f = model.frame(terms.f, newdata) #, na.action = na.action)
    x0 = model.matrix(terms.f, mf.f)
	list(y = y, X = X, x0 = x0)
}

idw0 = function(formula, data, newdata, y) {
	s = coordinates(data)
	s0 = coordinates(newdata)
	if (missing(y))
		y = extractFormula(formula, data, newdata)$y
	D = 1.0 / dst2(s0, s)
	sumD = apply(D, 1, sum)
	D %*% y / sumD
}

CHsolve = function(A, b) {
	# use CHsolve, as it assumes A to be symmetric instead of full
	ch = chol(A)
	backsolve(ch, forwardsolve(ch, b, upper = TRUE, trans=TRUE))
}

krige0 <- function(formula, data, newdata, model, beta, y) {

	lst = extractFormula(formula, data, newdata)
	X = lst$X
	x0 = lst$x0
	if (missing(y))
		y = lst$y

	s = coordinates(data)
	s0 = coordinates(newdata)
	if (is(model, "variogramModel")) {
		require(gstat)
		V = variogramLine(model, dist_vector = sqrt(dst2(s, s)),
			covariance=TRUE)
		V = matrix(V$gamma, nrow(s), nrow(s))
		v0 = variogramLine(model, dist_vector = sqrt(dst2(s, s0)),
			covariance=TRUE)
		v0 = matrix(v0$gamma, nrow(s), nrow(s0))
	} else {
		V = model(sqrt(dst2(s, s)))
		v0 = model(sqrt(dst2(s, s0)))
	}
	if (!missing(beta)) # sk:
		ret = x0 %*% beta + t(CHsolve(V, v0)) %*% (y - X %*% beta)
	else { # ok/uk:
		wts = CHsolve(V, cbind(v0, X))
		ViX = wts[,-(1:nrow(s0))]
		v0tVi = wts[,1:nrow(s0)]
		beta = solve(t(X) %*% ViX, t(ViX) %*% y)
		ret = x0 %*% beta + t(v0tVi) %*% (y - X %*% beta)
	}
	ret
}
