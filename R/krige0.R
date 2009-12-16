extractFormula = function(formula, data, newdata) {
	# extract y and X from data:
    m = model.frame(terms(formula), as(data, "data.frame"))
    y = model.extract(m, "response")
    if (length(y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)

	# extract x0 from newdata:
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
	D = 1.0 / (spDists(s0, s) ** 2)
	sumD = apply(D, 1, sum)
	D %*% y / sumD
}

CHsolve = function(A, b) {
	# solves A x = b for x if A is PD symmetric
	A = chol(A)
	backsolve(A, forwardsolve(A, b, upper.tri = TRUE, transpose = TRUE))
}

krige0 <- function(formula, data, newdata, model, beta, y) {

	stopifnot(identical(proj4string(data), proj4string(newdata)))
	lst = extractFormula(formula, data, newdata)
	X = lst$X
	x0 = lst$x0
	if (missing(y))
		y = lst$y
	ll = (!is.na(is.projected(data)) && !is.projected(data))

	s = coordinates(data)
	s0 = coordinates(newdata)
	if (is(model, "variogramModel")) {
		require(gstat)
		V = matrix(variogramLine(model, dist_vector = spDists(s, s, ll),
			covariance=TRUE)$gamma, nrow(s), nrow(s))
		v0 = matrix(variogramLine(model, dist_vector = spDists(s, s0, ll),
			covariance=TRUE)$gamma, nrow(s), nrow(s0))
	} else {
		V = model(data, data)
		v0 = model(data, newdata)
	}
	if (!missing(beta)) # sk:
		wts = CHsolve(V, v0)
	else { # ok/uk -- need to estimate beta:
		wts = CHsolve(V, cbind(v0, X))
		ViX = wts[,-(1:nrow(s0))]
		wts = wts[,1:nrow(s0)]
		beta = solve(t(X) %*% ViX, t(ViX) %*% y)
	}
	x0 %*% beta + t(wts) %*% (y - X %*% beta)
}
