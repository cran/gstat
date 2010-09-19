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

krige0 <- function(formula, data, newdata, model, beta, y, ..., 
		computeVar = FALSE, fullCovariance = FALSE) {

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
		V = variogramLine(model, dist_vector = spDists(s, s, ll),
			covariance=TRUE)
		v0 = variogramLine(model, dist_vector = spDists(s, s0, ll),
			covariance=TRUE)
		c0 = variogramLine(model, dist_vector = c(0), covariance=TRUE)$gamma
	} else {
		V = model(data, data)
		v0 = model(data, newdata)
		d0 = data[1, 1, drop=FALSE]
		c0 = as.numeric(model(d0,d0))
	}
	if (!missing(beta)) { # sk:
		skwts = CHsolve(V, v0)
		if (computeVar)
			var = c0 - t(v0) %*% skwts
	} else { # ok/uk -- need to estimate beta:
		skwts = CHsolve(V, cbind(v0, X))
		ViX = skwts[,-(1:nrow(s0))]
		skwts = skwts[,1:nrow(s0)]
		beta = solve(t(X) %*% ViX, t(ViX) %*% y)
		if (computeVar) {
			# (x0-X'C-1 c0)'(X'C-1X)-1 (x0-X'C-1 c0) -- precompute term 1+3:
			Q = t(x0) - t(ViX) %*% v0
			var = c0 - t(v0) %*% skwts + t(Q) %*% CHsolve(t(X) %*% ViX, Q)
		}
	}
	pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
	if (computeVar) {
		if (!fullCovariance)
			var = diag(var)
		list(pred = pred, var = var)
	} else
		pred
}
