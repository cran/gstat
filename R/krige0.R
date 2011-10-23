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
	A = chol(A, LINPACK=TRUE)
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
		V = model(data)
		v0 = model(data, newdata)
		if (computeVar) {
			if (is(newdata, "SpatialPolygons") || 
					is(newdata, "SpatialPolygons"))
				stop("unequal target support not implemented")
			else
				c0 = as.numeric(model(newdata[1, drop=FALSE]))
		}
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
			# here done the HARD, i.e. m x m way; first compute
			# (x0-X'C-1 c0)'(X'C-1X)-1 (x0-X'C-1 c0) 
			# -- precompute term 1+3:
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

krigeST <- function(formula, data, newdata, modelList, y, ..., 
		computeVar = FALSE, fullCovariance = FALSE) {

	if (is(data, "ST") && is(newdata, "ST")) {
		stopifnot(identical(proj4string(data@sp), proj4string(newdata@sp)))
		stopifnot(is(data, "STFDF"))
	}
	lst = extractFormula(formula, data, newdata)
	X = lst$X
	x0 = lst$x0
	if (missing(y))
		y = lst$y

	V = covfn.ST(data, model = modelList)
	stopifnot(!is.null(V$T) && !is.null(V$S))
	v0 = covfn.ST(data, newdata, model = modelList)
	stopifnot(!is.null(v0$T) && !is.null(v0$S))
	d0 = data[1, 1, drop=FALSE]
	c0 = as.numeric(covfn.ST(d0, d0, separate = FALSE, model = modelList))
	skwts = STsolve(V, v0, X)
	#ViX = skwts[,-(1:ncol(v0))]
	#skwts = skwts[,1:ncol(v0)]
	npts = prod(dim(newdata)[1:2])
	ViX = skwts[,-(1:npts)]
	skwts = skwts[,1:npts]
	beta = solve(t(X) %*% ViX, t(ViX) %*% y)
	pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
	if (computeVar) {
		# get (x0-X'C-1 c0)'(X'C-1X)-1 (x0-X'C-1 c0) -- precompute term 1+3:
		v0 = v0$Tm %x% v0$Sm
		Q = t(x0) - t(ViX) %*% v0
		var = c0 - t(v0) %*% skwts + t(Q) %*% CHsolve(t(X) %*% ViX, Q)
		if (!fullCovariance)
			var = diag(var)
		list(pred = pred, var = var)
	} else
		pred
}

STsolve = function(A, b, X) {
# V = A$T %x% A$S -- a separable covariance; solve A x = b for x
# kronecker: T %x% S vec(L) = vec(c0) <-->  S L T = c0
# solve for L: use Y = L T 
# S Y = c0 -> Y = solve(S, c0)
# L T = Y -> Tt Lt = Yt -> Lt = solve(Tt, Yt)
	Tm = chol(A$Tm, LINPACK=TRUE)
	Sm = chol(A$Sm, LINPACK=TRUE)
	STbacksolve = function(Tm, Cm, Sm) { 
		MyChSolve = function(A, b)
			backsolve(A, forwardsolve(A, b, upper.tri = TRUE, transpose = TRUE))
		# Y = MyChSolve(Sm, Cm)
		# L = MyChSolve(Tm, t(Y))
		# as.vector(t(L))
		as.vector(t(MyChSolve(Tm, t(MyChSolve(Sm, Cm)))))
	}
	# b comes separated:
	ret1 = apply(b$T, 2, function(x1) 
		apply(b$S, 2, function(x2)
			STbacksolve(Tm, matrix(x1 %x% x2, nrow(Sm), nrow(Tm)), Sm)))
	d = dim(ret1)
	dim(ret1) = c(d[1] / ncol(b$S), d[2] * ncol(b$S))
	# X comes full:
	ret2 = apply(X, 2, function(x) 
			STbacksolve(Tm, matrix(x, nrow(Sm), nrow(Tm)), Sm))
	#print(dim(ret1))
	#print(dim(ret2))
	cbind(ret1, ret2)
}

covfn.ST = function(x, y = x, separate = TRUE, model) {
	if (is(model$space, "variogramModel")) 
		Sm = variogramLine(model$space, covariance = TRUE, dist_vector = 
    		spDists(coordinates(x@sp), coordinates(y@sp)))
	else
		Sm = model$space(x, y)
	if (is(model$time, "variogramModel")) 
    	Tm = variogramLine(model$time, covariance = TRUE, dist_vector = 
    		abs(outer(as.numeric(index(x@time)), as.numeric(index(y@time)), "-")))
	else
		Tm = model$time(x, y)
    if (separate)
        list(Sm = Sm, Tm = Tm)
    else
        Tm %x% Sm # kronecker product
}
