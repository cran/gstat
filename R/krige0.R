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

idw0 = function(formula, data, newdata, y, idp = 2.0) {
	s = coordinates(data)
	s0 = coordinates(newdata)
	if (missing(y))
		y = extractFormula(formula, data, newdata)$y
	D = 1.0 / (spDists(s0, s) ^ idp)
	sumD = apply(D, 1, sum)
	D %*% y / sumD
}

CHsolve = function(A, b) {
	# solves A x = b for x if A is PD symmetric
	#A = chol(A, LINPACK=TRUE) -> deprecated
	A = chol(A) # but use pivot=TRUE?
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
		V = model(data, data, ...)
		v0 = model(data, newdata, ...)
		if (computeVar) {
			if (is(newdata, "SpatialLines") || 
					is(newdata, "SpatialPolygons"))
				stop("varying target support has not been implemented")
			c0 = as.numeric(model(newdata[1, drop=FALSE],
				newdata[1, drop=FALSE]))
			# ?check this: provide TWO arguments, so model(x,y) can target
			# eventually Y, instead of measurements Z=Y+e
			# with e measurement error term e
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
	stopifnot(inherits(modelList, "StVariogramModel"))
	stopifnot(inherits(data, c("STF", "STS", "STI")) & inherits(newdata, c("STF", "STS", "STI"))) 
	stopifnot(identical(proj4string(data@sp), proj4string(newdata@sp)))
  stopifnot(class(data@time) == class(newdata@time))
  
	if(is.null(attr(modelList,"temporal unit")))
	  warning("The spatio-temporal variogram model does not carry a time unit attribute: krisgeST cannot check whether the temporal distance metrics coincide.")

  separate <- length(data) > 1 & length(newdata) > 1 & inherits(data, "STF") & inherits(newdata, "STF")
  
	lst = extractFormula(formula, data, newdata)
	X = lst$X
	x0 = lst$x0
	if (missing(y))
		y = lst$y

	V = covfn.ST(data, model = modelList, separate=separate)
	v0 = covfn.ST(data, newdata, modelList)

	if (is(data,"STSDF"))
		d0 <- data[data@index[1,1],data@index[1,2],drop=F]
	else
    	d0 = data[1, 1, drop=FALSE]
	c0 = as.numeric(covfn.ST(d0, d0, modelList, separate = FALSE))
	if(modelList$stModel == "separable" & separate)
    skwts <- STsolve(V, v0, X)
  else 
    skwts <- CHsolve(V, cbind(v0,X))
	# ViX = skwts[,-(1:ncol(v0))]
	# skwts = skwts[,1:ncol(v0)]
	#npts = prod(dim(newdata)[1:2]) #-> does not work for STI
	npts = length(newdata)
	ViX = skwts[,-(1:npts)]
	skwts = skwts[,1:npts]
	beta = solve(t(X) %*% ViX, t(ViX) %*% y)
	pred = x0 %*% beta + t(skwts) %*% (y - X %*% beta)
	if (computeVar) {
		# get (x0-X'C-1 c0)'(X'C-1X)-1 (x0-X'C-1 c0) -- precompute term 1+3:
		if(is.list(v0)) # in the separable case
			v0 = v0$Tm %x% v0$Sm
		Q = t(x0) - t(ViX) %*% v0
		var = c0 - t(v0) %*% skwts + t(Q) %*% CHsolve(t(X) %*% ViX, Q)
		if (!fullCovariance)
			var = diag(var)
		res <- data.frame(pred = pred, var = var)
	} else
		res <- data.frame(pred)
  
  # wrapping the predictions in ST*DF again
	if (ncol(res) == 1)
	  names(res) = "var1.pred"
	if (ncol(res) == 2)
	  names(res) = c("var1.pred", "var1.var")
	addAttrToGeom(geometry(newdata), res)
}

STsolve = function(A, b, X) {
# V = A$T %x% A$S -- a separable covariance; solve A x = b for x
# kronecker: T %x% S vec(L) = vec(c0) <-->  S L T = c0
# solve for L: use Y = L T 
# S Y = c0 -> Y = solve(S, c0)
# L T = Y -> Tt Lt = Yt -> Lt = solve(Tt, Yt)
	#Tm = chol(A$Tm, LINPACK=TRUE)
	Tm = chol(A$Tm)
	#Sm = chol(A$Sm, LINPACK=TRUE)
	Sm = chol(A$Sm)
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

covfn.ST = function(x, y = x, model, ...) {
  switch(model$stModel,
         separable=covSeparable(x, y, model, ...),
         productSum=covProdSum(x, y, model),
         sumMetric=covSumMetric(x, y, model),
         simpleSumMetric=covSimpleSumMetric(x, y, model),
         metric=covMetric(x, y, model),
         stop(paste("Provided spatio-temporal model (",model$stModel,") is not supported.",sep="")))
}

## covariance models
####################

# separable covariance model
# covSep <- function(x, y, model, separate=TRUE) {
#   if (is(model$space, "variogramModel")) 
#     Sm = variogramLine(model$space, covariance = TRUE, dist_vector = 
#       spDists(coordinates(x@sp), coordinates(y@sp)))*model$sill
#   else
#     Sm = model$space(x, y)
#   stopifnot(!is.null(Sm))
#   if (is(model$time, "variogramModel")) 
#     Tm = variogramLine(model$time, covariance = TRUE, dist_vector = 
#       abs(outer(as.numeric(index(x@time)), as.numeric(index(y@time)), "-")))
#   else
#     Tm = model$time(x, y)
#   stopifnot(!is.null(Tm))
#   if (separate)
#     list(Sm = Sm, Tm = Tm)
#   else
#     Tm %x% Sm # kronecker product
# }

# tunit <- attr(modelList, "temporal unit")
# if(!is.null(tunit)) {
#   scaleFac <- switch(tunit, 
#                      days=24*60*60, 
#                      hours=60*60, 
#                      mins=60, 
#                      secs=1,
#                      -1)
#   if(scaleFac == -1) {
#     warning(tunit, "is an unknown time unit to krigeST. No attempt to correct the temporal scale is made.")
#     scaleFac <- 1
#   }
#   if(scaleFac > 1) {
#     cat(paste("The variogram model has a temporal unit (", tunit,
#               ") different than the temporal metric of krigeST (secs). ",
#               "An corrrection attempt has been made by the factor: ",
#               scaleFac, sep=""))
#     if(modelList$stModel %in% c("separable", "productSum", "sumMetric", "simpleSumMetric"))
#       modelList$time$range <- modelList$time$range*scaleFac
#     if(modelList$stModel %in% c("metric", "sumMetric", "simpleSumMetric"))
#       modelList$stAni <- modelList$stAni/scaleFac
#   }
# }


covSeparable <- function(x, y, model, separate) {  
  if(missing(separate))
    separate <- inherits(x, "STF") & inherits(y, "STF") & length(x) > 1 & length(y) > 1

  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)*model$sill
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
        
    if (separate)
      return(list(Sm = Sm, Tm = Tm))
    else
      return(Tm %x% Sm) # kronecker product
  } 
  
  # separate makes only sense if both of x and y inherit STF
  if(separate)
    stop("An efficient inversion by separating the covarinace model is only possible if both of \"x\" and \"y\" inherit \"STF\"")
  
  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    
    
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)*model$sill
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    return(Sm * Tm)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(coordinates(x@sp), coordinates(y@sp))
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  Sm = variogramLine(model$space, covariance = TRUE, dist_vector = sMat)*model$sill
  Tm = variogramLine(model$time, covariance = TRUE, dist_vector = tMat)
  
  return(Sm * Tm)  
}

## product-sum model, BG
covProdSum <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  
  # double check model for validity, i.e. k:
  k <- (sum(model$space$psill)+sum(model$time$psill)-model$sill)/(sum(model$space$psill)*sum(model$time$psill))
  if (k <= 0 | k > 1/max(model$space$psill[model$space$model!="Nug"], 
                         model$time$psill[model$time$model!="Nug"]))
    stop(paste("k (",k,") is non-positive or too large: no valid model!",sep=""))
  
  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds)
    vt = variogramLine(model$time, dist_vector = dt)
     
    return(model$sill-(vt %x% matrix(1,nrow(vs),ncol(vs)) + matrix(1,nrow(vt),ncol(vt)) %x% vs - k * vt %x% vs))
  } 

  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    vs = variogramLine(model$space, dist_vector = ds)
    vt = variogramLine(model$time, dist_vector = dt)
    
    return(model$sill-(vt + vs - k * vt * vs))
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(coordinates(x@sp), coordinates(y@sp))
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  vs = variogramLine(model$space, dist_vector = sMat)
  vt = variogramLine(model$time, dist_vector = tMat)
  
  return(model$sill-(vt + vs - k * vt * vs))
}

## sumMetric model
covSumMetric <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  
  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    h  = sqrt((matrix(1,nrow(dt),ncol(dt)) %x% ds)^2 
              + (model$stAni * dt %x% matrix(1,nrow(ds),ncol(ds)))^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(matrix(1,nrow(Tm),ncol(Tm)) %x% Sm + Tm %x% matrix(1,nrow(Sm),ncol(Sm)) + Mm)
  } 
  
  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    Sm = variogramLine(model$space, covariance = TRUE, dist_vector = ds)
    Tm = variogramLine(model$time, covariance = TRUE, dist_vector = dt)
    
    h  = sqrt(ds^2 + (model$stAni * dt)^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Sm + Tm + Mm)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(coordinates(x@sp), coordinates(y@sp))
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  Sm = variogramLine(model$space, covariance = TRUE, dist_vector = sMat)
  Tm = variogramLine(model$time, covariance = TRUE, dist_vector = tMat)
  
  h  = sqrt(sMat^2 + (model$stAni * tMat)^2)
  Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
  
  return(Sm + Tm + Mm)
}

## simple sumMetric model
covSimpleSumMetric <- function(x, y, model) {
  covSumMetric(x, y, vgmST("sumMetric", space=model$space, time=model$time,
                           joint=vgm(sum(model$joint$psill), 
                                     model$joint$model[model$joint$model != "Nug"],
                                     model$joint$range, model$nugget),
                           stAni=model$stAni))
}

## metric model
covMetric <- function(x, y, model) {
  stopifnot(inherits(x, c("STF", "STS", "STI")) & inherits(y, c("STF", "STS", "STI")))
  
  # the STF case
  if (inherits(x, "STF") & inherits(y, "STF")) {
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    h  = sqrt((matrix(1,nrow(dt),ncol(dt)) %x% ds)^2
              + (model$stAni * dt %x% matrix(1,nrow(ds),ncol(ds)))^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Mm)
  } 
  
  # the STI case
  if(inherits(x, "STI") | inherits(y, "STI")) {
    # make sure that now both are of type STI
    x <- as(x, "STI")
    y <- as(y, "STI")
    
    # calculate all spatial and temporal distances
    ds = spDists(coordinates(x@sp), coordinates(y@sp))
    dt = abs(outer(index(x@time), index(y@time), "-"))
    if(!is.null(attr(model,"temporal unit")))
      units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
    dt <- as(dt, "matrix")
    
    # compose the cov-matrix
    h  = sqrt(ds^2 + (model$stAni * dt)^2)
    Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
    
    return(Mm)
  }
  
  # the remaining cases, none of x and y is STI nor are both STF
  # make sure both are of type STS
  x <- as(x, "STS")
  y <- as(y, "STS")
  
  # calculate all spatial and temporal distances
  ds = spDists(coordinates(x@sp), coordinates(y@sp))
  dt = abs(outer(index(x@time), index(y@time), "-"))
  if(!is.null(attr(model,"temporal unit")))
    units(dt) <- attr(model, "temporal unit") # ensure the same temporal metric as in the variogram definition
  dt <- as(dt, "matrix")
  
  # re-arrange the spatial and temporal distances
  sMat <- matrix(NA, nrow(x@index), nrow(y@index))
  tMat <- matrix(NA, nrow(x@index), nrow(y@index))
  for(r in 1:nrow(x@index)) {
    sMat[r,] <- ds[x@index[r,1], y@index[,1]]
    tMat[r,] <- dt[x@index[r,2], y@index[,2]]
  }
  
  # compose the cov-matrix
  h  = sqrt(sMat^2 + (model$stAni * tMat)^2)
  Mm = variogramLine(model$joint, covariance = TRUE, dist_vector = h)
  
  return(Mm)
}
