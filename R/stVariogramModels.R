# constructiong spatio-temporal variogram models
vgmST <- function(stModel, ..., space, time, joint, sill, nugget, stAni) {
	stopifnot(is.character(stModel) && length(stModel)==1)
	vgmModel <- switch(stModel,
		separable = list(space = space, time = time, sill = sill),
		productSum = list(space = space, time = time, sill = sill, nugget = nugget),
		sumMetric = list(space = space, time = time, joint = joint, stAni = stAni),
		simpleSumMetric = list(space = space, time = time, 
			joint = joint, nugget = nugget, stAni = stAni),
		metric = list(joint = joint, stAni = stAni),
		stop(paste("model", stModel, "unknown")))
  vgmModel$stModel <- stModel
  class(vgmModel) <- c("StVariogramModel","list")
  vgmModel
}

# calculating spatio-temporal variogram surfaces
variogramSurface <- function(model, dist_grid, covariance=FALSE, ...) {
  if (!inherits(model, "StVariogramModel"))
    warning("\"model\" should be of class \"StVariogramModel\"; no further checks for a proper will made.")
  
  switch(model$stModel,
         separable=vgmSeparable(model, dist_grid, covariance, ...),
         productSum=vgmProdSum(model, dist_grid, covariance, ...),
         sumMetric=vgmSumMetric(model, dist_grid, covariance, ...),
         simpleSumMetric=vgmSimpleSumMetric(model, dist_grid, covariance, ...),
         metric=vgmMetric(model, dist_grid, covariance, ...),
         stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

# separable model: C_s * C_t
vgmSeparable <- function(model, dist_grid, covariance, separate=TRUE) {
  if(covariance)
    return(covSep(x=dist_grid$spacelag, y=dist_grid$timelag, model, separate=separate)) # to be integrated
  
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag)[,2]

  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, 
             model=model$sill*(vs+vt-vs*vt))
}

# productSum model: C_s*C_t + C_s + C_t
vgmProdSum <- function(model, dist_grid, covariance) {
  if(covariance)
    return(covProdSum(x=dist_grid$spacelag, y=dist_grid$timelag, model)) # to be integrated
  
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag)[,2]
  vt = variogramLine(model$time, dist_vector=dist_grid$timelag)[,2]

  k <- (sum(model$space$psill)+sum(model$time$psill)-model$sill)/(sum(model$space$psill)*sum(model$time$psill))
  
  if (k <= 0 | k > 1/max(model$space$psill[2],model$time$psill[2])) 
    k <- 10^6*abs(k) # distorting the model to let optim "hopefully" find suitable parameters
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=as.vector(vs+vt-k*vs*vt+model$nugget))
}

# sumMetric model: C_s + C_t + C_st (Gerard Heuvelink)
vgmSumMetric <- function(model, dist_grid, covariance) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag, covariance=covariance)[,2]
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag,  covariance=covariance)[,2]
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vst = variogramLine(model$joint,dist_vector=h, covariance=covariance)[,2]
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=(vs + vt + vst))
}

# simplified sumMetric model
vgmSimpleSumMetric <- function(model, dist_grid, covariance) {
  vs = variogramLine(model$space, dist_vector=dist_grid$spacelag, covariance=covariance)[,2]
  vt = variogramLine(model$time,  dist_vector=dist_grid$timelag,  covariance=covariance)[,2]
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vm = variogramLine(model$joint,dist_vector=h, covariance=covariance)[,2]
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=(vs + vt + vm + model$nugget))
}

vgmMetric <- function(model, dist_grid, covariance) {
  h = sqrt(dist_grid$spacelag^2 + (model$stAni * as.numeric(dist_grid$timelag))^2)
  vm = variogramLine(model$joint, dist_vector=h, covariance=covariance)[,2]
  data.frame(spacelag=dist_grid$spacelag, timelag=dist_grid$timelag, model=vm)
}

fit.StVariogram <- function(object, model, ..., wles=FALSE) {
  object <- na.omit(object)
  if (!inherits(object, "StVariogram"))
    stop("\"object\" must be of class \"StVariogram\"")
  if (!inherits(model, "StVariogramModel"))
    stop("\"model\" must be of class \"StVariogramModel\".")
  
  fitFun = function(par, trace = FALSE, ...) {
    resid = object$gamma - variogramSurface(insertPar(par,model), object[,c("spacelag","timelag")])$model
    if (trace)
      print(c(par, MSE = mean(resid^2)))
    if(wles)
      resid <- resid * object$np/sum(object$np)
    sqrt(mean(resid^2))
  }
  pars.fit <- optim(extractPar(model), fitFun, ...)
  
  ret <- insertPar(pars.fit$par, model)
  attr(ret,"optim.output") <- pars.fit
  ret
}

insertPar <- function(par, model) {
  switch(model$stModel,
         separable=insertParSeparable(par, model),
         productSum=insertParProdSum(par, model),
         sumMetric=insertParSumMetric(par, model),
         simpleSumMetric=insertParSimpleSumMetric(par,model),
         metric=insertParMetric(par,model),
         stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

extractPar <- function(model) {
  switch(model$stModel,
         separable=c(range.s=model$space$range[2], nugget.s=model$space$psill[1],
                     range.t=model$time$range[2],  nugget.t=model$time$psill[1],
                     sill= model$sill),
         productSum=c(sill.s = model$space$psill[2], range.s = model$space$range[2],
                      sill.t = model$time$psill[2],  range.t = model$time$range[2], 
                      sill=model$sill, nugget=model$nugget),
         sumMetric=c(sill.s = model$space$psill[2], range.s = model$space$range[2], nugget.s = model$space$psill[1], 
                     sill.t = model$time$psill[2], range.t = model$time$range[2], nugget.t = model$time$psill[1],
                     sill.st = model$joint$psill[2], range.st = model$joint$range[2], nugget.st = model$joint$psill[1],
                     anis = model$stAni),
# simplified sumMetric model
         simpleSumMetric=c(sill.s = model$space$psill[2], range.s = model$space$range[2], 
                           sill.t = model$time$psill[2], range.t = model$time$range[2],
                           sill.st = model$joint$psill[2], range.st = model$joint$range[2], 
                           nugget = model$nugget, anis = model$stAni),
         metric=c(sill = model$joint$psill[2], range = model$joint$range[2], nugget = model$joint$psill[1],
                  anis = model$stAni),
         stop("Only \"separable\", \"productSum\", \"sumMetric\", \"simpleSumMetric\" and \"metric\" are implemented."))
}

insertParSeparable <- function(par, model) {
  vgmST("separable",
        space=vgm(1-par[2],as.character(model$space$model[2]),par[1],par[2]),
        time= vgm(1-par[4],as.character(model$time$model[2]),par[3],par[4]),
        sill=par[5])
}

insertParProdSum <- function(par, model) {
  vgmST("productSum",
        space=vgm(par[1],as.character(model$space$model[2]),par[2],0),
        time= vgm(par[3],as.character(model$time$model[2]),par[4],0),
        sill=par[5], nugget=par[6])
}

insertParSumMetric <- function(par, model) {
  vgmST("sumMetric",
        space=vgm(par[1],as.character(model$space$model[2]),par[2],par[3]),
        time= vgm(par[4],as.character(model$time$model[2]),par[5],par[6]),
        joint=vgm(par[7],as.character(model$joint$model[2]),par[8],par[9]),
        stAni=par[10])
}

# simplified sumMetric model
insertParSimpleSumMetric <- function(par, model) {
  vgmST("simpleSumMetric",
        space=vgm(par[1],as.character(model$space$model[2]),par[2],0),
        time= vgm(par[3],as.character(model$time$model[2]),par[4],0),
        joint=vgm(par[5],as.character(model$joint$model[2]),par[6],0),
        nugget=par[7], stAni=par[8])
}

insertParMetric <- function(par, model) {
  vgmST("metric",
        joint=vgm(par[1], as.character(model$joint$model[2]), par[2], par[3]),
        stAni=par[4])
}
