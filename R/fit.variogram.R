"fit.variogram" <-
function (object, model, fit.sills = TRUE, fit.ranges = TRUE, 
    fit.method = 7, print.SSE = FALSE, debug.level = 1) 
{
    if (missing(object)) 
        stop("nothing to fit to")
	if (!inherits(object, "variogram"))
		stop("object should be of class variogram")
    if (missing(model)) 
        stop("no model to fit")
    if (!inherits(model, "variogram.model"))
        stop("model should be of class variogram.model (use vgm)")
    if (fit.method == 5)
    	stop("use function fit.variogram.reml() to use REML")
    if (length(fit.sills) < length(model$model)) 
        fit.sills = rep(fit.sills, length(model$model))
    if (length(fit.ranges) < length(model$model)) 
        fit.ranges = rep(fit.ranges, length(model$model))
    fit.ranges = fit.ranges & (model$model != "Nug")
    .Call("gstat_init", as.integer(debug.level), PACKAGE = "gstat")
    .Call("gstat_load_ev", object$np, object$dist, object$gamma, 
		PACKAGE = "gstat")
    load.variogram.model(model)
    ret = .Call("gstat_fit_variogram", as.integer(fit.method), 
        as.integer(fit.sills), as.integer(fit.ranges), PACKAGE = "gstat")
    .Call("gstat_exit", 0, PACKAGE = "gstat")
    model$psill = ret[[1]]
    model$range = ret[[2]]
	attr(model, "singular") = as.logical(ret[[3]]);
    if (print.SSE) 
        print(paste("SSErr: ", ret[[4]]))
    model
}
