"fit.variogram" <-
function (object, model, fit.sills = TRUE, fit.ranges = TRUE, 
    fit.method = 7, print.SSE = FALSE, debug.level = 1) 
{
    if (missing(object)) 
        stop("nothing to fit to")
	if (class(object) != "variogram")
		stop("object should be of class variogram")
    if (missing(model)) 
        stop("no model to fit")
    if (class(model) != "variogram.model") 
        stop("model should be of class variogram.model (use vgm)")
    if (length(fit.sills) < length(model$model)) 
        fit.sills = rep(fit.sills, length(model$model))
    if (length(fit.ranges) < length(model$model)) 
        fit.ranges = rep(fit.ranges, length(model$model))
    fit.ranges = fit.ranges & (model$model != "Nug")
    .Call("gstat_init", as.integer(debug.level))
    .Call("gstat_load_ev", object$np, object$dist, object$gamma)
    load.variogram.model(model)
    ret = .Call("gstat_fit_variogram", as.integer(fit.method), 
        as.integer(fit.sills), as.integer(fit.ranges))
    .Call("gstat_exit", 0)
    model$psill = ret[[1]]
    model$range = ret[[2]]
    if (print.SSE) 
        print(paste("SSErr: ", ret[[3]]))
    model
}
