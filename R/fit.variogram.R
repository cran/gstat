"fit.variogram" <-
function (object, model, fit.sills = TRUE, fit.ranges = TRUE, 
    fit.method = 7, print.SSE = FALSE, debug.level = 1) 
{
    if (missing(object)) 
        stop("nothing to fit to")
	if (!inherits(object, "gstatVariogram"))
		stop("object should be of class variogram")
	if (length(unique(object$id)) > 1)
		stop("to use fit.variogram, variogram object should be univariable")
    if (missing(model)) 
        stop("no model to fit")
    if (!inherits(model, "variogramModel"))
        stop("model should be of class variogramModel (use vgm)")
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
	direct = attr(object, "direct")
	if (!is.null(direct)) {
		id = unique(object$id)
		if (direct[direct$id == id, "is.direct"] && any(model$psill < 0))
			stop("partial sill fitted to direct variogram is negative")
	}
    if (print.SSE) 
        print(paste("SSErr: ", ret[[4]]))
    model
}
