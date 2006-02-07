# $Id: fit.variogram.q,v 1.6 2006-02-10 19:01:07 edzer Exp $

"fit.variogram" <-
function (object, model, fit.sills = TRUE, fit.ranges = TRUE, 
    fit.method = 7, debug.level = 1, warn.if.neg = FALSE) 
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
    .Call("gstat_init", as.integer(debug.level))
    .Call("gstat_load_ev", object$np, object$dist, object$gamma)
    load.variogram.model(model)
    ret = .Call("gstat_fit_variogram", as.integer(fit.method), 
        as.integer(fit.sills), as.integer(fit.ranges))
    .Call("gstat_exit", 0)
    model$psill = ret[[1]]
    model$range = ret[[2]]
	attr(model, "singular") = as.logical(ret[[3]]);
	attr(model, "SSErr") = ret[[4]]
	direct = attr(object, "direct")
	if (!is.null(direct)) {
		id = unique(object$id)
		if (direct[direct$id == id, "is.direct"] && any(model$psill < 0)) {
			if (warn.if.neg)
				warning("partial sill or nugget fixed at zero value")
			fit.sills = model$psill > 0
			model$psill[model$psill < 0] = 0.0
			return(fit.variogram(object, model, fit.sills = fit.sills, fit.ranges =
				fit.ranges, fit.method = fit.method, debug.level = debug.level,
				warn.if.neg = warn.if.neg))
		}
	}
    model
}
