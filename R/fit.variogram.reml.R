"fit.variogram.reml" <-
function (formula, locations, data, model, debug.level = 1, set)
{
    if (missing(formula)) 
        stop("nothing to fit to")
    if (class(formula) != "formula") 
        stop("formula should be of class formula")
    if (missing(model)) 
        stop("no model to fit")
    if (!inherits(model, "variogram.model"))
        stop("model should be of class variogram.model (use vgm)")
    fit.sills = rep(TRUE, length(model$model))
    fit.ranges = rep(FALSE, length(model$model))
    .Call("gstat_init", as.integer(debug.level)
    	, PACKAGE = "gstat"
	)
    ret = gstat.formula(formula, locations, data)
    ret$y <- residuals(lm(formula, data))
    .Call("gstat_new_data", as.vector(ret$y), as.vector(ret$locations),
         as.vector(ret$X), as.integer(1), as.vector(numeric(0)),
         as.integer(-1), as.numeric(-1), as.integer(1)
	 , PACKAGE = "gstat"
	 )
    load.variogram.model(model)
    if (!missing(set))
    	gstat.load.set(set)
    ret = .Call("gstat_fit_variogram", as.integer(5), 
        as.integer(fit.sills), as.integer(fit.ranges)
	, PACKAGE = "gstat"
	)
    .Call("gstat_exit", 0, PACKAGE = "gstat")
    model$psill = ret[[1]]
    model$range = ret[[2]]
    model
}
