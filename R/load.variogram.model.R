# $Id: load.variogram.model.q,v 1.7 2008-11-12 10:04:22 edzer Exp $

"load.variogram.model" <- function(model, ids = c(0, 0)) {
	if (missing(model))
		stop("model is missing");
	if (!inherits(model, "variogramModel"))
		stop("model should be of mode variogramModel (use function vgm)")
	if (any(model$range < 0.0)) {
		print(model)
		stop("variogram range can never be negative")
	}
	anis = c(model$ang1, model$ang2, model$ang3, model$anis1, model$anis2)
	if (is.null(attr(model, "table")))
		covtable = numeric(0)
	else  {
		covtable = attr(model, "table")
		if (dim(model)[1] > 1 || model$model != "Tab")
			stop("table can only have one single model")
	}
	.Call(gstat_load_variogram, 
		as.integer(ids),
		as.character(model$model),
		as.numeric(model$psill),
		as.numeric(model$range),
		as.numeric(model$kappa),
		as.numeric(anis),
		covtable)
}
