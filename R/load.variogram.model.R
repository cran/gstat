"load.variogram.model" <- function(model, ids = c(0, 0)) 
{
	if (missing(model))
		stop("model is missing");
	if (!inherits(model, "variogramModel"))
		stop("model should be of mode variogramModel (use function vgm)")
	anis = c(model$ang1, model$ang2, model$ang3, model$anis1, model$anis2)
	.C("Cgstat_load_variogram", 
		as.integer(ids),
		as.integer(length(model$model)),
		as.character(model$model),
		as.numeric(model$psill),
		as.numeric(model$range),
		as.numeric(model$kappa),
		as.numeric(anis)
		, PACKAGE = "gstat"
		)
}
