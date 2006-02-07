# $Id: gstat.formula.q,v 1.7 2006-02-10 19:01:07 edzer Exp $

"gstat.formula" <-
function (formula, data)
{
    m = model.frame(terms(formula), as(data, "data.frame"))
    Y = model.extract(m, response)
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")

	if (gridded(data))
		grid = gridparameters(data)
	else
		grid = numeric(0)

    list(y = Y, locations = coordinates(data), X = X, call = call,
        has.intercept = has.intercept, grid = as.double(unlist(grid)))
}
