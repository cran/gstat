"gstat.formula" <-
function (formula, locations, data)
{
	grid = numeric(0)
    if (has.coordinates(locations)) { # ignore data: locations arg has coords+data
        data = as.data.frame(locations)
		grid = try.gridparameters(locations)
        locations = coordinates(locations)
    } else if (has.coordinates(data)) { # ignore locations: data arg name was used
        locations = coordinates(data)
		grid = try.gridparameters(data)
        data = as.data.frame(data)
    } else { # resolve locations formula from data.frame:
        if (!inherits(locations, "formula"))
            stop("locations argument should be a formula, such as ~x+y")
		m = model.frame(terms(locations), data)
        Terms = attr(m, "terms")
        attr(Terms, "intercept") = 0
        if ((yvar = attr(Terms, "response")) > 0)
            stop("no response allowed in locations formula")
		# retrieve coord columns from model frame:
        locations = model.matrix(Terms, m)
    }
	# now extract main formula from data:
    m = model.frame(terms(formula), data)
    Y = model.extract(m, response)
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")

    list(y = Y, locations = as.matrix(locations), X = X, call = call,
        has.intercept = has.intercept, grid = as.double(unlist(grid)))
}
