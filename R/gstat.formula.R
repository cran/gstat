"gstat.formula" <-
function (formula, locations, data)
{
    if (has.coordinates(locations)) {
        data = as.data.frame(locations)
        locations = coordinates(locations)
    } else if (has.coordinates(data)) {
        locations = coordinates(data)
        data = as.data.frame(data)
    } else { # resolve formula from data.frame:
        if (!inherits(locations, "formula"))
            stop("locations argument should be a formula, such as ~x+y")
		m = model.frame(terms(locations), data)
        Terms = attr(m, "terms")
        attr(Terms, "intercept") = 0
        if ((yvar = attr(Terms, "response")) > 0)
            stop("no response allowed in locations formula")
        locations = model.matrix(Terms, m)
    }
    m = model.frame(terms(formula), data)
    Y = model.extract(m, response)
    if (length(Y) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)
    has.intercept = attr(Terms, "intercept")

    list(y = Y, locations = as.matrix(locations), X = X, call = call,
        has.intercept = has.intercept)
}
