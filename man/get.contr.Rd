% $Id: get.contr.Rd,v 1.5 2007-11-16 12:59:47 edzer Exp $
\name{get.contr}
\alias{get.contr}
\title{ Calculate contrasts from multivariable predictions }
\description{ Given multivariable predictions and prediction (co)variances,
calculate contrasts and their (co)variance }
\usage{ get.contr(data, gstat.object, X, ids = names(gstat.object$data)) }
\arguments{
\item{data}{data frame, output of \link{predict.gstat} }
\item{gstat.object}{object of class \code{gstat}, used to
extract ids; may be missing if \code{ids} is used }
\item{X}{ contrast vector or matrix; the number of variables in
\code{gstat.object} should equal the number of elements in \code{X}
if \code{X} is a vector, or the number of rows in \code{X} if \code{X}
is a matrix. }
\item{ids}{ character vector with (selection of) id names, present in data }
}

\details{ 
From data, we can extract the $(n \times 1)$ vector with multivariable
predictions, say $y$, and its $(n \times n)$ covariance matrix $V$. Given
a contrast matrix in $X$, this function computes the contrast vector is
$C=X'y$ and $Var(C)=X'V X$.
}

\value{
a data frame containing for each row in \code{data} the generalized
least squares estimates (named beta.1, beta.2, ...), their
variances (named var.beta.1, var.beta.2, ...) and covariances
(named cov.beta.1.2, cov.beta.1.3, ...)
}

\references{ 
\url{http://www.gstat.org/}
}
\author{ Edzer J. Pebesma }
\seealso{\link{predict.gstat}}


\keyword{ models }
