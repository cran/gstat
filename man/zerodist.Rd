\name{zerodist}
\alias{zerodist}
\title{ find point pairs with equal spatial coordinates }
\description{ find point pairs with equal spatial coordinates }
\usage{ zerodist(x, y, z, zero = 0.0) }
\arguments{
\item{x}{ vector with x-coordinate }
\item{y}{ vector with y-coordinate (may be missing) }
\item{z}{ vector with z-coordinate (may be missing) }
\item{zero}{ value to be compared to for establishing when a distance
is considered zero (default 0.0) }
}
\value{pairs of row numbers with identical coordinates, numeric(0) if no
such pairs are found }
\note{ Duplicate observations sharing identical spatial locations
result in singular covariance matrices in kriging situations.
This function may help identifying spatial duplications, so they
can be removed.  A matrix with all pair-wise distances is calculated,
so if x, y and z are large this function is slow }

\keyword{dplot}
