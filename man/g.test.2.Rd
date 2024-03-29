% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{g.test.2}
\alias{g.test.2}
\title{Log-likelihood test of independence.}
\usage{
g.test.2(obsx, obsy)
}
\arguments{
\item{obsx}{a numeric vector of positive counts, with at least one non-zero value.}

\item{obsy}{a numeric vector of positive counts of same length as obsx, with at least one non-zero value.

The order of values in the two vectors should be the same.
No corrections are applied.
No input checks are applied, as RATs needs to run this millions of times.}
}
\description{
For two sets of observations.
}
\details{
General utility function.
}
