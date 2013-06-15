\name{dzim.plot}
\alias{dzim.plot}

\title{Trace Plots from DZIM}

\description{
Function to display trace plots from a dynamic zero-inflated model.
}

\usage{
dzim.plot(object, k.inv = FALSE, sigma.sq = FALSE, ...) 
}

\arguments{
\item{object}{objective from \code{\link{dzim}} or \code{\link{dzim.fit}}.}
\item{k.inv}{logical; indicating whether an inverse transformation is needed for the dispersion parameter.}
\item{sigma.sq}{logical; indicating whether a square transformation is needed for the standard deviation parameter.}
\item{...}{additional arguments.}
}

\seealso{
\code{\link{dzim}},
\code{\link{dzim.fit}},
\code{\link{dzim.filter}}, 
\code{\link{dzim.smooth}},
\code{\link{dzim.control}},
\code{\link{dzim.sim}}
}

\keyword{regression}

