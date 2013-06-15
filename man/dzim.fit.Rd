\name{dzim.fit}
\alias{dzim.fit}

\title{Fitter Function for Dynamic Zero-Inflated Models}

\description{
\code{\link{dzim.fit}} is the basic computing engine called by \code{\link{dzim}} used to fit dynamic zero-inflated models. 
This should usually \emph{not} be used directly unless by experienced users. 
}

\usage{
dzim.fit(y, X, offset = rep(0, n), control = dzim.control(...), ...)
}

\arguments{
\item{y}{response variable.}
\item{X}{design matrix.}
\item{offset}{offset variable.}
\item{control}{control arguments.}
\item{...}{additional arguments.}
}

\seealso{
\code{\link{dzim}},
\code{\link{dzim.control}}, 
\code{\link{dzim.filter}}, 
\code{\link{dzim.smooth}},
\code{\link{dzim.sim}},
\code{\link{dzim.plot}}
}

\keyword{regression}

