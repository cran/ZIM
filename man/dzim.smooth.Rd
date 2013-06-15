\name{dzim.smooth}
\alias{dzim.smooth}

\title{Particle Smoothing for DZIM}

\description{
Function to implement the particle smoothing method proposed by Gordsill et al. (2004).
}

\usage{
dzim.smooth(y, X, w, para, control)
}

\arguments{
\item{y}{response variable.}
\item{X}{design matrix.}
\item{w}{\code{log(w)} is used as an offset variable in the linear predictor.}
\item{para}{model parameters.}
\item{control}{control arguments.}
}

\seealso{
\code{\link{dzim}},
\code{\link{dzim.fit}}, 
\code{\link{dzim.filter}},
\code{\link{dzim.control}}, 
\code{\link{dzim.sim}}, 
\code{\link{dzim.plot}}
}

\references{
Gordsill, S. J., Doucet, A., and West, M. (2004). Monte Carlo smoothing for nonlinear time series. 
\emph{Journal of the American Statistical Association}, \bold{99}, 156-168. 
}

\keyword{regression}

