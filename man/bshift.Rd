\name{bshift}
\alias{bshift}

\title{Backshift Operator}

\description{
Apply the backshift or lag operator to a time series objective.
}

\usage{
bshift(x, k = 1)
}

\arguments{
\item{x}{univariate or multivariate time series.}
\item{k}{number of lags.}
}

\seealso{
\code{\link{lag}}, \code{\link[TSA]{zlag}} 
}

\examples{
x <- arima.sim(model = list(ar = 0.8, sd = 0.5), n = 120)
bshift(x, k = 12)
}

\keyword{misc}

