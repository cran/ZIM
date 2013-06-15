\name{pvalue}
\alias{pvalue}

\title{Function to Compute P-value.}

\description{
Function to compute p-value based on a t-statistic.
}

\usage{
pvalue(t, df = Inf, alternative = c("two.sided", "less", "greater"))
}

\arguments{
\item{t}{t-statistic.}
\item{df}{degree of freedoms.}
\item{alternative}{type of alternative.}
}

\examples{
pvalue(1.96, alternative = "greater")
}

\keyword{misc}

