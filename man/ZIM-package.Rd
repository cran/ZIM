\name{ZIM-package}
\alias{ZIM}
\docType{package}

\title{
Statistical Models for Count Time Series with Excess Zeros
}

\description{
Fits observation-driven and parameter-driven models for count time series with excess zeros. 
}

\details{
\tabular{ll}{
Package: \tab ZIM        \cr
Type:    \tab Package    \cr
Version: \tab 1.0        \cr
Date:    \tab 2013-06-15 \cr
License: \tab GPL-2      \cr
}
The package \code{ZIM} contains functions to fit statistical models for count time series 
with excess zeros (Yang, 2012; Yang et al., 2013). The main function for fitting observation-driven models 
is \code{\link{zim}}, and the main function for fitting parameter-driven models is \code{\link{dzim}}. 
}

\author{
Ming Yang, Gideon K. D. Zamba, and Joseph E. Cavanaugh

Maintainer: Ming Yang <myang@sdac.harvard.edu>
}

\note{
The observation-driven models for zero-inflated count time series can also be fit using the function 
\code{\link[pscl]{zeroinfl}} from the \code{pscl} package (Zeileis et al., 2008). 
Fitting parameter-driven models is based on sequential Monte Carlo (SMC) methods, which are
computer intensive and could take several hours to estimate the model parameters.
}

\references{
Yang, M. (2012). Statistical models for count time series with excess zeros. 
\emph{PhD Dissertation}, The University of Iowa (\url{http://ir.uiowa.edu/etd/3019/}). 

Yang, M., Zamba, G. K. D., and Cavanaugh, J. E. (2013). Markov regression models for count time series
with excess zeros: A partial likelihood approach. \emph{Statistical Methodology}, \bold{14}:26-38.

Zeileis, A., Kleiber, C., and Jackman, S. (2008). Regression models for count data in \code{R}.  
\emph{Journal of Statistical Software}, \bold{27}(8). 
}

\keyword{package}
