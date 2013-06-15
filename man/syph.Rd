\name{syph}
\alias{syph}
\docType{data}

\title{
Example: Syphilis Series
}

\description{
Weekly number of syphilis cases in the United States from 2007 to 2010.
}

\format{
  A data frame with 209 observations on the following 69 variables.
  \describe{
    \item{\code{year}}{Year}
    \item{\code{week}}{Week}
    \item{\code{a1}}{\bold{United States}}
    \item{\code{a2}}{\bold{New England}}
    \item{\code{a3}}{Connecticut}
    \item{\code{a4}}{Maine}
    \item{\code{a5}}{Massachusetts}
    \item{\code{a6}}{New Hampshire}
    \item{\code{a7}}{Rhode Island}
    \item{\code{a8}}{Vermont}
    \item{\code{a9}}{\bold{Mid. Atlantic}}
    \item{\code{a10}}{New Jersey}
    \item{\code{a11}}{New York (Upstate)}
    \item{\code{a12}}{New York City}
    \item{\code{a13}}{Pennsylvania}
    \item{\code{a14}}{\bold{E.N. Central}}
    \item{\code{a15}}{Illinois}
    \item{\code{a16}}{Indiana}
    \item{\code{a17}}{Michigan}
    \item{\code{a18}}{Ohio}
    \item{\code{a19}}{Wisconsin}
    \item{\code{a20}}{\bold{W.N. Central}}
    \item{\code{a21}}{Iowa}
    \item{\code{a22}}{Kansas}
    \item{\code{a23}}{Minnesota}
    \item{\code{a24}}{Missouri}
    \item{\code{a25}}{Nebraska}
    \item{\code{a26}}{North Dakota}
    \item{\code{a27}}{South Dakota}
    \item{\code{a28}}{\bold{S. Atlantic}}
    \item{\code{a29}}{Delaware}
    \item{\code{a30}}{District of Columbia}
    \item{\code{a31}}{Florida}
    \item{\code{a32}}{Georgia}
    \item{\code{a33}}{Maryland}
    \item{\code{a34}}{North Carolina}
    \item{\code{a35}}{South Carolina}
    \item{\code{a36}}{Virginia}
    \item{\code{a37}}{West Virginia}
    \item{\code{a38}}{\bold{E.S. Central}}
    \item{\code{a39}}{Alabama}
    \item{\code{a40}}{Kentucky}
    \item{\code{a41}}{Mississippi}
    \item{\code{a42}}{Tennessee}
    \item{\code{a43}}{\bold{W.S. Central}}
    \item{\code{a44}}{Arkansas}
    \item{\code{a45}}{Louisana}
    \item{\code{a46}}{Oklahoma}
    \item{\code{a47}}{Texas}
    \item{\code{a48}}{\bold{Moutain}}
    \item{\code{a49}}{Arizona}
    \item{\code{a50}}{Colorado}
    \item{\code{a51}}{Idaho}
    \item{\code{a52}}{Montana}
    \item{\code{a53}}{Nevada}
    \item{\code{a54}}{New Mexico}
    \item{\code{a55}}{Utah}
    \item{\code{a56}}{Wyoming}
    \item{\code{a57}}{\bold{Pacific}}
    \item{\code{a58}}{Alaska}
    \item{\code{a59}}{California}
    \item{\code{a60}}{Hawaii}
    \item{\code{a61}}{Oregon}
    \item{\code{a62}}{Washington}
    \item{\code{a63}}{American Samoa}
    \item{\code{a64}}{C.N.M.I.}
    \item{\code{a65}}{Guam}
    \item{\code{a66}}{Peurto Rico}
    \item{\code{a67}}{U.S. Virgin Islands}
  }
}

\note{
C.N.M.I.: Commonwealth of Northern Mariana Islands.
}

\source{
CDC Morbidity and Mortality Weekly Report (\url{http://www.cdc.gov/MMWR/}).
}

\examples{
data(syph)
plot(ts(syph$a33), main = "Maryland")
}

\keyword{datasets}
