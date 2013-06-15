dzip <-
function(x, lambda, omega, log = FALSE) {
  d <- omega * (x == 0) + (1 - omega) * dpois(x, lambda) 
  if(log == FALSE) {
    d
  } else if(log == TRUE) {
    log(d)
  }              
}

pzip <-
function(q, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  p <- omega * (q >= 0) + (1 - omega) * ppois(q, lambda)
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- log(p)
  }
  p
}

qzip <-
function(p, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- exp(p)
  }
  qpois(pmax(0, (p - omega) / (1 - omega)), lambda)
}

rzip <-
function(n, lambda, omega) {
  ifelse(rbinom(n, 1, omega), 0, rpois(n, lambda))
}



