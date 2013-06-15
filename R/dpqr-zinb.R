dzinb <-
function(x, k, lambda, omega, log = FALSE) {
  d <- omega * (x == 0) + (1 - omega) * dnbinom(x, k, mu = lambda)
  if(log == FALSE) {
    d
  } else if(log == TRUE) {
    log(d)
  }
}

pzinb <-
function(q, k, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  p <- omega * (q >= 0) + (1 - omega) * pnbinom(q, k, mu = lambda)
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- log(p)
  }
  p
}

qzinb <-
function(p, k, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- exp(p)
  }
  qnbinom(pmax(0, (p - omega) / (1 - omega)), k, mu = lambda)
}

rzinb <-
function(n, k, lambda, omega) {
  ifelse(rbinom(n, 1, omega), 0, rnbinom(n, k, mu = lambda))
}

