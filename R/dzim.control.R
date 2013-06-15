dzim.control <- function(dist = c("poisson", "nb", "zip", "zinb"),
  trace = FALSE, start = NULL, order = 1, 
  mu0 = rep(0, order), Sigma0 = diag(1, order), 
  N = 1000, R = 1000, niter = 500) {
  dist <- match.arg(dist)
  list(dist = dist, trace = trace, start = start,
    order = order, mu0 = mu0, Sigma0 = Sigma0,
    N = N, R = R, niter = niter)
}
