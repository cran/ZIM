dzim.sim <-  function(X, w, omega, k, beta, phi, sigma, mu0, Sigma0) {
  n <- NROW(X)
  p <- length(phi)
  s0 <- mvrnorm(1, mu0, Sigma0)
  s <- matrix(NA, n, p)
  u <- rep(NA, n) 
  v <- rep(NA, n) 
  y <- rep(NA, n)
  Phi <- suppressWarnings(rbind(phi, cbind(diag(1, p - 1), 0))) 
  s[1, ] <- Phi %*% s0 +  c(rnorm(1, 0, sigma), rep(0, p - 1))
  u[1] <- rbinom(1, 1, omega)
  v[1] <- ifelse(k == Inf, 1, rgamma(1, shape = k, scale = 1 / k))
  lambda <- w[1] * exp(sum(X[1, ] * beta) + s[1, 1])                        
  y[1] <- rpois(1, (1 - u[1]) * v[1] * lambda)
  for(t in 2:n) {
    s[t, ] <- Phi %*% s[t - 1, ] + c(rnorm(1, 0, sigma), rep(0, p - 1))
    u[t] <- rbinom(1, 1, omega)
    v[t] <- ifelse(k == Inf, 1, rgamma(1, shape = k, scale = 1 / k))
    lambda <- w[t] * exp(sum(X[t, ] * beta) + s[t, 1])
    y[t] <- rpois(1, (1 - u[t]) * v[t] * lambda)
  }
  list(s0 = s0, s = s, u = u, v = v, y = y)
}
