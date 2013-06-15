dzim.filter <- function(y, X, w, para, control) {
  n <- NROW(X)
  pX <- NCOL(X)
  p <- control$order
  N <- control$N
  R <- control$R
  mu0 <- control$mu0
  Sigma0 <- control$Sigma0
  para <- switch(control$dist, "poisson" = c(0, Inf, para), 
    "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
  omega <- para[1]
  k <- para[2]
  beta <- para[2 + 1:pX]
  phi <- para[2 + pX + 1:p]
  sigma <- para[length(para)]
  sp <- array(NA, dim = c(N, n, p)) 
  up <- matrix(NA, N, n)
  vp <- matrix(NA, N, n)
  qf <- matrix(NA, N, n) 
  sf <- array(NA, dim = c(N, n, p)) 
  uf <- matrix(NA, N, n) 
  vf <- matrix(NA, N, n)
  Phi <- suppressWarnings(rbind(phi, cbind(diag(1, p - 1), 0)))  
  qf0 <- rep(1, N)
  sf0 <- mvrnorm(N, mu0, Sigma0)
  for(i in 1:N) {
    sp[i, 1, ] <- Phi %*% sf0[i, ] + c(rnorm(1, 0, sigma), rep(0, p - 1))
  }
  up[, 1] <- rbinom(N, 1, omega)
  vp[, 1] <- ifelse(rep(k == Inf, N), rep(1, N), rgamma(N, shape = k, scale = 1 / k))
  lambda <- w[1] * exp(sum(X[1, ] * beta) + sp[, 1, 1])
  qf[, 1] <- dpois(y[1], (1 - up[, 1]) * vp[, 1] * lambda)
  index <- suppressWarnings(sample(1:N,
    size = N, replace = TRUE, prob = qf[, 1]))
  sf[, 1, ] <- sp[index, 1, ]
  uf[, 1] <- up[index, 1]
  vf[, 1] <- vp[index, 1]
  for(t in 2:n) {
    for(i in 1:N) {
      sp[i, t, ] <- Phi %*% sf[i, t - 1, ] + c(rnorm(1, 0, sigma), rep(0, p - 1))    
    }
    up[, t] <- rbinom(N, 1, omega)
    vp[, t] <- ifelse(rep(k == Inf, N), rep(1, N), rgamma(N, shape = k, scale = 1 / k))
    lambda <- w[t] * exp(sum(X[t, ] * beta) + sp[, t, 1])
    qf[, t] <- dpois(y[t], (1 - up[, t]) * vp[, t] * lambda)
    index <- suppressWarnings(sample(1:N,
      size = N, replace = TRUE, prob = qf[, t]))
    sf[, t, ] <- sp[index, t, ]
    uf[, t] <- up[index, t]
    vf[, t] <- vp[index, t]
  }
  loglik <- sum(log(colMeans(qf)))
  list(qf0 = qf0, sf0 = sf0, sp = sp, up = up, vp = vp, 
    qf = qf, sf = sf, uf = uf, vf = vf, loglik = loglik)
}
