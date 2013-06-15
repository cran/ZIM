dzim.smooth <- function(y, X, w, para, control) {
  n <- NROW(X)
  pX <- NCOL(X)
  p <- control$order
  N <- control$N
  R <- control$R
  mu0 <- control$mu0
  Sigma0 <- control$Sigma0
  pf <- dzim.filter(y, X, w, para, control)
  para <- switch(control$dist, "poisson" = c(0, Inf, para), 
    "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
  omega <- para[1]
  k <- para[2]
  beta <- para[2 + 1:pX]
  phi <- para[2 + pX + 1:p]
  sigma <- para[length(para)] 
  qs <- matrix(NA, N, n) 
  ss <- array(NA, dim = c(R, n, p)) 
  us <- matrix(NA, R, n) 
  vs <- matrix(NA, R, n)
  qs0 <- rep(NA, N)
  ss0 <- matrix(NA, R, p) 
  for(r in 1:R) {
    qs[, n] <- pf$qf[, n]
    index <- sample(1:N, size = 1, prob = qs[, n])
    ss[r, n, ] <- pf$sf[index, n, ]
    us[r, n] <- pf$uf[index, n]
    vs[r, n] <- pf$vf[index, n]
    for(t in (n - 1):1) {
      qs[, t] <- pf$qf[, t] *
        dnorm(ss[r, t + 1, 1], as.matrix(pf$sf[, t, ]) %*% phi, sigma) *
        dbinom(us[r, t + 1], 1, omega) *
        ifelse(k == Inf, 1, dgamma(vs[r, t + 1], shape = k, scale = 1 / k))
      index <- sample(1:N, size = 1, prob = qs[, t])
      ss[r, t, ] <- pf$sf[index, t, ]
      us[r, t] <- pf$uf[index, t]
      vs[r, t] <- pf$vf[index, t]    
    }
    qs0 <- pf$qf0 *
        dnorm(ss[r, 1, 1], as.matrix(pf$sf0) %*% phi, sigma) *
        dbinom(us[r, 1], 1, omega) *
        ifelse(k == Inf, 1, dgamma(vs[r, 1], shape = k, scale = 1 / k))
    index <- sample(1:N, size = 1, prob = qs0)    
    ss0[r, ] <- pf$sf0[index, ]  
  }
  list(pf = pf, qs = qs, ss = ss, us = us, vs = vs, qs0 = qs0, ss0 = ss0)  
}
