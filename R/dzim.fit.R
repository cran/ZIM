dzim.fit <- function(y, X, offset = rep(0, n), control = dzim.control(...), ...) {
  deriv <- function(para) {
    R <- control$R
    para <- switch(control$dist, "poisson" = c(0, Inf, para), 
      "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
    omega <- para[1]
    k <- para[2]
    beta <- para[2 + 1:pX]
    phi <- para[2 + pX + 1:p]
    sigma <- para[length(para)]
    m <- length(para)
    score <- rep(0, m)
    info.com <- matrix(0, m, m)
    info.mis <- matrix(0, m, m)
    info.obs <- matrix(0, m, m)
    U <- array(NA, dim = c(m, m, R))
    V <- array(NA, dim = c(m, m, R))
    W <- array(NA, dim = c(m, m, R))
    for(r in 1:R) {
      temp <- rep(0, m)
      v <- rep(0, m)
      K <- matrix(0, m, m)
      M <- matrix(0, m, m)  
      v[1] <- ifelse(omega == 0, 0, 
        ps$us[r, 1] / omega - (1 - ps$us[r, 1]) / (1 - omega))
      v[2] <- ifelse(k == Inf, 0, 
        1 + log(k) - digamma(k) + log(ps$vs[r, 1]) - ps$vs[r, 1])
      v[2 + 1:pX] <- (1 - ps$us[r, 1]) * (y[1] - 
        ps$vs[r, 1] * w[1] * exp(sum(X[1, ] * beta) + ps$ss[r, 1, 1])) * X[1, ] 
      v[2 + pX + 1:p] <- (1 / sigma^2) * 
        (ps$ss[r, 1, 1] - sum(phi * ps$ss0[r, ])) * ps$ss0[r, ]  
      v[m] <- (ps$ss[r, 1, 1] - 
        sum(phi * ps$ss0[r, ]))^2 / sigma^3 - 1 / sigma            
      K <- v %*% t(v)
      M[1, 1] <- ifelse(omega == 0, 0, 
        ps$us[r, 1] / omega^2 + (1 - ps$us[r, 1]) / (1 - omega)^2)
      M[2, 2] <- ifelse(k == Inf, 0, 
        trigamma(k) - 1 / k)
      M[2 + 1:pX, 2 + 1:pX] <- (1 - ps$us[r, 1]) * ps$vs[r, 1] * w[1] * 
        exp(sum(X[1, ] * beta) + ps$ss[r, 1, 1]) * X[1, ] %*% t(X[1, ])
      M[2 + pX + 1:p, 2 + pX + 1:p] <- (1 / sigma^2) * 
        ps$ss0[r, ] %*% t(ps$ss0[r, ])
      M[m, m] <- (3 / sigma^4) *  
        (ps$ss[r, 1, 1] - sum(phi * ps$ss0[r, ]))^2 - 1 / sigma^2
      M[2 + pX + 1:p, m] <- (2 / sigma^3) *
        (ps$ss[r, 1, 1] - sum(phi * ps$ss0[r, ])) * ps$ss0[r, ]               
      for(t in 2:n) {
        temp[1] <- ifelse(omega == 0, 0, 
          ps$us[r, t] / omega - (1 - ps$us[r, t]) / (1 - omega))
        temp[2] <- ifelse(k == Inf, 0, 
          1 + log(k) - digamma(k) + log(ps$vs[r, t]) - ps$vs[r, t])
        temp[2 + 1:pX] <- (1 - ps$us[r, t]) * (y[t] - 
          ps$vs[r, t] * w[t] * exp(sum(X[t, ] * beta) + ps$ss[r, t, 1])) * X[t, ] 
        temp[2 + pX + 1:p] <- (1 / sigma^2) * 
          (ps$ss[r, t, 1] - sum(phi * ps$ss[r, t - 1, ])) * ps$ss[r, t - 1, ] 
        temp[m] <- (ps$ss[r, t, 1] - 
          sum(phi * ps$ss[r, t - 1, ]))^2 / sigma^3 - 1 / sigma 
        v[1] <- v[1] + temp[1]
        v[2] <- v[2] + temp[2]
        v[2 + 1:pX] <- v[2 + 1:pX] + temp[2 + 1:pX]
        v[2 + pX + 1:p] <- v[2 + pX + 1:p] + temp[2 + pX + 1:p]
        v[m] <- v[m] + temp[m]
        K <- K + temp %*% t(temp)   
        M[1, 1] <- M[1, 1] + ifelse(omega == 0, 0, 
          ps$us[r, t] / omega^2 + (1 - ps$us[r, t]) / (1 - omega)^2)
        M[2, 2] <- M[2, 2] + ifelse(k == Inf, 0, 
          trigamma(k) - 1 / k)
        M[2 + 1:pX, 2 + 1:pX] <- M[2 + 1:pX, 2 + 1:pX] + 
          (1 - ps$us[r, t]) * ps$vs[r, t] * w[t] * 
          exp(sum(X[t, ] * beta) + ps$ss[r, t, 1]) * X[t, ] %*% t(X[t, ])
        M[2 + pX + 1:p, 2 + pX + 1:p] <- M[2 + pX + 1:p, 2 + pX + 1:p] + 
          (1 / sigma^2) * ps$ss[r, t - 1, ] %*% t(ps$ss[r, t - 1, ])
        M[m, m] <- M[m, m] + (3 / sigma^4) *  
          (ps$ss[r, t, 1] - sum(phi * ps$ss[r, t - 1, ]))^2 - 1 / sigma^2
        M[2 + pX + 1:p, m] <- M[2 + pX + 1:p, m] + (2 / sigma^3) *
          (ps$ss[r, t, 1] - sum(phi * ps$ss[r, t - 1, ])) * ps$ss[r, t - 1, ]       
      }
      M[m, 2 + pX + 1:p] <- t(M[2 + pX + 1:p, m])
      score <- score + v / R
      U[, , r] <- M
      V[, , r] <- v %*% t(v)
      W[, , r] <- K
    }
    gradient <- switch(control$dist, "poisson" = score[-(1:2)], 
      "nb" = score[-1], "zip" = score[-2], "zinb" = score)
    J <- apply(W, c(1, 2), mean)
    J <- switch(control$dist, "poisson" = J[-(1:2), -(1:2)], 
      "nb" = J[-1, -1], "zip" = J[-2, -2], "zinb" = J)
    J <- (J + t(J)) / 2
    info.com <- apply(U, c(1, 2), mean)
    info.mis <- apply(V, c(1, 2), mean) - score %*% t(score)
    for(prop in seq(0, 1, 0.01)) {
      info <- info.com - (1 - prop) * info.mis
      info <- switch(control$dist, "poisson" = info[-(1:2), -(1:2)], 
        "nb" = info[-1, -1], "zip" = info[-2, -2], "zinb" = info)
      info <- (info + t(info)) / 2 
      if(all(eigen(info)$values > 0)) {
        break
      }     
    }
    list(gradient = gradient, J = J, info = info)
  }  
  ar.fit <- function(ss, ss0) {
    ss.dim <- dim(ss)
    R <- ss.dim[1]
    n <- ss.dim[2]
    p <- ss.dim[3]
    A <- array(0, dim = c(R, p, p))
    b <- matrix(0, R, p)
    for(r in 1:R) {
        A[r, , ] <- ss0[r, ] %*% t(ss0[r, ])
        b[r, ] <- ss[r, 1, 1] * ss0[r, ]
      for(t in 2:n) {
        A[r, , ] <- A[r, , ] + ss[r, t - 1, ] %*% t(ss[r, t - 1, ])
        b[r, ] <- b[r, ] + ss[r, t, 1] * ss[r, t - 1, ]
      }
    }
    A <- apply(A, c(2, 3), mean)
    b <- colMeans(b)
    c <- colMeans(ss[, , 1]^2)
    phi <- solve(A) %*% b
    sigma <- sqrt((sum(c) - t(b) %*% solve(A) %*% b) / n)
    list(phi = c(phi), sigma = c(sigma))
  }
  nb.fit <- function(vs) {                     
    h <- function(p) {
      k <- p / (1 - p)
      e <- colMeans(vs)
      f <- colMeans(log(vs))
      abs(sum(1 + log(k) - digamma(k) + f - e))
    }
    p <- optimize(h, c(0, 1))$minimum
    p / (1 - p)
  }
  em <- function(para) {
    para <- switch(control$dist, "poisson" = c(0, Inf, para), 
      "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
    omega <- para[1]
    k <- para[2]
    beta <- para[2 + 1:pX]
    phi <- para[2 + pX + 1:p]
    sigma <- para[length(para)]
    d <- colMeans(ps$us)
    g <- colMeans((1 - ps$us) * ps$vs * exp(ps$ss[, , 1]))  
    omega <- ifelse(sum(control$dist == c("poisson", "nb")), 0, mean(d))
    k <- ifelse(sum(control$dist == c("poisson", "zip")), Inf, nb.fit(ps$vs))
    beta <- glm.fit(X, y, weights = 1 - d, 
     offset = ifelse(1 - d, log(g * w / (1 - d)), 0),
     family = poisson(), start = beta)$coef
    ar <- ar.fit(ps$ss, ps$ss0)  
    phi <- ar$phi
    sigma <- ar$sigma
    switch(control$dist,
      "poisson" = c(beta, phi, sigma),
      "nb" = c(k, beta, phi, sigma),
      "zip" = c(omega, beta, phi, sigma),
      "zinb" = c(omega, k, beta, phi, sigma))
  }
  n <- NROW(X)
  pX <- NCOL(X)
  w <- exp(offset)
  p <- control$order
  iter <- 0
  if(is.null(control$start)) {
    if(control$dist == "poisson") {
      fit1 <- glm(y ~ 0 + X + offset(offset), family = poisson)
      omega.old <- 0
      k.old <- Inf
      beta.old <- fit1$coef
    } else if(control$dist == "nb") {
      fit2 <- glm.nb(y ~ 0 + X + offset(offset))
      omega.old <- 0
      k.old <- fit2$theta
      beta.old <- fit2$coef
    } else if(control$dist == "zip") {
      fit3 <- zim(y ~ 0 + X + offset(offset) | 1, dist = "zip")
      omega.old <- plogis(fit3$para[pX + 1])
      k.old <- Inf
      beta.old <- fit3$para[1:pX]
    } else if(control$dist == "zinb") {
      fit4 <- zim(y ~ 0 + X + offset(offset) | 1, dist = "zinb")
      omega.old <- plogis(fit4$para[1 + pX + 1])  
      k.old <- exp(fit4$para[1])
      beta.old <- fit4$para[1 + 1:pX]    
    }
    phi.old <- rep(0, control$order)
    sigma.old <- 1    
  }
  para.old <- switch(control$dist,
    "poisson" = c(beta.old, phi.old, sigma.old),
    "nb" = c(k.old, beta.old, phi.old, sigma.old),
    "zip" = c(omega.old, beta.old, phi.old, sigma.old),
    "zinb" = c(omega.old, k.old, beta.old, phi.old, sigma.old))
  para.trace <- NULL
  loglik.trace <- NULL
  for(iter in 1:control$niter) {
    ps <- dzim.smooth(y, X, w, para.old, control)
    para.new <- em(para.old)
    para.old <- para.new
    loglik.old <- ps$pf$loglik 
    para.trace <- rbind(para.trace, para.old)
    loglik.trace <- c(loglik.trace, loglik.old) 
    if(control$trace == TRUE) {
      cat("iter =", iter, "\t loglik =", loglik.old, "\n")
    }
  }
  para <- para.new
  para.names <- c("omega", "k", paste("beta", as.character(1:pX - 1), sep = ""),
    paste("ar", as.character(1:control$order), sep = ""), "sigma")
  para.names <- switch(control$dist, "poisson" = para.names[-(1:2)], 
    "nb" = para.names[-1], "zip" = para.names[-2], "zinb" = para.names)
  colnames(para.trace) <- para.names
  names(para) <- para.names
  control$N <- 10 * control$N
  control$R <- 10 * control$R  
  ps <- dzim.smooth(y, X, w, para.new, control)
  deriv <- deriv(para.new)
  para.tmp <- switch(control$dist, "poisson" = c(0, Inf, para), 
    "nb" = c(0, para), "zip" = c(para[1], Inf, para[-1]), "zinb" = para)
  omega <- para.tmp[1]
  k <- para.tmp[2]
  beta <- para.tmp[2 + 1:pX]
  phi <- para.tmp[2 + pX + 1:p]
  sigma <- para.tmp[length(para)]
  lambda <- exp(offset + X %*% beta) * colMeans(exp(ps$ss[, , 1]))
  mu <- (1 - omega) * lambda
  names(deriv$gradient) <- names(para)
  rownames(deriv$J) <- names(para)
  colnames(deriv$J) <- names(para)
  rownames(deriv$info) <- names(para)
  colnames(deriv$info) <- names(para)
  list(omega = omega, k = k, beta = beta, phi = phi, sigma = sigma, 
    lambda = as.vector(lambda), mu = as.vector(mu),
    para = para, ps = ps, deriv = deriv,
    para.trace = ts(para.trace), loglik.trace = ts(loglik.trace))
}                               
