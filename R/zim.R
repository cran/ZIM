zim <-
function(formula, data, subset, na.action, weights = 1, offset = 0, 
  control = zim.control(...),  ...) {
  call <- match.call()
  if(missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action",
    "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf[[1]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE 
  if(any(as.character(formula[[3]]) == "|")) {
    formulaX <- . ~ .
    formulaZ <- . ~ .
    formulaX[[2]] <- formula[[2]]
    formulaZ[[2]] <- formula[[2]]
    formulaX[[3]] <- formula[[3]][[2]]
    formulaZ[[3]] <- formula[[3]][[3]]
    formula[[3]][[1]] <- as.name("+")
    mf$formula <- formula
  } else {
    formulaX <- formula
    formulaZ <- formula
    mf$formula <- formula
  }
  mf <- eval(mf, parent.frame())    
  y <- round(model.response(mf, "numeric"))
  X <- model.matrix(terms(formulaX, data = data), mf)
  Z <- model.matrix(terms(formulaZ, data = data), mf)
  nobs <- length(y)
  pX <- NCOL(X)
  pZ <- NCOL(Z)
  y0 <- (y == 0)
  weights <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  if(is.null(weights)) {
    weights <- rep(1, nobs)
  }
  if(is.null(offset)) {
    offset <- rep(0, nobs)
  }    
  fit <- zim.fit(y, X, Z, weights = weights, offset = offset, control = control)
  fit$call <- call
  fit$control <- control
  fit$na.action <- attr(mf, "na.action")
  fit$y <- y
  fit$X <- X
  fit$Z <- Z                         
  if(control$dist == "zip") {
    lambda <- exp(X %*% fit$para[1:pX] + offset)
    omega <- plogis(Z %*% fit$para[pX + 1:pZ]) 
    p0 <- dzip(0, lambda, omega)
    J <- matrix(NA, 1 + pX + pZ, 1 + pX + pZ)
    J[1, 1] <- sum(weights * (lambda^2 * (2 * (1 - omega) -
      omega * lambda^2 * (1 - omega / p0)))) / 4
    J[1, 1 + 1:pX] <- t(X) %*% (weights * (omega *
      lambda^3 * (1 - omega / p0))) / 2
    J[1, 1 + pX + 1:pZ] <- t(Z) %*% (weights * (omega *
      lambda^2 * (1 - omega / p0))) / 2
    J[1 + 1:pX, 1 + 1:pX] <- t(X) %*% (weights * as.vector(lambda *
      ((1 - omega) - omega * lambda * (1 - omega / p0))) * X)
    J[1 + 1:pX, 1 + pX + 1:pZ] <- (-1) * t(X) %*% (weights *
      as.vector(omega * lambda * (1 - omega / p0)) * Z)
    J[1 + pX + 1:pZ, 1 + pX + 1:pZ] <- t(Z) %*% (weights *
      as.vector(omega^2 * (1 / p0 - 1)) * Z)
    J[1 + 1:pX, 1] <- t(J[1, 1 + 1:pX])
    J[1 + pX + 1:pZ, 1] <- t(J[1, 1 + pX + 1:pZ])
    J[1 + pX + 1:pZ, 1 + 1:pX] <- t(J[1 + 1:pX, 1 + pX + 1:pZ])
    J <- (J + t(J)) / 2
    S <- sum(weights * ((y - lambda)^2 - y - y0 * omega * lambda^2 / p0)) / 2  
    fit$score.test <- S * sqrt(control$inv(J)[1, 1])
    fit$p.value <- pnorm(fit$score.test, lower.tail = FALSE)
    fit$k <- Inf
    fit$lambda <- lambda
    fit$omega <- omega     
  } else {
    fit$k <- exp(fit$para[1])
    fit$lambda <- exp(X %*% fit$para[1 + 1:pX] + offset)
    fit$omega <- plogis(Z %*% fit$para[1 + pX + 1:pZ]) 
  }
  fit$fitted.values <- fit$lambda * (1 - fit$omega)
  fit$residuals <- (y - fit$fitted.values) / sqrt(fit$fitted.values * 
    (1 + fit$omega * fit$lambda + fit$lambda / fit$k))             
  class(fit) <- "zim"
  fit
}

