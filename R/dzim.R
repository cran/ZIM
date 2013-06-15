dzim <-
function(formula, data, subset, na.action, weights = 1, offset = 0, 
  control = dzim.control(...), ...) {
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
  mf$formula <- formula
  mf <- eval(mf, parent.frame())    
  y <- round(model.response(mf, "numeric"))
  X <- model.matrix(terms(formula, data = data), mf)
  n <- NROW(X)
  weights <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  if(is.null(weights)) {
    weights <- rep(1, n)
  }
  if(is.null(offset)) {
    offset <- rep(0, n)
  }    
  fit <- dzim.fit(y, X, offset = offset, control = control)
  fit$call <- call
  fit$control <- control
  fit$na.action <- attr(mf, "na.action")
  fit$y <- y
  fit$X <- X  
  fit$aic <- (-2) * fit$ps$pf$loglik + 2 * length(fit$para)  
  fit$bic <- (-2) * fit$ps$pf$loglik + log(n) * length(fit$para) 
  if(is.numeric(try(solve(fit$deriv$info)))) {
    fit$se <- sqrt(diag(solve(fit$deriv$info)))
    if(any(is.na(fit$se))) {
      fit$se <- rep(NA, length(fit$para))    
    }
    fit$tic <- (-2) * fit$ps$pf$loglik + 
      2 * sum(diag(fit$deriv$J %*% solve(fit$deriv$info)))    
  } else {
    fit$se <- rep(NA, length(fit$para)) 
    fit$tic <- NA
  }
  names(fit$se) <- names(fit$para)      
  class(fit) <- "dzim"
  fit
}

