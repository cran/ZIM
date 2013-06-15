dzim.plot <- function(object, k.inv = FALSE, sigma.sq = FALSE, ...) {
  pX <- NCOL(object$X)
  control <- object$control
  para.trace <- object$para.trace
  I <- NROW(para.trace)
  J <- NCOL(para.trace)
  if(k.inv == TRUE) {
    if(control$dist == "nb")   para.trace[, 1] <- 1 / para.trace[, 1]
    if(control$dist == "zinb") para.trace[, 2] <- 1 / para.trace[, 2]
  }
  if(sigma.sq == TRUE) {
    para.trace[, NCOL(para.trace)] <- para.trace[, NCOL(para.trace)]^2
  }
  scale <- apply(para.trace[round(I/2):I, ], 2, sd)
  scale <- ifelse(scale == 0, 1, scale)
  para.trace.plot <- (para.trace -  matrix(rep(para.trace[1, ], I), nrow = I, byrow = TRUE)) %*% solve(diag(scale))
  col <- c(rep("green", pX), rep("blue", control$order), "purple")
  lty <- c(1:pX, 1:control$order, 1)
  col <- switch(control$dist, "poisson" = col, 
    "nb" = c("red", col), "zip" = c("gold", col), "zinb" = c("gold", "red", col))
  lty <- switch(control$dist, "poisson" = lty, 
    "nb" = c(1, lty), "zip" = c(1, lty), "zinb" = c(1, 1, lty))  
  plot.ts(para.trace.plot, plot.type = "single", lty = lty, col = col, yaxt = "n", xlab = "Iteration", ylab = " ", ...)
  points(1, 0, pch = 20)
}

