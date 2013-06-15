print.dzim <- 
function(x, digits = max(3, getOption("digits") - 3), ...) {
  pX <- NCOL(x$X)
  z.value <- x$para / x$se 
  z.prob <- pvalue(z.value)
  coef <- data.frame(x$para, x$se, z.value, z.prob)  
  coefX <- switch(x$control$dist, 
    "poisson" = coef[1:pX, ], "nb" = coef[1 + 1:pX, ], 
    "zip" = coef[1 + 1:pX, ], "zinb" = coef[2 + 1:pX, ])
  coefAR<- switch(x$control$dist, 
    "poisson" = coef[-(1:pX), ], "nb" = coef[-(1:(1 + pX)), ], 
    "zip" = coef[-(1:(1 + pX)), ], "zinb" = coef[-(1:(2 + pX)), ])
  coefAR <- coefAR[1:x$control$order, ]
  colnames(coefX) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)") 
  colnames(coefAR) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(coefX) <- colnames(x$X)
  rownames(coefAR) <- paste("ar", 1:x$control$order, sep = "")  
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if(x$control$dist == "nb") {                   
    cat(paste("(Dispersion parameter for negative binomial taken to be ",
      round(x$para[1], 4), ")", sep = ""), "\n\n")
  } else if(x$control$dist == "zip") {
    cat(paste("(Zero-inflation parameter taken to be ",
      round(x$para[1], 4), ")", sep = ""), "\n\n")
  } else if(x$control$dist == "zinb") {
    cat(paste("(Zero-inflation parameter taken to be ",
      round(x$para[1], 4), ")", sep = ""), "\n\n")
    cat(paste("(Dispersion parameter for negative binomial taken to be ",
      round(x$para[2], 4), ")", sep = ""), "\n\n")
  } 
  cat("Coefficients (log-linear): \n")     
  printCoefmat(coefX, signif.legend = FALSE)
  cat("\n")
  cat("Coefficients (autoregressive): \n")
  printCoefmat(coefAR, signif.legend = FALSE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n") 
    cat(paste("(Standard deviation parameter taken to be ",
      round(x$para[length(x$para)], 4), ")", sep = ""), "\n\n")  
  cat("Criteria for assessing goodness of fit \n") 
  cat("loglik:", x$ps$pf$loglik, "\n")
  cat("aic:", x$aic, "\n")
  cat("bic:", x$bic, "\n") 
  cat("tic:", x$tic, "\n") 
  cat("\n")
  invisible(x)
}





