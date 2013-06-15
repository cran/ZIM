print.zim <- 
function(x, digits = max(3, getOption("digits") - 3), ...) {
  pX <- NCOL(x$X)
  pZ <- NCOL(x$Z)
  z.value <- x$para / x$se 
  z.prob <- pvalue(z.value)
  coef <- data.frame(x$para, x$se, z.value, z.prob)
  coefX <- coef[(x$control$dist == "zinb") + 1:pX, ]
  coefZ <- coef[(x$control$dist == "zinb") + pX + 1:pZ, ]
  colnames(coefX) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)") 
  colnames(coefZ) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")  
  rownames(coefX) <- colnames(x$X)
  rownames(coefZ) <- colnames(x$Z)     
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")  
  cat("Coefficients (log-linear): \n")   
  printCoefmat(coefX, signif.legend = FALSE)
  cat("\n")
  cat("Coefficients (logistic): \n")     
  printCoefmat(coefZ, signif.legend = FALSE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n")  
  if(x$control$dist == "zip") {                   
    cat("Test for overdispersion (H0: ZIP vs. H1: ZINB) \n")
    cat("score.test:", round(x$score.test, 4), "\n")
    cat("p.value:", format.pval(x$p.value), "\n\n")
  } else {
    cat(paste("(Dispersion parameter for negative binomial taken to be ",
      round(x$k, 4), ")", sep = ""), "\n\n")
  }
  cat("Criteria for assessing goodness of fit \n") 
  cat("loglik:", x$loglik, "\n")
  cat("aic:", x$aic, "\n")
  cat("bic:", x$bic, "\n") 
  cat("tic:", x$tic, "\n") 
  cat("\n")
  cat("Number of", x$control$method, "iterations:", x$iter, "\n")
  cat("Maximum absolute gradient:", max(abs(x$gradient)), "\n")
  cat("\n")
  invisible(x)
}





