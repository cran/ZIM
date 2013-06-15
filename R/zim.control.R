zim.control <-
function(dist = c("zip", "zinb"), method = c("EM-NR", "EM-FS"), 
  type = c("solve", "ginv"), robust = FALSE, trace = FALSE, start = NULL, 
  minit = 10, maxit = 10000, epsilon = 1e-8) {
  dist <- match.arg(dist)
  method <- match.arg(method)
  type <- match.arg(type)
  inv <- switch(type, "solve" = solve, "ginv" = ginv)
  list(dist = dist, method = method, type = type, inv = inv, robust = robust, 
    trace = trace, start = start, minit = minit, maxit = maxit, epsilon = epsilon)
}

