pvalue <- function(t, df = Inf, alternative = c("two.sided", "less", "greater")) {
  alternative <- match.arg(alternative)
  if(alternative == "less") {
    p <- pt(t, df = df)
  } else if(alternative == "greater") {
    p <- pt(t, df = df, lower.tail = FALSE)
  } else {
    p <- 2 * pt(-abs(t), df = df)
  }
  p
}
