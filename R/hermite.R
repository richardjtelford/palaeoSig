#' @importFrom stats dnorm

hermite <- function(x, theta) {
  herm <- function(kk, y) {
    hh <- list(kk + 1)
    hh[[0 + 1]] <- rep(1, length(y))
    hh[[1 + 1]] <- -y
    if (kk > 1) {
      for (k in 2:kk) {
        hh[[k + 1]] <- -1 / sqrt(k) * y * hh[[k + 1 - 1]] - sqrt((k - 1) / k) *
          hh[[k + 1 - 2]]
      }
    }
    hh
  }
  hh <- herm(length(theta) - 1, x)
  out <- 0
  for (k in seq_along(theta)) {
    # k is hermite_k plus 1
    out <- out + theta[k] * hh[[k]]
  }
  out
}
