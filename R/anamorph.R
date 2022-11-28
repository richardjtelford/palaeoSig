#' @importFrom stats predict lm qnorm approxfun
#' @importFrom graphics plot lines points
#' @export

`anamorph` <- function(x, k, plot = FALSE) {
  theta <- hermite.coef(x, k)
  len <- length(x)
  gaus <- qnorm(seq(0, 1, length = len + 2)[-c(1, len + 2)])
  x_g <- hermite(gaus, theta)
  x_g_s <- sort(x_g)
  # extrapolation
  x_g_low <- x_g[1:10]
  x_g_high <- x_g[(len - 10):len]
  gaus_low <- gaus[1:10]
  gaus_high <- gaus[(len - 10):len]
  low <- predict(lm(x_g_low ~ gaus_low),
    newdata = data.frame(gaus_low = -10)
  )
  high <- predict(lm(x_g_high ~ gaus_high),
    newdata = data.frame(gaus_high = 10)
  )
  r <- range(x)
  if (plot) {
    plot(c(-10, gaus, 10), c(low, x_g, high),
      type = "l",
      xlim = range(gaus) + c(-1, +1), ylim = r + diff(r) / 4 * c(-1, 1)
    )
    lines(c(-10, gaus, 10), c(low, x_g_s, high), type = "l", col = 2)
    points(gaus, sort(x), col = 2)
  }
  xtog <- approxfun(c(x_g_s, low, high), c(gaus, -10, 10))
  gtox <- approxfun(c(gaus, -10, 10), c(x_g_s, low, high))
  list(xtog = xtog, gtox = gtox)
}
