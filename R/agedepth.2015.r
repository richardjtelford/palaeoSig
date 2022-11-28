#' @export
#'
"agelme" <- function(depup, depdo, bpup, bpdo, use,
                     weights = c(1, rep(0, length(depup) - 1)),
                     vspan = 1, k = length(depup) - 1,
                     m = 2, diagnostic = FALSE) {
  # Cagedepth.fun version 1.4R, 26.09.05,
  # is written by Einar Heegaard, Bjerknes Centre for Climate Research,
  # University of Bergen, Allegaten 55, 5007 Bergen, Norway
  # e.mail: einar.heegaard@bio.uib.no
  # Use of this function is your responsibility
  if (missing(use)) {
    use <- TRUE
  }
  data <- data.frame(
    depthup = depup, depthdo = depdo,
    cageup = bpup, cagedo = bpdo, use = use
  )
  depup <- depup[use]
  depdo <- depdo[use]
  bpup <- bpup[use]
  bpdo <- bpdo[use]
  weights <- weights[use]
  n <- length(weights)
  x <- (depup + ((depdo - depup) / 2)) - min(depup)
  y <- ((bpup + bpdo) / 2) - min(bpup)
  sd <- abs(((bpup + bpdo) / 2) - bpup)
  fit_w <- 1 / sd
  fit_w[weights == 1] <- 1
  fit_con <- gam(y ~ s(x, k = k, m = m),
    quasi(link = "identity", variance = "constant"),
    weights = drop(fit_w), scale = -1
  )
  fit_mu <- gam(y ~ s(x, k = k, m = m),
    quasi(link = "identity", variance = "mu"),
    weights = drop(fit_w), scale = -1
  )
  fit_pl <- list(constant = fit_con, mu = fit_mu)
  if (diagnostic) {
    par(mfrow = c(3, 4))
    v <- unlist(lapply(fit_pl, function(v) {
      v[[2]]
    }))
    x1 <- range(v)
    x2 <- range(sqrt(abs(v)))
    for (i in 1:2) {
      x3 <- fit_pl[[i]][[3]]
      x4 <- fit_pl[[i]][[2]]
      x5 <- fit_pl[[i]]$y
      plot(x3, x4, ylim = x1, xlab = "Fitted", ylab = "Residuals")
      lines(loess.smooth(x3, x4, span = vspan))
      abline(h = 0, lty = 2)
      plot(x3, sqrt(abs(x4)),
        ylim = x2,
        ylab = "Sqrt(abs(Residuals))", xlab = "Fitted"
      )
      lines(loess.smooth(x3, sqrt(abs(x4)), span = vspan))
      abline(h = 0, lty = 2)
      plot(x3, x5, ylab = "Observed", xlab = "Fitted")
      lines(c(0, max(x5)), c(0, max(x5)), lty = 2)
      lines(loess.smooth(x3, x5, span = vspan))
      qqnorm(x4)
      qqline(x4, lty = 2)
    }
  }
  age_res <- list(
    tdf = c(sum(fit_con$hat), sum(fit_mu$hat)),
    weights = fit_w,
    RES = data.frame(
      Constvar = sum(fit_con[[2]]^2) / 1000,
      Muvar = sum(fit_mu[[2]]^2) / 1000
    ),
    Models = fit_pl,
    data = data
  )
  class(age_res) <- "agelme"
  age_res
}

#' @export
"predict.agelme" <- function(object, v = 1, depth, ...) {
  # Cagenew.fun version 1.4R,26.09.05,
  # written by Einar Heegaard, Bjerknes Centre for Climate Research,
  # University of Bergen, Allegaten 55, 5007 Bergen, Norway
  # e.mail: einar.heegaard@bio.uib.no
  # Use of this function is your responibility
  depup <- object$data[object$data$use, 1]
  depdo <- object$data[object$data$use, 2]
  bpup <- object$data[object$data$use, 3]
  bpdo <- object$data[object$data$use, 4]
  fit_m <- object$Models[[v]]
  xa <- data.frame(fit_w = object$weights)
  xd <- depth - min(depup)

  yp <- predict(fit_m, newdata = data.frame(x = xd), type = "response")
  ypm <- predict(fit_m, se.fit = TRUE)

  sd <- abs(((bpup + bpdo) / 2) - bpup)
  x <- (depup + ((depdo - depup) / 2))
  sd2 <- sqrt(sd^2 + ypm$se.fit^2)
  sdp <- approx(x, sd2, depth, rule = 2)
  sdp <- sdp$y
  y1 <- yp - (1.96 * sdp)
  y2 <- yp + (1.96 * sdp)
  vt <- data.frame(
    depth = depth,
    Estage = yp + min(bpup),
    Lowlim = y1 + min(bpup),
    Upplim = y2 + min(bpup),
    Tsd = sdp
  )
  vt <- list(
    v = ifelse(v == 1, "Constant variance", "Constant variance"),
    fit = vt,
    data = object$data
  )
  class(vt) <- "fittedAgelme"
  vt
}
