#' @title Estimation of the relationship between Calibrated age and depth
#' @description
#' Estimates the relationship of Calibrated age and depth for palaeorecords.
#' The function uses a smooth spline from the mgcv library by Simon Wood.
#' It produces predicted confidence interval for the relationship approximating
#' a mixed effect model, as there are two levels of uncertainty,
#' i.e. within dated object and between dated objects.
#'
#' @param depup The upper depths of the dated slides
#' @param depdo The lower depths of the dated slides
#' @param bpup The younger calibrated ages of the dated slides
#' @param bpdo The older calibrated ages of the dated slides
#' @param use Logical vector of dates to include in the model.
#' Default is to use all.
#' @param weights Weights to be used for the estimation, default is fixed
#' top-layer followed by inverse variance of within dated object
#' @param vspan The span to be used for the diagnostic plots, default span = 1
#' @param k Number of base function to start the shrinkage in the
#' gam estimation procedure
#' @param m The order of penalty for the term, i.e. the degree of continuity at
#' the knots (default, m = 2 gives cubic smooth spline)
#' @param diagnostic Logical, should diagnostic plots be made.
#'
#' @details
#' Note that the fixation of the top layer is done by a weight = 1, whereas the
#' other weights follows inverse variance within object.
#'
#' The diagnostic plots is used to check the quality of the estimation and to
#' see if there is a need for an assumption of between object variance
#' proportional to mean.
#' The latter however is rarely encountered for palaeodata.
#'
#' @returns
#' A list of class 'agelme' with the following items:
#'
#' \item{tdf }{Degrees of freedom used by the cubic smooth spline, a vector with first value for constant variance and second vector for variance equal to mu.}
#' \item{weights }{A vector of the weights used by the cubic smooth spline }
#' \item{Constant}{A matrix with the numerical results for the dated points using a constant variance}
#' \item{Muvar}{A matrix with the numerical results for the dated points using variance equal to mu}
#' \item{RES}{A vector of the Residual sum of squares}
#' \item{Models}{A list with the models from the cubic smooth spline, constant and mu variance, respectively}
#' \item{Data}{A data.frame including the data used for the estimation}
#'
#' @references
#' Heegaard, E., Birks, HJB. & Telford, RJ. 2005. Relationships between
#' calibrated ages and depth in stratigraphical sequences: an estimation
#' procedure by mixed-effect regression. _The Holocene_ **15**: 612-618
#' \doi{10.1191/0959683605hl836rr}
#'
#' @author Einar Heegaard
#' @importFrom mgcv gam s
#' @examples
#' data(STOR)
#'
#' fit.mod <- with(STOR, agelme(depthup, depthdo, cageup, cagedo))
#'
#' # Predicting using the constant variance model,
#' # for each cm between 70 and 400 cm.
#' fit.pre <- predict(fit.mod, 1, 70:400)
#' plot(fit.pre)
#' @export

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


#############################################
#' @title Predicts the Calibrated age for agelme models
#' @description
#' This function uses the output from 'agelme' to predict the Calibrated ages
#' for specified depths.
#'
#' @param object An `agelme` model
#' @param v Using constant (1) or mu (2) variance
#' @param depth A vector of the depths to be predicted
#' @param \dots Other arguments, currently unused.
#' @returns
#'  A list with three items
#'  \describe{
#'  \item{v}{ Whether constant variance or mu variance used.}
#'  \item{fit}{ A data.frame with the model predictions, containing columns:
#' \describe{
#'   \item{Depth}{ The depths for the predicted ages}
#'   \item{Estage}{ Predicted age}
#'   \item{Lowlim}{ Lower 95 % confidence interval}
#'   \item{Upplim}{ Upper 95 % confidence interval}
#'   \item{Tsd}{ Total standard deviation}
#' }
#' }
#' }
#' \item{data}{ A data.frame containing the age and depth information of the radiocarbon dates.}
#' @references
#' Heegaard, E., Birks, HJB. & Telford, RJ. 2005.
#' Relationships between calibrated ages and depth in stratigraphical sequences:
#' an estimation procedure by mixed-effect regression. _The Holocene_ **15**: 612-618
#' \doi{10.1191/0959683605hl836rr}
#' @author Einar Heegaard
#' @examples
#' data(STOR)
#'
#' fit.mod <- with(STOR, agelme(depthup, depthdo, cageup, cagedo))
#'
#' # Predicting using the constant variance model,
#' # for each cm between 70 and 400 cm.
#' fit.pre <- predict(fit.mod, 1, 70:400)
#' plot(fit.pre)
#' @export

"predict.agelme" <- function(object, v = 1, depth, ...) {
  # Cagenew.fun version 1.4R,26.09.05,
  # written by Einar Heegaard, Bjerknes Centre for Climate Research,
  # University of Bergen, Allegaten 55, 5007 Bergen, Norway
  # e.mail: einar.heegaard@bio.uib.no
  # Use of this function is your responsibility
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
