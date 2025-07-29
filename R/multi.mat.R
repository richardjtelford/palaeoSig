#' @title MAT for multiple variables
#' @description
#' MAT for many environmental variables simultaneously.
#' More efficient than calculating them separately for each variable.
#' @param training.spp Community data
#' @param envs Environmental variables - or simulations
#' @param core.spp Optional fossil data to make predictions for
#' @param noanalogues Number of analogues to use
#' @param method distance metric to use
#' @param run Return LOO predictions or predictions for fossil data
#' @returns If `run = "both"`, a list with two elements:
#' \item{jack}{ Matrix of leave-one-out cross-validation predictions for the
#' calibration set}
#' \item{core}{ Matrix of predictions for the fossil data}
#' Otherwise, one of these matrices is returned.
#'
#' @references
#' Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in
#' spatially structured environments. \emph{Quaternary Science Reviews}
#' \bold{28}: 1309--1316. \doi{10.1016/j.quascirev.2008.12.020}
#' @author Richard Telford \email{Richard.Telford@bio.uib.no}
#' @examples
#' data(arctic.env)
#' data(arctic.pollen)
#'
#' mMAT <- multi.mat(arctic.pollen, arctic.env[, 9:67], noanalogues = 5)

#' @keywords multivariate


#' @export

`multi.mat` <- function(training.spp, envs, core.spp, noanalogues = 10,
                        method = "sq-chord", run = "both") {
  ests <- function(d.mat, nRow, nSamp) {
    lapply(noanalogues, function(ana) {
      res <- sapply(1:nRow, function(s) {
        analogues <- (1:nSamp)[order(d.mat[, s], decreasing = FALSE)][1:ana]
        nDWmean <- colMeans(envs[analogues, , drop = FALSE])
        c(nDWmean = nDWmean)
      })
      t(res)
    })
  }
  if (missing(core.spp)) run <- "jack"
  nSamp <- nrow(training.spp)
  if (run == "both" || run == "core") {
    spp <- rbind(training.spp, core.spp)
  } else {
    spp <- training.spp
  }
  dist.mat <- make.dist(spp, method = method)
  diag(dist.mat) <- Inf
  if (!run == "core") {
    jack.dist.mat <- dist.mat[1:nSamp, 1:nSamp]
    jack <- ests(
      d.mat = jack.dist.mat, nRow = nrow(training.spp),
      nSamp = nSamp
    )
  }
  if (!run == "jack") {
    core.dist.mat <- dist.mat[1:nSamp, -(1:nSamp)]
    core <- ests(d.mat = core.dist.mat, nRow = nrow(core.spp), nSamp = nSamp)
  }
  if (run == "both") {
    list(jack = jack, core = core)
  } else if (run == "jack") {
    jack
  } else if (run == "core") {
    core
  }
}
