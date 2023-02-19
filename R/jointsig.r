#' @name jointsig
#' @aliases jointsig
#' @aliases plot.js
#' @title Test if two variables jointly control changes in fossil data
#' @description Generates synthetic variables with different proportion of two
#' environmental variables, and tests how much variance in the fossil data
#' reconstructions of these synthetic variables explain.

#' @param spp Data frame of modern training set species data, transformed as
#' required, for example with sqrt
#' @param fos Data frame of fossil species data, with same species codes and
#' transformations as spp
#' @param var1 Training set environmental variable 1
#' @param var2 Training set environmental variable 2
#' @param n number of random training sets used to generate the null model
#' @param method Which significance test to use.
#' Current option are randomTF and obs.cor.
#' The latter may give strange results - use with caution.
#' @param r How many synthetic variables to make. More is better but slower

#' @param \dots Other arguments to the significance test
#' (some of these are required)

#' @details With `method="randomTF"`, the function calculates the proportion of
#' variance in the fossil data explained by transfer function reconstructions of
#' synthetic variables.
#' The synthetic variables are composed of two environmental variables, weighted
#' between -1 and +1, so to represent a circle.
#' This is compared with a null distribution of the proportion of variance
#' explained by reconstructions based on random environmental variables.
#' Any transfer function in the rioja library can be used.
#' With method="obs.cor", the aim is the same, but the function reports the
#' correlation between the species weighted average optima on the synthetic
#' variables and the species first axis scores.
#' This option has some pathological behaviour and should probably be avoided.
#'
#' @return
#'  A list with components
#' \itemize{
#'  \item{PCA}{ The unconstrained ordination of the fossil data.}
#'  \item{preds}{ A list of the containing the reconstructions for each
#'  environmental variable.}
#'  \item{MAX}{ Proportion of the variance explained by the first axis of the
#'  unconstrained ordination. This is the maximum amount that a reconstruction
#'  of a single variable can explain.}
#'  \item{EX}{ The proportion of the variance in the fossil data explained by
#'  each reconstruction.}
#'  \item{sim.ex}{ The proportion of variance explained by each of the random
#'  environmental variables.}
#'  \item{sig}{ The p-value of each reconstruction.}
#' }
#'
#' @references Unpublished method - use with caution. Can give spurious results
#' with weighted averaging.
#' @author Richard Telford \email{richard.telford@bio.uib.no}
#' @seealso \code{\link{randomTF}},\code{\link{obs.cor}}
#' @examples
#' require(rioja)
#' data(SWAP)
#' data(RLGH)
#'
#' rlgh.js <- jointsig(
#'   spp = sqrt(SWAP$spec),
#'   fos = sqrt(RLGH$spec),
#'   var1 = SWAP$pH,
#'   var2 = sample(SWAP$pH),
#'   method = "randomTF",
#'   n = 49, r = 32, fun = WA, col = 1
#' )
#' # nonsense second variable
#'
#' plot(rlgh.js, c("acid", "alkaline"), c("down", "up"))
#' @keywords multivariate htest hplot



#' @export

jointsig <- function(spp, fos, var1, var2,
                     method = "randomTF", n = 99, r = 32, ...) {
  if (r %% 4 != 0) {
    warning("r not divisible by 4")
  }
  v1 <- sin(seq(0, 2 * pi, length = r + 1))[-(r + 1)]
  v2 <- cos(seq(0, 2 * pi, length = r + 1))[-(r + 1)]
  syn_env <- mapply(function(v1, v2) {
    scale(var1) * v1 + scale(var2) * v2
  }, v1 = v1, v2 = v2)
  syn_env <- as.data.frame(syn_env)
  if (method == "randomTF") {
    res <- randomTF(spp, fos, ..., env = syn_env, n = n)
  } else if (method == "obs.cor") {
    warning("obscor can pathological results in jointsig")
    res <- list()

    res$EX <- sapply(syn_env, function(e) {
      obs.cor(spp, fos, ..., env = e, n = 0)$ob$res$wc
    }) ## edit to match new obsc
    res$sim.ex <- obs.cor(spp, fos, ..., env = syn_env[, 1], n = n)$sim[, 2]
  }
  res$v1 <- v1
  res$v2 <- v2
  class(res) <- "js"
  return(res)
}
