#' @title Weighted correlation between weighted averaging optima and constrained
#' ordination species scores.
#' @description  obs.cor calculates the weighted correlation between the species
#' weighted average optima and the axis one species scores of an ordination
#' constrained by the WA reconstruction.
#' @param spp Data frame of modern training set species data, transformed if
#' required, for example with \code{sqrt}
#' @param env Vector of a single environmental variable
#' @param fos Data frame of fossil species data. Species codes and
#' transformations should match those in spp.
#' @param ord Constrained ordination method to use.
#' \code{\link[vegan]{rda}} is the default,
#' \code{\link[vegan]{cca}} should also work.
#' \code{\link[vegan]{capscale}} won't work without modifications to the code
#' (or a wrapper).
#' @param n number of random training sets. More is better.
#' @param min.occur Minimum number of occurrences of species in the species and
#' fossil data.
#' @param autosim Optional data frame of random values. This is useful if the
#' training set is spatially autocorrelated and the supplied data frame contains
#' autocorrelated random variables.
#' If \code{autosim} is missing, and \code{permute} is \code{FALSE}, the
#' transfer functions are trained on random variables drawn from a uniform
#' distribution.
#' @param permute logical value. Generate random environmental variables by
#' permuting existing variable.
#' Only possible if there is only one environmental variable and \code{autosim}
#' is missing.


#' @details  Obs.cor calculates the (weighted) correlation between the species
#' WA optima in the calibration set and their ordination axis one scores in the
#' fossil data. Seven different weights for the species are implemented.
#' \itemize{
#'   \item{"abun.fos" - weight by the mean abundance in the fossil data.}
#'    \item{"abun.calib" - weight by the mean abundance in the calibration data}
#'    \item{"abun.joint" - weight by the product of the mean abundance in the
#'    fossil and calibration data}
#'    \item{"n2.fos" - weight by the effective number of occurrences
#'    (Hill's N2) of each species in the fossil data}
#'    \item{"n2.calib" - weight by the effective number of occurrences
#'    (Hill's N2) of each species in the calibration data}
#'    \item{"n2.joint" - weight by the product of n2.calib and n2.fos}
#'    \item{"unweighted" - all species receive same weight.
#'    This is unlikely to be the best option but is included for completeness.}
#'    }
#'    It is unclear which of these weights is likely to be best:
#'    research is in progress.
#'    A square root transformation of the species data is often useful. n = 99
#'     is too small in practice to give a smooth histogram of the null model.
#'     n = 999 is better.

#' @return
#' obs.cor returns an obscor object, which is a list
#' \itemize{
#'   \item{ob}{ Observed correlation. Data.frame with columns Optima, RDA1 and
#'    abun containing the species optima, ordination axis 1 scores, and
#'    abundance used to weight the species respectively and a vector containing
#'    the weighted and unweighted correlations between species optima and
#'    ordination scores.}
#'   \item{sim}{ Matrix with the correlation between species weighted average
#'   optima and ordination scores in the first column and the weighted
#'   correlation in the second column. Each row represents a different random
#'   environmental variable.}
#'   \item{sigs}{ p-value for the observed correlation between species weighted
#'   average optima and ordination scores for each of the weights.}
#'  }


#' @references Telford, R. J. and Birks, H. J. B. (2011) A novel method for
#' assessing the statistical significance of quantitative reconstructions
#' inferred from biotic assemblages. \emph{Quaternary Science Reviews}
#' \bold{30}: 1272--1278.
#' \doi{10.1016/j.quascirev.2011.03.002}
#' @author Richard Telford \email{richard.telford@uib.no}
#' @note The test of the weighted correlation between species optima and
#' ordination axis scores is more powerful, especially with a small number of
#' fossil observations, that the test of variance explained in
#' \code{\link{randomTF}} but is only applicable to WA and will have a large
#' type II error if there are few species.

#' @seealso \code{\link{randomTF}}, \code{\link[rioja]{WA}},
#' \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}}
#' @examples
#' require(rioja)
#' data(SWAP)
#' data(RLGH)
#' rlgh.obs <- obs.cor(
#'   spp = sqrt(SWAP$spec),
#'   env = SWAP$pH,
#'   fos = sqrt(RLGH$spec),
#'   n = 49 # low number for speed
#' )
#' rlgh.obs$sig
#' plot(rlgh.obs, which = 1)
#' plot(rlgh.obs, which = 2)
#'
#' require(ggplot2)
#' autoplot(rlgh.obs, which = 1)
#' autoplot(rlgh.obs, which = 2, variable_names = "pH")
#' @keywords multivariate htest hplot

#' @importFrom vegan rda scores
#' @importFrom rioja WA
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stats cov.wt predict coef
#' @export


# to do
# allow multiple environmental variables to be considered -
# need to update plot funs for which = 1

obs.cor <- function(spp, env, fos, ord = rda, n = 99, min.occur = 1,
                    autosim, permute = FALSE) {
  # Check env is data.frame or vector
  if (!is.data.frame(env) && !is.vector(env)) {
    stop("env must be a data.frame containing one or more environemental
         variables, or vector containing a single environemental variable")
  }

  # force data to be a data.frame
  if (!is.data.frame(env)) {
    env <- data.frame(env = env)
  }
  rownames(spp) <- seq_len(nrow(spp))

  # check env and spp have same number of rows
  if (!identical(nrow(spp), nrow(env))) {
    stop("spp and env must have same number of rows")
  }

  # permute and autosim don't play together
  if (isTRUE(permute) && !missing(autosim)) {
    stop("permute does not make sense if autosim is provided")
  }

  # check only one variable if permute is true
  if (isTRUE(permute) && length(env) > 1) {
    stop("permute is only possible with one environmental variable")
  }

  # check min.occur is >0
  if (!is.numeric(min.occur) && min.occur > 0) {
    stop("min.occur must be numeric and > 0")
  }

  # remove taxa rarer than min.occur
  spp <- spp[, colSums(spp > 0) >= min.occur]
  fos <- fos[, colSums(fos > 0) >= min.occur]

  spp <- spp[, order(colnames(spp))]
  fos <- fos[, order(colnames(fos))]

  shared_fos <- names(fos) %in% names(spp)
  shared_spp <- names(spp) %in% names(fos)

  if (sum(shared_fos) == 0) {
    stop("No taxa in common to fossil and modern data sets")
  }

  # fit models to observed data
  mod <- WA(spp, env)
  pred <- predict(mod, fos)$fit[, 1]
  RDA <- ord(fos ~ pred)
  optima <- coef(mod)
  sco <- scores(RDA, display = "spec", choice = 1)

  abundances <- tibble(
    abun.fos = colMeans(fos[, shared_fos]),
    abun.calib = colMeans(spp[, shared_spp]),
    abun.joint = .data$abun.fos * .data$abun.calib,
    n2.fos = Hill.N2(fos[, shared_fos], margin = 2),
    n2.calib = Hill.N2(spp[, shared_spp], margin = 2),
    n2.joint = .data$n2.fos * .data$n2.calib
  )

  optima <- optima[shared_spp, , drop = FALSE]
  sco <- sco[shared_fos, , drop = FALSE]
  x <- data.frame(optima, sco, abundances)
  wcs <- sapply(abundances, function(abun) {
    abs(cov.wt(x[, 1:2], wt = abun, cor = TRUE)$cor[1, 2])
  })

  res_obs <- list(x = x, res = c(wcs, unweighted = abs(cor(x[, 1:2])[1, 2])))

  # simulations using random data
  # make random environmental variables
  if (!missing(autosim)) {
    # check autosim has correct size
    if (nrow(autosim) != nrow(env)) {
      stop("autosim must have same number of rows as env")
    }
    rnd <- autosim
  } else if (isTRUE(permute)) {
    rnd <- replicate(n = n, sample(env[[1]]), simplify = TRUE)
  } else {
    rnd <- matrix(runif(nrow(spp) * n), ncol = n)
  }

  # iterate over random variables
  res_sim <- apply(rnd, 2, function(sim) {
    mod <- WA(spp, sim)
    pred <- predict(mod, fos)$fit[, 1]
    RDA <- ord(fos ~ pred)
    optima <- mod$coef
    sco <- scores(RDA, display = "spec", choice = 1)
    optima_scores <- data.frame(
      optima = optima[intersect(rownames(optima), rownames(sco)), ,
        drop = FALSE
      ],
      sco = sco[intersect(rownames(sco), rownames(optima)), , drop = FALSE]
    )

    wcs <- sapply(abundances, function(abun) {
      abs(cov.wt(optima_scores[, 1:2], wt = abun, cor = TRUE)$cor[1, 2])
    })

    c(wcs, unweighted = abs(cor(x[, 1:2])[1, 2]))
  })
  res_sim <- as.data.frame(t(res_sim))
  sigs <- mapply(
    function(ob, sim) {
      mean(ob <= c(ob, sim))
    },
    ob = res_obs$res,
    sim = res_sim
  )

  res <- list(ob = res_obs, sim = res_sim, sigs = sigs)
  class(res) <- "obscor"
  return(res)
}
