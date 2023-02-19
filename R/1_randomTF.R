#' @name randomTF
#' @title Proportion of variance in the fossil data explained by an
#' environmental reconstruction
#' @description Calculate the proportion of variance in the fossil data
#' explained by an environmental reconstruction with a constrained ordination.
#' This value is compared with a null distribution calculated as the proportion
#' of variance in the fossil data explained by reconstructions from transfer
#' functions trained on random data.

#' @param spp Data frame of modern training set species data, transformed as
#' required for example with \code{sqrt}
#' @param env Data frame of training set environmental variables or vector with
#' single environmental variable
#' @param fos Data frame of fossil species data, with same species codes and
#' transformations as `spp`
#' @param n number of random training sets. More is better.
#' @param fun Transfer function method.
#' Additional arguments to \code{fun} can be passed with \code{...}
#' @param col Some transfer functions return more than one column of results,
#' for example with different \code{\link[rioja]{WAPLS}} components.
#' \code{col} selects which column of the reconstructions to use.
#' See the relevant transfer function method help file.
#' @param condition Optional data frame of reconstructions to partial out when
#' testing if multiple independent reconstructions are possible.
#' @param autosim Optional data frame of random values.
#' This is useful if the training set is spatially autocorrelated and the
#' supplied data frame contains autocorrelated random variables.
#' If \code{autosim} is missing, and \code{permute} is \code{FALSE}, the
#' transfer functions are trained on random variables drawn from a uniform
#' distribution.
#' @param ord Constrained ordination method to use. \code{\link[vegan]{rda}} is
#' the default, \code{\link[vegan]{cca}} should also work.
#' \code{\link[vegan]{capscale}} won't work without modifications to the code
#' (or a wrapper).
#' @param permute logical value. Generate random environmental variables by
#' permuting existing variable. Only possible if there is only one environmental
#' variable and \code{autosim} is missing.
#' @param models list of models made by \code{randomTF} with argument
#' \code{make_models = TRUE}
#' @param make_models logical, should a list of transfer functions trained on
#' random data be returned
#' @param \dots Other arguments to the transfer function. For example to change
#' the distance metric in \code{\link[rioja]{MAT}}.
#' Also extra arguments to plot.

#' @details The function calculates the proportion of variance in the fossil
#' data explained by the transfer function reconstruction.
#' This is compared with a null distribution of the proportion of variance
#' explained by reconstructions based on random environmental variables.
#' Reconstructions can be partialled out to test if multiple reconstructions are
#' statistically significant. If the environment is spatially autocorrelated, a
#' red-noise null should be used instead of the default white noise null.
#' Red noise environmental variables can be generated with the \pkg{gstat}
#' package.
#'
#' Any transfer function in the \pkg{rioja} package can be used. Other methods
#' (e.g. random forests) can be used by making a wrapper function.
#'
#' If reconstructions from several sites are to be tested using the same
#' training set it can be much faster to train the models on random
#' environmental data once and then use them repeatedly.
#' This can be done with \code{make_models = TRUE} and then running
#' \code{randomTF} again giving the resultant models to the \code{models}
#'  argument.
#'  \code{make_models} does not work with MAT.
#'
#' For some transfer function methods, including `WA`, the code can be made
#' somewhat faster by coercing the modern and fossil species data to matrices
#' (`spp <- as.matrix(spp)`), otherwise `WA` has to do this repeatedly.
#' With `MAT`, this should not be done as it might cause an error.
#'
#' @return
#' A list with components
#' \itemize{
#'    \item{PCA}{ The unconstrained ordination of the fossil data.}
#'    \item{preds}{ A list of the containing the reconstructions for each
#'    environmental variable.}
#'    \item{MAX}{ Proportion of the variance explained by the first axis of the
#'     unconstrained ordination.
#'     This is the maximum amount that a reconstruction of a single variable can
#'     explain.}
#'    \item{EX}{ The proportion of the variance in the fossil data explained by
#'     each reconstruction.}
#'    \item{sim.ex}{ The proportion of variance explained by each of the random
#'     environmental variables.}
#'    \item{sig}{ The p-value of each reconstruction.}
#' }
#' If \code{make_models = TRUE}, a list of transfer function models is returned.
#'
#'    \code{autoplot.palaeoSig} returns a \code{ggplot2} object

#' @references Telford, R. J. and Birks, H. J. B. (2011) A novel method for
#' assessing the statistical significance of quantitative reconstructions
#' inferred from biotic assemblages. \emph{Quaternary Science Reviews}
#' \bold{30}: 1272--1278.
#' \doi{10.1016/j.quascirev.2011.03.002}
#' @author Richard Telford \email{richard.telford@uib.no}
#' @note If there are only a few fossil levels, \code{\link{obs.cor}} might have
#' more power.
#' If there are few taxa, tests on \code{\link[rioja]{MAT}} reconstructions have
#' more statistical power than those based on \code{\link[rioja]{WA}}.
#' @seealso \code{\link{obs.cor}}, \code{\link[rioja]{WA}},
#' \code{\link[rioja]{MAT}}, \code{\link[rioja]{WAPLS}},
#' \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}}
#' @examples
#' require(rioja)
#' data(SWAP)
#' data(RLGH)
#' rlghr <- randomTF(
#'   spp = sqrt(SWAP$spec), env = data.frame(pH = SWAP$pH),
#'   fos = sqrt(RLGH$spec), n = 49, fun = WA, col = "WA.inv"
#' )
#' rlghr$sig
#' plot(rlghr, "pH")
#'
#' require("ggplot2")
#' autoplot(rlghr, "pH")

#' @keywords multivariate htest hplot

#' @importFrom rioja MAT
#' @importFrom purrr map map_dbl
#' @importFrom tibble lst
#' @importFrom vegan rda
#' @importFrom stats formula predict runif
#' @export

randomTF <- function(spp, env, fos, n = 99, fun, col,
                     condition = NULL, autosim, ord = rda,
                     permute = FALSE,
                     models,
                     make_models = FALSE, ...) {
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

  # check condition is data.frame and fos have same number of rows
  partial <- !is.null(condition)
  if (partial) {
    if (!is.data.frame(condition)) {
      stop("condition must be a data.frame of reconstructions to partial out")
    }
    if (!identical(nrow(fos), nrow(condition))) {
      stop("fos and condition must have the same number of rows")
    }
  }

  # make_models only?
  if (!missing(make_models)) {
    make_models <- isTRUE(make_models)
  } else {
    make_models <- FALSE
  }

  # MAT and make_models don't work well together
  if (identical(fun, MAT) && make_models) {
    stop("MAT and make_models don't work together because
         a shortcut is used to speed up MAT")
  }

  if (make_models && !missing(models)) {
    stop("If make_models is true, no not provide models")
  }

  if (!missing(models)) {
    if (!inherits(models, "model_list")) {
      stop("models must be a model_list made by running
           randomTF with make_models = TRUE")
    }
  }

  # if MAT, for speed, drop training set samples that are never analogues.
  if (identical(fun, MAT)) {
    mod1 <- predict(MAT(spp, env[[1]], ...), fos)
    analogues <- unique(as.vector(as.numeric(mod1$match.name)))
    spp <- spp[analogues, ]
    env <- env[analogues, , drop = FALSE]
    rownames(spp) <- seq_len(nrow(spp))
    if (!missing(autosim)) {
      autosim <- autosim[analogues, , drop = FALSE]
    }
  }

  # find inertia explained by first axis of unconstrained ordination
  if (!make_models) {
    # only if not in make_model mode
    if (!partial) {
      PC <- ord(fos)
    } else {
      conditions <- paste(names(condition), collapse = "+")
      form1 <- formula(paste("fos ~ 1 + Condition(", conditions, ")"))
      PC <- ord(form1, data = condition)
    }
    MAX <- PC$CA$eig[1] / PC$tot.chi

    # Find inertia explained by reconstructions
    obs <- lapply(env, FUN = function(ev) {
      mod <- fun(spp, ev, ...)
      Pred <- predict(mod, fos)
      if (is.list(Pred)) {
        col # force evaluation - otherwise lazy evaluation fails
        p <- Pred$fit[, col]
      } else {
        p <- Pred
      }
      if (!partial) {
        RDA <- ord(fos ~ p)
      } else {
        form <- formula(paste("fos ~ p + Condition(", conditions, ")"))
        RDA <- ord(form, data = condition)
      }

      list(
        EX = RDA$CCA$tot.chi / RDA$tot.chi,
        pred = p,
        EIG1 = RDA$CA$eig[1] / RDA$tot.chi,
        mod = Pred
      )
    })
  }

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

  # if MAT, can take shortcut as always same analogues chosen
  if (identical(fun, MAT)) {
    selected_analogues <- apply(obs[[1]]$mod$match.name, 2, as.numeric)
    p <- apply(selected_analogues, 1, function(n) {
      colMeans(rnd[n, ])
    })
    sim_ex <- apply(p, 1, function(pp) {
      if (!partial) {
        r <- ord(fos ~ pp)
      } else {
        form <- formula(paste("fos ~ pp + Condition(", conditions, ")"))
        r <- ord(form, data = condition)
      }
      r$CCA$tot.chi / r$tot.chi
    })
  } else {
    if (missing(models)) {
      # pre-calculated models not provided
      models <- apply(rnd, 2, function(sim) {
        m <- fun(spp, sim, ...)
        m
      })
    }
    if (make_models) {
      class(models) <- "model_list"
      return(models)
    }
    sim_ex <- sapply(models, function(m) {
      p <- predict(m, fos)
      if (is.list(p)) {
        p <- p$fit[, col]
      }
      if (!partial) {
        r <- ord(fos ~ p)
      } else {
        form <- formula(paste("fos ~ p + Condition(", conditions, ")"))
        r <- ord(form, data = condition)
      }
      r$CCA$tot.chi / r$tot.chi
    })
  }

  res <- lst(
    PCA = PC,
    preds = map(obs, "pred"),
    MAX = MAX,
    EX = map_dbl(obs, "EX"),
    eig1 = map_dbl(obs, "EIG1"),
    sim.ex = sim_ex,
    sig = map_dbl(.data$EX, function(e) mean(e <= c(e, sim_ex)))
  )
  class(res) <- "palaeoSig"
  return(res)
}
