#' @name rne
#' @aliases rne
#' @aliases plot.RNE
#'
#' @title Random, neighbour, environment deletion analysis for transfer function
#'  models
#' @description Calculates effect of deleting sites from training set at random,
#' from a geographic neighbourhood, or from an environmental neighbourhood.
#' A simple graphical technique for gauging the effect of spatial
#' autocorrelation on the transfer function model.
#'
#' @param   y Community data, or distance object, or distance matrix
#' @param   env Environmental variable
#' @param   geodist Matrix of geographical distances between sites
#' @param   fun Transfer function
#' @param   neighbours  Neighbourhood radii
#' @param   subsets  Proportion of sites to retain in random deletion
#' @param   ... Arguments passed to fun

#' @details Finds the leave-one-out transfer function performance if sites are
#' deleted at random, from a
#' neighbourhood zone, or by deleting environmentally close sites.
#'
#' Prior to version 2.1, this function would repeat the random removal 10 times
#' to reduce variance in results.
#' This is no longer done as the variance is small for large training sets,
#' it took a long time, and treats one treatment of the data differently.

#' @return
#' Returns an RNE object, list with two components
#' \itemize{
#'  \item{random }{Performance with random deletion.}
#'  \item{neighbour }{Performance with deletion by neighbourhood,
#'  or environment}
#'  }
#' @references Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer
#' functions in spatially structured environments.
#' \emph{Quaternary Science Reviews} \bold{28}: 1309--1316.
#' \doi{10.1016/j.quascirev.2008.12.020}
#' @author Richard Telford \email{Richard.Telford@bio.uib.no}

#' @examples
#' require(rioja)
#' require(sf)
#' data(arctic.env)
#' data(arctic.pollen)
#'
#' # using just the first 100 sites so that code runs quickly (about 15 seconds for all 828 sites)
#'
#' # convert environmental data into an sf object
#' arctic.env <- st_as_sf(
#'   x = arctic.env,
#'   coords = c("Longitude", "Latitude"),
#'   crs = 4326
#' )
#'
#' # find great circle distances and remove units
#' arctic.dist <- st_distance(arctic.env[1:100, ]) |>
#'   units::set_units("km") |>
#'   units::set_units(NULL)
#'
#' # rne
#' arctic.rne <- rne(
#'   y = arctic.pollen[1:100, ], env = arctic.env$tjul[1:100],
#'   geodist = arctic.dist, fun = MAT, neighbours = c(0, 200),
#'   subsets = c(1, .5), k = 5
#' )
#'
#' plot(arctic.rne)

#' @keywords multivariate

#' @importFrom stats cor predict
#' @importFrom rioja crossval performance
#' @export

rne <- function(y, env, geodist, fun, neighbours,
                subsets = c(1, 0.75, 0.5, 0.25, 0.1), ...) {
  dots <- list(...)
  if (inherits(geodist, "dist")) {
    geodist <- as.matrix(geodist)
  }
  # ensure diagonal of geographic distances is zero
  diag(geodist) <- 0

  rne <- list()
  nr <- nrow(y)

  # fit model
  mod <- do.call(
    "fun",
    c(list(y = y, x = env), dots, lean = FALSE)
  )

  # deletion at random


  rne$random <- t(sapply(subsets, function(ss) {
    message(paste("random subset = ", ss))

    # make pseudo-distance matrix, 1 to keep, 0, to drop
    rdist <- lapply(seq_len(nr), function(n) {
      d <- rep(0, times = nr)
      retain <- sample(seq_len(nr)[-n], size = round((nr - 1) * (ss)))
      d[retain] <- 1
      d
    })
    rdist <- do.call(what = "rbind", args = rdist)

    # cross validate
    mod_random_cv <- crossval(mod,
      cv.method = "h-block",
      h.cutoff = 0.5,
      h.dist = rdist
    )
    mod_random_r2 <- performance(mod_random_cv)$crossval[, "R2"]

    # extract only required k for MAT models
    if (identical(fun, MAT)) {
      # pick max k
      k <- mod$k
      mod_random_r2 <- mod_random_r2[c(k, k * 2)]
    }

    c(prop = ss, r2 = mod_random_r2)
  }))
  message(rne$random)

  # deletion from geographic/environmental neighbourhood
  rne$neighbour <- lapply(neighbours, function(neighbour) {
    message(paste("neighbourhood = ", neighbour, "km"))
    en <- sapply(seq_len(nr), function(n) {
      sum(geodist[n, ] > neighbour)
    })
    effn <- (nr - mean(en)) / (nr - 1)

    # cross validate
    mod_space_cv <- crossval(mod,
      cv.method = "h-block",
      h.cutoff = neighbour,
      h.dist = geodist
    )
    mod_space_r2 <- performance(mod_space_cv)$crossval[, "R2"]

    # extract only required k for MAT models
    if (identical(fun, MAT)) {
      # pick max k
      k <- mod$k
      mod_space_r2 <- mod_space_r2[c(k, k * 2)]
    }


    # delete by environmental distance
    # make pseudo distance matrix  - 1 to include in model, 0 to exclude

    edist <- as.matrix(dist(env))
    edist <- lapply(seq_len(nr), function(n) {
      e <- edist[n, ]
      # set actual to -1 so never matched
      e[n] <- -1
      as.numeric(rank(-e, ties.method = "random") <= (en[n]))
    })
    edist <- do.call(what = "rbind", args = edist)

    # cross validate
    mod_environment_cv <- crossval(mod,
      cv.method = "h-block",
      h.cutoff = 0.5,
      h.dist = edist
    )
    mod_environment_r2 <- performance(mod_environment_cv)$crossval[, "R2"]

    # extract only required k for MAT models
    if (identical(fun, MAT)) {
      # pick max k
      k <- mod$k
      mod_environment_r2 <- mod_environment_r2[c(k, k * 2)]
    }

    list(
      neighbour = neighbour,
      effn = effn,
      hb.r2 = mod_space_r2,
      eb.r2 = mod_environment_r2
    )
  })

  class(rne) <- "RNE"
  return(rne)
}
