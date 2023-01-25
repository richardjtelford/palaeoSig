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
#' @param   nrep integer, number of times to delete sites at random

#' @details Finds the leave-one-out transfer function performance if sites are
#'  deleted at random (repeated 10 times to reduce variance in results), from a
#'   neighbourhood zone, or by deleting environmentally close sites.

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
#' require(fields)
#' require(rioja)
#' data(arctic.env)
#' data(arctic.pollen)
#'
#' # using just the first 20 sites so that code runs in an reasonable time
#' arctic.dist <- rdist.earth(
#'   x1 = arctic.env[1:20, c("Longitude", "Latitude")],
#'   miles = FALSE
#' )
#' arctic.rne <- rne(
#'   y = arctic.pollen[1:20, ], env = arctic.env$tjul[1:20],
#'   geodist = arctic.dist, fun = MAT, neighbours = c(0, 200),
#'   subsets = c(1, .5), nrep = 2, k = 5
#' )
#'
#' plot(arctic.rne)

#' @keywords multivariate


#' @importFrom stats cor predict
#' @export

rne <- function(y, env, geodist, fun, neighbours,
                subsets = c(1, 0.75, 0.5, 0.25, 0.1), nrep = 10, ...) {
  dots <- list(...)
  if (inherits(geodist, "dist")) {
    geodist <- as.matrix(geodist)
  }
  rne <- list()
  nr <- nrow(y)
  
  # deletion at random
  rne$random <- t(sapply(subsets, function(ss) {
    print(paste("random subset = ", ss))
    r2 <- replicate(nrep, {
      est <- sapply(seq_len(nr), function(n) {
        retain <- sample(seq_len(nr - 1), size = round((nr - 1) * (ss)))
        y2 <- y[-n, ][retain, ]
        keepcols <- colSums(y2) != 0
        mod <- do.call(
          "fun",
          c(list(y = y2[, keepcols], x = env[-n][retain]), dots)
        )
        predict(mod, y[n, keepcols, drop = FALSE])$fit
      })
      if (is.null(dim(est))) {
        cor(est, env)^2
      } else {
        apply(est, 1, cor, env)^2
      }
    })
    if (is.null(dim(r2))) {
      r2 <- mean(r2) # average across replicates
    } else {
      r2 <- rowMeans(r2) # average across replicates
    }
    print(r2)
    c(prop = ss, r2 = r2)
  }))
  print(rne$random)
  
  # deletion from geographic/environmental neighbourhood
  rne$neighbour <- lapply(neighbours, function(neighbour) {
    print(paste("neighbourhood = ", neighbour, "km"))
    en <- sapply(seq_len(nr), function(n) {
      sum(geodist[n, ] >= neighbour)
    })
    effn <- (nr - mean(en)) / (nr - 1)
    hb <- sapply(seq_len(nr), function(n) {
      y1 <- y[-n, ]
      env1 <- env[-n]
      exneigh <- geodist[n, -n] >= neighbour
      y2 <- y1[exneigh, ]
      keepcols <- colSums(y2) > 0
      mod <- do.call(
        "fun",
        c(list(y = y2[, keepcols], x = env1[exneigh]), dots)
      )
      predict(mod, y[n, keepcols, drop = FALSE])$fit
    })
    if (is.null(dim(hb))) {
      hbr <- cor(hb, env)^2
    } else {
      hbr <- apply(hb, 1, cor, env)^2
    }
    
    # delete by environmental distance
    eb <- sapply(seq_len(nr), function(n) {
      y1 <- y[-n, ]
      env1 <- env[-n]
      neigh <- which(
        rank(-abs(env1 - env[n]), ties.method = "random") <= (en[n])
      )
      y2 <- y1[neigh, ]
      keepcols <- colSums(y2) != 0
      mod <- do.call(
        "fun",
        c(list(y = y2[, keepcols], x = env1[neigh]), dots)
      )
      predict(mod, y[n, keepcols, drop = FALSE])$fit
    })
    if (is.null(dim(eb))) {
      ebr <- cor(eb, env)^2
    } else {
      ebr <- apply(eb, 1, cor, env)^2
    }

    list(
      neighbour = neighbour,
      effn = effn,
      hb.r2 = hbr,
      eb.r2 = ebr
    )
  })

  class(rne) <- "RNE"
  return(rne)
}
