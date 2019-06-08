
#' @importFrom rioja MAT
#' @importFrom purrr map map_dbl
#' @importFrom tibble lst
#' @export

randomTF <- function(spp, env, fos, n = 99, fun, col,
                     condition, autosim, ord = rda, ...){
  #reconstruct random variables
  if (!is.list(env)){
    env <- list(env = env)
  }
  rownames(spp) <- 1:nrow(spp)
  
  #if MAT, for speed, drop training set samples that are never analogues.
  if (identical(fun, MAT)) {
    mod1 <- predict(MAT(spp, env[[1]], ...), fos)
    analogues <- unique(as.vector(as.numeric(mod1$match.name)))
    spp <- spp[analogues, ]
    env <- env[analogues, , drop = FALSE]
    rownames(spp) <- 1:nrow(spp)
  }
  
  partial <- !missing(condition)
  
  #find inertia explained by first axis of unconstrained ordination
  if (!partial) {
    PC <- ord(fos)
  } else{
    conditions <- paste(names(condition), collapse = "+")
    form1 <- formula(paste("fos ~ 1 + Condition(", conditions, ")"))
    PC <- ord(form1, data = condition)
  }
  MAX <- PC$CA$eig[1] / PC$tot.chi
  
  # Find inertia explained by reconstructions
  obs <- lapply(env, function(ev) {
    Mod <- fun(spp, ev, ...)
    Pred <- predict(Mod, fos)
    if (is.list(Pred)) {
      p <- Pred$fit[, col]
    }
    else {
      p <- Pred
    }
    if (!partial) {
      RDA <- ord(fos ~ p)
    } else{
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

  # simulations using random data  
  if (missing(autosim)) {
    rnd <- matrix(runif(nrow(spp) * n), ncol = n)
  }
  else {
    rnd <- autosim
  }
  if (identical(fun, MAT)) {
    #if MAT, can take shortcut as always same analogues chosen
    selected_analogues <- apply(obs[[1]]$mod$match.name, 2, as.numeric)
    p <- apply(selected_analogues, 1, function(n){
        colMeans(rnd[n, ])})
    sim.ex <- apply(p, 1, function(pp) {
      if (!partial) {
        r <- ord(fos ~ pp)
      } else{
        form <- formula(paste("fos ~ pp + Condition(", conditions, ")"))
        r <- ord(form, data = condition)
      }
      r$CCA$tot.chi / r$tot.chi
    })
  }
  else{
    sim.ex <- apply(rnd, 2, function(sim) {
      m <- fun(spp, sim, ...)
      p <- predict(m, fos)
      if (is.list(p))
        p <- p$fit[, col]
      if (!partial) {
        r <- ord(fos ~ p)
      } else{
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
      sim.ex = sim.ex,
      sig = map_dbl(EX, function(E) mean(E <= c(E, sim.ex)))
    )
  class(res) <- "palaeoSig"
  return(res)
}

