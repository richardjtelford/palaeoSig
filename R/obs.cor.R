#obs.cor <-
#function(spp,env,fos, ord=rda){
#  mod<-WA(spp,env)
#  pred<-predict(mod, fos)$fit[,1]
#  RDA<-ord(fos~pred)
#  optima<-mod$coef
#  sco<-scores(RDA, display="spec", choice=1)
#  abun<-t(t(colMeans(spp)))
#  
#  optima<-optima[intersect(rownames(optima),rownames(sco)),, drop=FALSE]
#  abun<-abun[intersect(rownames(abun),rownames(sco)),, drop=FALSE]
#  sco<-sco[intersect(rownames(sco),rownames(optima)),, drop=FALSE]
#  x<-data.frame(optima, sco, abun=abun)
#  res<-list(x=x,res=list(wc=abs(cov.wt(x[,1:2], wt=sqrt(x$abun), cor=TRUE)$cor[1,2]),cc=abs(cor(x[,1:2])[1,2])))
#  class(res)<-"obscor"
#  return(res)
#}
#code alternative weighting strategies
#fossil abundance
#mean training abundance
#product of fossil and training abundance
#N2 of fossil/training/joint

#' @importFrom vegan rda scores
#' @importFrom rioja WA
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stats cov.wt predict coef
#' @export

obs.cor <- function (spp, env, fos, ord = rda, n = 99, min.occur = 1){
  spp <- spp[, colSums(spp > 0) >= min.occur]
  fos <- fos[, colSums(fos > 0) >= min.occur]
  
  spp <- spp[, order(colnames(spp))]
  fos <- fos[, order(colnames(fos))]
  
  shared_fos <- names(fos) %in% names(spp)
  shared_spp <- names(spp) %in% names(fos)
  
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
    wcs <- sapply(abundances, function(abun){
      abs(cov.wt(x[, 1:2], wt = abun, cor = TRUE)$cor[1, 2])
      })
    
    res.obs <- list(x = x, res = c(wcs, unweighted = abs(cor(x[, 1:2])[1, 2])))

    res.sim <- replicate(n, {
        mod <- WA(spp, runif(length(env)))
        pred <- predict(mod, fos)$fit[, 1]
        RDA <- ord(fos ~ pred)
        optima <- mod$coef
        sco <- scores(RDA, display = "spec", choice = 1)
        optima <- optima[intersect(rownames(optima), rownames(sco)), , drop = FALSE]
        sco <- sco[intersect(rownames(sco), rownames(optima)), , drop = FALSE]
        x <- data.frame(optima, sco)
        wcs <- sapply(abundances,function(abun)abs(cov.wt(x[, 1:2], wt = abun, cor = TRUE)$cor[1, 2]))

        c(wcs, unweighted = abs(cor(x[, 1:2])[1, 2]))
    })
  res.sim <- as.data.frame(t(res.sim))
  sigs <- mapply(
    function(ob, sim){mean(ob <= c(ob, sim))}, 
    ob = res.obs$res, 
    sim = res.sim)
    
  res <- list(ob = res.obs, sim = res.sim, sigs = sigs)
  class(res) <- "obscor"
  return(res)
}
