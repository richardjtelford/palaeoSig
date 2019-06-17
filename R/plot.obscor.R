#' @importFrom graphics plot hist abline box
#' @importFrom stats quantile
#' @export

plot.obscor <- function(x, xlab, ylab, f = 5, which = 1, label = "env", 
                        abun = "abun.calib", p.val = 0.05, ...) {
    weightings <- c("abun.fos", "abun.calib", "abun.joint", "n2.fos",
                    "n2.calib", "n2.joint", "unweighted")
    w <- pmatch(abun, weightings)
    if (is.na(w)) {
        stop("Unknown abundance weighting")
    }
    w <- weightings[w]
    

    if (which == 1) {
      if (missing(xlab)) {
        xlab <- "WA optima"
      }
      if (missing(ylab)) {
        ylab <- "RDA scores"
      }
      if (w == "unweighted") {
        a <- rep(1, nrow(x$ob$x))
      }
      else{
        a <- x$ob$x[[w]]
        a <- a / max(a) * f
      }
    plot(x = x$ob$x$Optima, y = x$ob$x$RDA1, cex = a, 
         xlab = xlab, ylab = ylab, ...)
  } else if(which == 2){
    if (missing(xlab)) {
      xlab <- ifelse(w == "unweighted", "Correlation", "Weighted correlation")
    }
    sim <- x$sim[[w]]
    ob <- x$ob$res[w]
    hist(sim, xlim = range(c(sim, ob)), xlab = xlab, 
         col = "grey80", border = NA, ...)
    abline(v = ob, col = 1)
    abline(v = quantile(sim, prob = 1 - p.val), col =  2, lty = 3)
    
    text(ob, par()$usr[4] * 0.9, label = label, srt = 90, pos = 2)
    box()
  }else stop("which==what")
}


#' @importFrom ggplot2 autoplot ggplot geom_point scale_size_area
#' @importFrom stats quantile
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export

autoplot.obscor <- function(x, which = 1, label = "env", 
                        abun = "abun.calib", p.val = 0.05, ...) {
  
  weightings <- c("abun.fos", "abun.calib", "abun.joint", "n2.fos",
                  "n2.calib", "n2.joint", "unweighted")
  w <- pmatch(abun, weightings)
  if (is.na(w)) {
    stop("Unknown abundance weighting")
  }
  abun <- weightings[w]
  
  if(which == 1){
    
    x$ob$x %>% 
      mutate(unweighted = 1) %>% 
      ggplot(aes(x = .data$Optima, y = .data$RDA1, size = .data[[abun]])) +
      geom_point(alpha  =0.3) +
      scale_size_area() +
      labs(x = "WA optima", y = "RDA scores", size = "Abundance")

  } else if(which == 2){
      xlab <- ifelse(w == "unweighted", "Correlation", "Weighted correlation")

    sim <- x$sim[[w]]
    ob <- x$ob$res[w]

    autoplot_sig(sim_bin, lines_to_label, width,  xlab = xlab)
  } else stop("Unknown plot")
}

