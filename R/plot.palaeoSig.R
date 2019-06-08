#' @importFrom graphics hist lines text strwidth
#' @importFrom stats quantile
#' @importFrom TeachingDemos spread.labs
#' 
plot.palaeoSig <- function(x, vnames, top = 0.7, adj = c(0,0.5), p.val = 0.05, ...){
  if (missing(vnames)) {
    vnames <- names(x$EX)
  }
  with(x,{
    hist(sim.ex, breaks = seq(min(sim.ex), max(sim.ex), length = 20), 
         xlim = c(0, MAX[1] * 1.1), 
         main = "", xlab = "Proportion variance explained", 
         col= "grey80", border = NA, ...)
    tops <- par()$usr[4] * top
    sapply(EX, function(z)
      lines(rep(z, 2), c(0, tops)))
    #abline(v=MAX, col=1, lwd=2, lty=3)
    lines(rep(MAX, 2), c(0, tops), col = 1, lwd = 2, lty = 3)
    lines(rep(quantile(sim.ex, probs = 1 - p.val), 2), c(0, tops), 
          col = 2, lwd = 1, lty = 3)
    #abline(v=EX, col=1)     
    putEX <- spread.labs(EX, 1.2 * strwidth('A', cex = .8))

    text(putEX, par()$usr[4] * .71, label = vnames, 
         srt = 90, adj = adj, cex = .8)
  })
}

#' @importFrom   ggplot2 autoplot ggplot aes geom_col geom_linerange geom_text scale_colour_identity scale_linetype_identity labs
#' @importFrom tibble tibble
#' @importFrom rlang .data

autoplot.palaeoSig <- function(x, variable_names, nbins = 20, top = 0.7, p.val = 0.05){
  if (missing(variable_names)) {
    variable_names <- names(x$EX)
  }
  
  breaks <- seq(min(x$sim.ex), max(x$sim.ex), length = nbins + 1)
  id <- cut(x$sim.ex, breaks = breaks, include.lowest = TRUE)
  
  sim_bin <- tibble(  
    mid_point = (breaks[-length(breaks)] + breaks[-1]) / 2,
    n = as.vector(table(id))
  )

  width <- diff(sim_bin$mid_point[1:2])
  
  lines_to_add <- tibble(
    label = c("PC1", paste("p =", p.val), variable_names), 
    value = c(x$MAX, quantile(x$sim.ex, probs = 1 - p.val), x$EX), 
    max = max(sim_bin$n) * top, 
    linetype = c("dashed", "dotted", rep("solid", length(variable_names))),
    colour = c("black", "red", rep("black", length(variable_names)))
  )
  
  autoplot_sig(sim_bin, lines_to_add, xlab = "Proportion variance explained")
}    

autoplot_sig <- function(sim_bin, lines_to_add, xlab ){
  g <- ggplot(sim_bin, aes(x = .data$mid_point, y = .data$n)) +
    geom_col(fill = "grey70", width = width) +
    geom_linerange(data = lines_to_add, 
                  aes(x = .data$value, ymin = 0, ymax = .data$max, linetype = .data$linetype, colour = .data$colour), 
                  inherit.aes = FALSE) +
    geom_text(data = lines_to_add, 
                   aes(x = .data$value, y = .data$max, label = .data$label,
                       colour = .data$colour),
           # method = "top.qp",
            angle = 90, hjust = 0) +
    scale_colour_identity() +
    scale_linetype_identity() +
    xlim(0, NA) +
    labs(x = xlab, y = "Frequency")
  
  return(g)
}

