#' @describeIn randomTF Plot palaeoSig object
#' @param x Output from randomTF
#' @param variable_names Names of environmental variables. If missing, taken from \code{env} data.frame.
#' @param top Proportion of the figure below the environmental name labels.
#' @param adj Adjust the position that the environmental names are plotted at.
#' @param p_val P value to draw a line vertical line at (with which=2)  
#' 
#' @importFrom graphics hist lines text strwidth
#' @importFrom stats quantile
#' @importFrom TeachingDemos spread.labs
#' @method plot palaeoSig
#' @export
#' 
plot.palaeoSig <- function(x, variable_names, top = 0.7, 
                           adj = c(0, 0.5), p_val = 0.05, ...){
  if (missing(variable_names)) {
    variable_names <- names(x$EX)
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
    lines(rep(quantile(sim.ex, probs = 1 - p_val), 2), c(0, tops), 
          col = 2, lwd = 1, lty = 3)
    #abline(v=EX, col=1)     
    putEX <- spread.labs(EX, 1.2 * strwidth('A', cex = .8))

    text(putEX, par()$usr[4] * .71, label = variable_names, 
         srt = 90, adj = adj, cex = .8)
  })
}

#' @describeIn randomTF autoplot function for palaeoSig object
#' @param nbins integer giving number of bins for the histogram
#' @importFrom   ggplot2 autoplot ggplot aes geom_col geom_linerange geom_text scale_colour_identity scale_linetype_identity labs
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom ggrepel geom_text_repel
#' @method autoplot palaeoSig
#' @export


autoplot.palaeoSig <- function(x, variable_names, 
                               nbins = 20, top = 0.7, p_val = 0.05){
  if (missing(variable_names)) {
    variable_names <- names(x$EX)
  }
  
  x_fort <- fortify_palaeosig(
    sim = x$sim, 
    variable_names = variable_names, 
    p_val = p_val,
    nbins = nbins,
    top = top,
    PC1 = x$MAX,
    EX = x$EX
  )
  
  autoplot_sig(x_fort, xlab = "Proportion variance explained", xmin = 0)
}  

#' @importFrom tibble tibble lst
#' @importFrom rlang .data

fortify_palaeosig <- function(sim, variable_names, p_val, nbins, 
                              top, PC1 = NA, EX){  
  
  breaks <- seq(min(sim), max(sim), length = nbins + 1)
  id <- cut(sim, breaks = breaks, include.lowest = TRUE)
  
  sim_bin <- tibble(  
    mid_point = (breaks[-length(breaks)] + breaks[-1]) / 2,
    n = as.vector(table(id))
  )
  
  width <- diff(sim_bin$mid_point[1:2])
  
  lines_to_add <- tibble(
    label = c("PC1", paste("p =", p_val), variable_names), 
    value = c(PC1, quantile(sim, probs = 1 - p_val), EX), 
    max = max(sim_bin$n) * top, 
    linetype = c("dashed", "dotted", rep("solid", length(variable_names))),
    colour = c("black", "red", rep("black", length(variable_names)))
  ) %>% 
    filter(!is.na(.data$value))
  
  result <- lst(sim_bin, lines_to_add, width)
  return(result)
  
}    

#' @importFrom ggplot2 xlim

autoplot_sig <- function(x, xlab, xmin){
  g <- ggplot(x$sim_bin, aes(x = .data$mid_point, y = .data$n)) +
    geom_col(fill = "grey70", width = x$width) +
    geom_linerange(data = x$lines_to_add, 
                   aes(x = .data$value, ymin = 0, ymax = .data$max, 
                       linetype = .data$linetype, colour = .data$colour), 
                   inherit.aes = FALSE) +
    geom_text_repel(data = x$lines_to_add, 
                    aes(x = .data$value, y = .data$max, label = .data$label,
                        colour = .data$colour),
                    angle = 90, hjust = .5, vjust = 0, direction = "x") +
    scale_colour_identity() +
    scale_linetype_identity() +
    xlim(xmin, NA_real_) +
    labs(x = xlab, y = "Frequency")
  
  return(g)
}

