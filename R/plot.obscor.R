#' @describeIn obs.cor Plots for obscor object
#' @param x An obscor object.
#' @param xlab X-axis label if the default is unsatisfactory.
#' @param ylab Y-axis label if the default is unsatisfactory.
#' @param f Scale factor for the abundances, the maximum cex of points for the
#' which=1 plot.
#' @param which Which type of plot. which = 1 gives a plot of RDA scores against
#' species optima. which = 2 gives a histogram showing the null distribution of
#' correlations between RDA scores and species optima, together with the
#' observed correlation.
#' @param variable_names Name of environmental variable (only 1 currently) for
#' the label on the observed correlation with which = 2
#' @param abun Which species weighting required for plots. See details
#' @param p_val P value to draw a line vertical line at (with which=2)

#' @importFrom graphics plot hist abline box
#' @importFrom stats quantile
#' @method plot obscor
#' @export

plot.obscor <- function(x, xlab, ylab, f = 5, which = 1,
                        variable_names = "env",
                        abun = "abun.calib", p_val = 0.05, ...) {
  weightings <- c(
    "abun.fos", "abun.calib", "abun.joint", "n2.fos",
    "n2.calib", "n2.joint", "unweighted"
  )
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
    } else {
      a <- x$ob$x[[w]]
      a <- a / max(a) * f
    }
    plot(
      x = x$ob$x$Optima, y = x$ob$x$RDA1, cex = a,
      xlab = xlab, ylab = ylab, ...
    )
  } else if (which == 2) {
    if (missing(xlab)) {
      xlab <- ifelse(w == "unweighted", "Correlation", "Weighted correlation")
    }
    sim <- x$sim[[w]]
    ob <- x$ob$res[w]
    hist(sim,
      xlim = range(c(sim, ob)), xlab = xlab,
      col = "grey80", border = NA, ...
    )
    abline(v = ob, col = 1)
    abline(v = quantile(sim, prob = 1 - p_val), col = 2, lty = 3)

    text(ob, par()$usr[4] * 0.9, label = variable_names, srt = 90, pos = 2)
    box()
  } else {
    stop("which==what")
  }
}


#' @describeIn obs.cor Identify species on obs.cor plot
#' @param labels Labels for the points in identify. By default, the species
#' names from intersection of colnames(spp) and colnames(fos) are used.
#' @param \dots Other arguments to plot or identify
#' @importFrom graphics identify
#'
identify.obscor <- function(x, labels, ...) {
  if (missing(labels)) {
    labels <- rownames(x$ob$x)
  }
  identify(x$ob$x[, 1:2], labels = labels, ...)
}




#' @describeIn obs.cor autoplot for obscor object
#' @param top Proportion of the figure below the environmental name labels.
#' @param nbins integer giving number of bins for the histogram
#' @importFrom ggplot2 autoplot ggplot geom_point scale_size_area
#' @importFrom stats quantile
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @method autoplot obscor
#' @export

autoplot.obscor <- function(x, which = 1, variable_names = "env",
                            abun = "abun.calib", p_val = 0.05,
                            nbins = 20, top = 0.7, ...) {
  weightings <- c(
    "abun.fos", "abun.calib", "abun.joint", "n2.fos",
    "n2.calib", "n2.joint", "unweighted"
  )
  w <- pmatch(abun, weightings)
  if (is.na(w)) {
    stop("Unknown abundance weighting")
  }
  abun <- weightings[w]

  if (which == 1) {
    x$ob$x %>%
      mutate(unweighted = 1) %>%
      ggplot(aes(x = .data$Optima, y = .data$RDA1, size = .data[[abun]])) +
      geom_point(alpha = 0.3) +
      scale_size_area() +
      labs(x = "WA optima", y = "RDA scores", size = "Abundance")
  } else if (which == 2) {
    xlab <- ifelse(w == "unweighted", "Correlation", "Weighted correlation")

    x_fort <- fortify_palaeosig(
      sim = x$sim[, abun],
      variable_names = variable_names,
      p_val = p_val,
      nbins = nbins,
      top = top,
      EX = x$ob$res[abun]
    )

    autoplot_sig(x_fort, xlab = xlab, xmin = NA_real_)
  } else {
    stop("Unknown plot")
  }
}
