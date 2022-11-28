#' centipede_plot
#' @description Plot of species WA optima and tolerance
#' @param x A tolerance weighted weighted-average model from
#' \code{\link[rioja]{WA}}
#' @param spp data.frame of species data used to train the WA model
#' @param minN2 numeric giving minimum N2 for inclusion in plot
#' @param mult numeric multiplier for the tolerances
#' @details Extracts and sorts \code{\link[rioja]{WA}} optima and tolerances and
#'  generates a ggplot.
#'  Tends only to work well when there are a reasonable number of taxa,
#'  otherwise it is difficult to read the names on the axis.
#'  Rare taxa can be excluded with the \code{minN2} argument.
#'  The \code{tol.cut} argument in \code{\link[rioja]{WA}} may need to be set to
#'  prevent very small tolerances in rare taxa.
#' This function is very similar to the \code{\link[analogue]{caterpillar}}
#'  plot, but produces a ggplot
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' library(rioja)
#' data(SWAP)
#' mod <- WA(SWAP$spec, SWAP$pH, tolDW = TRUE)
#' coef(mod)
#' centipede_plot(mod, spp = SWAP$spec, minN2 = 20)
#' @importFrom ggplot2 ggplot geom_point geom_errorbar aes coord_flip
#' @importFrom dplyr filter inner_join mutate
#' @importFrom tibble as_tibble enframe
#' @importFrom magrittr %>%
#' @importFrom forcats fct_reorder
#' @importFrom assertr verify has_all_names
#' @importFrom rioja Hill.N2
#' @importFrom rlang .data
#' @export

centipede_plot <- function(x, spp, minN2 = 1, mult = 1) {
  # check WA object
  stopifnot(inherits(x, "WA"))

  # calculate N2
  n2 <- Hill.N2(spp) %>%
    enframe(name = "Taxon", value = "n2")

  # extract optima & tolerance
  opt_tol <- coef(x) %>%
    as_tibble(rownames = "Taxon") %>%
    verify(has_all_names("Optima", "Tolerances")) %>%
    inner_join(n2, by = "Taxon") %>%
    filter(.data$n2 >= minN2) %>%
    mutate(
      Taxon = factor(.data$Taxon),
      Taxon = fct_reorder(.data$Taxon, .data$Optima),
      ymin = .data$Optima - .data$Tolerances * mult,
      ymax = .data$Optima + .data$Tolerances * mult
    )

  # plot
  g <- ggplot(opt_tol, aes(
    x = .data$Taxon, y = .data$Optima,
    ymin = .data$ymin, ymax = .data$ymax
  )) +
    geom_errorbar() +
    geom_point() +
    coord_flip()

  return(g)
}
