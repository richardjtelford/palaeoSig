#' Calculate the effective number of species in the fossil data
#' @description
#' Gives a measure of the species diversity in the fossil data.
#' @param spp data.frame of Species data
#' @details
#' Uses [rioja::Hill.N2()] from the rioja package. 
#' @returns 
#' Minimum, first quartile and median effective number of species
#' @references 
#' Hill, M. O. (1973) Diversity and evenness: a unifying notation and its 
#' consequences. \emph{Ecology} \bold{54}: 427--432. 
#' @author  Richard Telford
#' @note
#' If the effective number of species is small, WA based reconstructions are 
#' unlikely to be significant, and MAT based reconstructions should be tested 
#' instead.
#' @seealso [rioja::Hill.N2()]
#' @examples
#' require(rioja)
#' data(RLGH)
#' Hill.N2.core(RLGH$spec)
#' @keywords attribute
#' @importFrom rioja Hill.N2
#' @importFrom stats quantile
#' @export

Hill.N2.core <- function(spp) {
  hn2 <- Hill.N2(spp, margin = 1)
  quantile(hn2, prob = c(0, .25, .5))
}
