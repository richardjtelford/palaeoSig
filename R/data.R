#' Atlantic core-top foram assemblages
#'
#' A dataset containing over 1000 foram assemblages from the Atlantic from
#' Kucera et al (2005) and the 50m SST for the warmest season. Rare taxa and
#' co-located assemblages are removed.
#'
#' @format A data frame with 1093 rows and 33 variables.
#' summ50 is 50m water temperature of the warmest season
#' @usage data(Atlantic)
#' @source \doi{10.1594/PANGAEA.227322}
#' \url{https://www.nodc.noaa.gov/cgi-bin/OC5/woa13/woa13.pl?parameter=t}
#' @keywords datasets
"Atlantic"




#' @title Storsandsvatnet
#' @description Storsandsvatnet is a lake in western Norway. From the sediments
#' a core was obtained, and 11 samples was submitted for radiocarbon dating.
#' The data contain the depths of the slides dated and the younger and older
#' calibrated ages for each slide.
#' @format
#'  A data.frame with 11 observations on the following 4 variables.
#'    \describe{
#'        \item{depthup}{The upper border of the dated slide}
#'        \item{depthdo}{The lower border of the dated slide}
#'        \item{cageup}{The younger calibrated age of the dated slide}
#'        \item{cagedo}{The older calibrated age of the dated slide}
#'        }
#' @details  The calibrated ages is obtained by calibration of the radiocarbon
#' dates.
#' The borders represent mean calibrated age +/- 1 SD of calibrated age.
#' @source The data are unpublished and provided by H. John B. Birks
#' <john.birks@bio.uib.no> and Sylvia M. Peglar
#' @usage data(STOR)
#' @references  Heegaard, E., Birks, HJB. & Telford, RJ. 2005. Relationships
#' between calibrated ages and depth in stratigraphical sequences: an estimation
#' procedure by mixed-effect regression. The Holocene 15: 612-618 \doi{10.1191/0959683605hl836rr}
#' @keywords datasets
"STOR"


#' @name arctic.pollen
#' @title Arctic Pollen and associated environmental data
#' @description  Arctic pollen percent data and associated environmental data
#' @format
#'  \describe{
#'    \item{arctic.pollen}{A data frame with 828 observations on the percentage
#'    of 39 pollen taxa}
#'    \item{arctic.env}{Environmental data for the pollen sites}
#'  }
#' @source Data extracted from North American Pollen Database and
#' New \emph{et al.} (2002) by Fréchette \emph{et al.} (2008).
#' Following Fréchette (Pers. Comm.), three duplicate sites have been deleted.
#' @references  Fréchette, B., de Vernal, A., Guiot, J., Wolfe, A. P., Miller,
#' G. H., Fredskild, B., Kerwin, M. W. and Richard, P. J. H. (2008)
#' Methodological basis for quantitative reconstruction of air temperature and
#' sunshine from pollen assemblages in Arctic Canada and Greenland.
#' \emph{Quaternary Science Reviews} \bold{27}, 1197--1216
#' \doi{10.1016/j.quascirev.2008.02.016}
#' @usage data(arctic.pollen)
#' data(arctic.env)
#' @keywords datasets
#' @rdname arctic.pollen
"arctic.env"
