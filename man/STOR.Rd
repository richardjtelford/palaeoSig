\name{STOR}
\alias{STOR}
\docType{data}
\title{ Storsandsvatnet }
\description{
  Storsandsvatnet is a lake in western Norway. From the sediments a core was obtained, and 11 samples was submitted for radiocarbon dating. The data contain the depths of the slides dated and the younger and older calibrated ages for each slide.
}
\usage{data(STOR)}
\format{
  A data frame with 11 observations on the following 4 variables.
  \describe{
    \item{depthup}{The upper border of the dated slide}
    \item{depthdo}{The lower border of the dated slide}
    \item{cageup}{The younger calibrated age of the dated slide}
    \item{cagedo}{The older calibrated age of the dated slide}
  }
}
\details{
  The calibrated ages is obtained by calibration of the radiocarbon dates. The borders represent mean calibrated age +/- 1 SD of calibrated age.
}
\source{
  The data is unpublished and provided by H. John B. Birks <john.birks@bio.uib.no> and Sylvia M. Peglar
}
\references{
  Heegaard, E., Birks, HJB. & Telford, RJ. 2005. Relationships between calibrated ages and depth in stratigraphical sequences: an estimation procedure by mixed-effect regression. The Holocene 15: 612-618
}
\examples{
data(STOR)

}
\keyword{datasets}
