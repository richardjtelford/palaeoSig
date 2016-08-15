\name{coverage.plot}
\alias{coverage.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coverage of fossil taxa in modern calibration set }
\description{
  A simple diagnostic plot showing the coverage of fossil taxa in modern calibration set
}
\usage{
coverage.plot(mod,fos, rare=5, identify=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{ Modern calibration set species data }
  \item{fos}{ Fossil species data }
  \item{rare}{ Value of Hill's N2 below which species are highlighted}
  \item{identify}{ Whether to identify selected taxa.}
}
\details{
  Finds the maximum abundance of fossil taxa and plots this against the maximum abundance the taxa in the modern calibration set. Taxa with a Hill's N2 less than \code{rare} in the calibration set are highlighted in blue. Taxa absent from the calibration set are highlighed in red. If there are many taxa above the 1:1 line, or important fossil taxa have a low N2 in the calibration set, reconstructions should be interpreted with caution.
}

\value{
  An invisible \code{data.frame} with the modern and fossil maximum abundances and N2 in the calibration set. 
}

\author{ Richard Telford \email{Richard.Telford@bio.uib.no}  }
\examples{

require(rioja)
data(SWAP)
data(RLGH)
coverage.plot(mod=SWAP$spec, fos=RLGH$spec, identify=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
