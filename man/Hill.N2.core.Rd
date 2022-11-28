\name{Hill.N2.core}
\alias{Hill.N2.core}

\title{
Calculate the effective number of species in the fossil data
}
\description{
Gives a measure of the species diversity in the fossil data.
}
\usage{
Hill.N2.core(spp)
}
\arguments{
  \item{spp}{
Species data
}
}
\details{
Uses \code{\link{Hill.N2}} from the rioja package 

}
\value{
Minimum, first quartile and median effective number of species
}
\references{
Hill, M. O. (1973) Diversity and evenness: a unifying notation and its consequences. \emph{Ecology} \bold{54}: 427--432. 

}
\author{
Richard Telford
}
\note{
If the effective number of species is small, WA based reconstructions are unlikely to be significant, and MAT based reconstructions should be tested instead.
}


\seealso{
\code{\link{Hill.N2}}
}
\examples{
require(rioja)
data(RLGH)

Hill.N2.core(RLGH$spec)
}
\keyword{ attribute }

