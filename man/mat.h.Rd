\name{mat.h}
\alias{mat.h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MAT with minimum geographic distance for analogues }
\description{
h-block cross-validation for MAT. 
}
\usage{
mat.h(y, x, noanalogues, geodist, thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Community data or distance matrix }
  \item{x}{ Environmental variable }
  \item{noanalogues}{ Number of analogues }
  \item{geodist}{ Geographical distances }
  \item{thresh}{ Threshold distance below which sites are excluded during cross validation }
}
\details{
Analogues closer than the geographic distance specified by thresh are not excluded, reducing the effect of spatial autocorrelation on the performance statistics. 
H-block cross-validation is also implemented in the \pkg{rioja} package. This function can be faster if used many times as it can take a pre-calculated distance matrix.
  
}
\value{
  Predictions for each site.
}
\references{ Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. DOI: 10.1016/j.quascirev.2008.12.020 }
\author{ Richard Telford \email{Richard.Telford@bio.uib.no} }

\examples{

require(fields)
data(arctic.env)
data(arctic.pollen)
arctic.dist <- rdist.earth(arctic.env[,c("Longitude", "Latitude")], miles=FALSE)


h0 <- mat.h(arctic.pollen, arctic.env$tjul, noanalogues = 10, arctic.dist, thresh=-1)
h100 <- mat.h(arctic.pollen, arctic.env$tjul, noanalogues = 10, arctic.dist, thresh=100)

cor(h0, arctic.env$tjul)
cor(h100, arctic.env$tjul)

sqrt(mean((h0 - arctic.env$tjul)^2))
sqrt(mean((h100 - arctic.env$tjul)^2))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
