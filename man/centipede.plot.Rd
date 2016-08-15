\name{centipede.plot}
\alias{centipede.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plots species optima and tolerances }
\description{
  Plots species optima and tolerances as calculated by a \code{\link{WA}} model.
}
\usage{
centipede.plot(x,keep=TRUE, xlab="", xlim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{\link{WA}} model. }
  \item{keep}{ A logical vector of which taxa to include in the plot. Defaults to all taxa. }
  \item{xlab}{ The xlab of the plot}
  \item{xlim}{ xlims for the plot. Defaults to a sensible range.}
  \item{\dots}{ Extra arguments to plot, for example of the \code{cex}, \code{col}, and \code{pch} of the points.}
}
\details{
  Extracts, sorts and plots \code{\link{WA}} optima and tolerances. Tends only to work well when there are a reasonable number of taxa, otherwise it is difficult to read the names on the axis. If rare taxa are not excluded via the \code{keep} argument, the \code{tol.cut} argument in \code{\link{WA}} may need to be set to prevent very small tolerances of these taxa. 
}

\value{
  None
}

\author{ Richard Telford \email{Richard.Telford@bio.uib.no}  }
\examples{

require(rioja)
data(SWAP)
mod<-WA(SWAP$spec, SWAP$pH, tolDW=TRUE)
keep<-colSums(SWAP$spec)>40#only common taxa
opar<-par(no.readonly = TRUE)
par(mar=c(5.1,6.1,2.1,1.1))#extra space on left side of plot for names
centipede.plot(x=mod, keep=keep, pch=20, cex.axis=.7)
par(opar)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
