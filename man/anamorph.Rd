\name{anamorph}
\alias{anamorph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Anamorph }
\description{
  Creates functions that transform arbitrary distributions into a Gaussian distributions, and vice versa.
}
\usage{
anamorph(x, k, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data }
  \item{k}{ number of Hermite polynomials }
  \item{plot}{ plot }
}
\details{
  Increasing k can give a better fit.
}
\value{
  Returns two function in a list
  \item{xtog }{Function to transform arbitrary variable x into a Gaussian distribution}
  \item{gtox }{The back transformation}
}
\references{ Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} 3rd edition, Springer-Verlag, Berlin. }
\author{ Richard Telford \email{Richard.Telford@bio.uib.no}  }
\examples{
set.seed(42)
x <- c(rnorm(50, 0, 1), rnorm(50, 6, 1))
hist(x)
ana.fun <- anamorph(x, 30, plot = TRUE)
xg <- ana.fun$xtog(x)
qqnorm(xg)
qqline(xg)
all.equal(x, ana.fun$gtox(xg))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

