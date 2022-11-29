\name{multi.mat}
\alias{multi.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MAT for multiple variables }
\description{
  MAT for many environmental variables simultaneously. More efficient than calculating them separately for each variable.
}
\usage{
multi.mat(training.spp, envs, core.spp, noanalogues = 10, method = "sq-chord",
    run = "both")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training.spp}{ Community data }
  \item{envs}{ Environmental variables - or simulations }
  \item{core.spp}{ Optional fossil data to make predictions for}
  \item{noanalogues}{ Number of analogues to use }
  \item{method}{ distance metric to use }
  \item{run}{ Return LOO predictions or predictions for fossil data  }
}
\value{
  Matrix of predictions
}
\references{ Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. \doi{10.1016/j.quascirev.2008.12.020}  }
\author{ Richard Telford \email{Richard.Telford@bio.uib.no}  }
\examples{

data(arctic.env)
data(arctic.pollen)

mMAT <- multi.mat(arctic.pollen, arctic.env[,9:67], noanalogues = 5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
