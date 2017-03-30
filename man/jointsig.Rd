\name{jointsig}
\alias{jointsig}
\alias{plot.js}
\title{
Test if two variables jointly control changes in fossil data
}
\description{
Generates synthetic variables with different proportion of two environmental variables, and tests how much variance in the fossil data reconstructions of these synthetic variables explain.
}

\usage{
jointsig(spp, fos, var1, var2, method = "randomTF", n = 99, r = 32, ...)
\method{plot}{js}(x, names.v1, names.v2, ...)
}

\arguments{
  \item{spp}{
Data frame of modern training set species data, transformed as required for example with sqrt
}
  \item{fos}{
Data frame of fossil species data, with same species codes and transformations as spp
}
  \item{var1}{
Training set environmental variable 1.
}
  \item{var2}{
Training set environmental variable 2.
}
  \item{n}{
number of random training sets used to generate the null model
}
  \item{method}{
Which significance test to use. Current option are randomTF and obs.cor. The latter may give strange results - use with caution.
}
  \item{r}{
How many synthetic variables to make. More is better but slower.
}
  \item{x}{
Output from jointsig
}
  \item{names.v1}{
Vector length 2 with names of the end members of the first environmental variable, e.g., c("cold", "warm") for temerature.
}
  \item{names.v2}{
Ditto for the second variable.
}
  \item{\dots}{
Other arguments to the significance test (some of these are required) or plot.
}

}
\details{
With method="randomTF", the function calculates the proportion of variance in the fossil data explained by transfer function reconstructions of synthetic variables. The synthetic variables are composed of two environmental variables, weighted betweeen -1 and +1, so to represent a circle. This is compared with a null distribution of the proportion of variance explained by reconstructions based on random environmental variables. Any transfer function in the rioja library can be used. 
With method="obs.cor", the aim is the same, but the function reports the correlation between the species weighted average optima on the synthetic variables and the species first axis scores. This option has some pathological behaviour and should probably be avoided.
}
\value{
A list with components
 \item{PCA}{The unconstrained ordination of the fossil data.}
\item{preds}{A list of the containing the reconstructions for each environmental variable.}
\item{MAX}{Proportion of the variance explained by the first axis of the unconstrained ordination. This is the maximum amount that a reconstruction of a single variable can explain.}
\item{EX}{The proportion of the variance in the fossil data explained by each reconstruction.}
\item{sim.ex}{The proportion of variance explained by each of the random environmental variables.}
\item{sig}{The p-value of each reconstruction.}



}
\references{
Unpublished method - use with caution. Can give spurious results with weighted averaging.
}
\author{
Richard Telford \email{richard.telford@bio.uib.no}
}


\seealso{
\code{\link{randomTF}},\code{\link{obs.cor}} 
}
\examples{
require(rioja)
data(SWAP)
data(RLGH)

rlgh.js <- jointsig(spp = sqrt(SWAP$spec), fos = sqrt(RLGH$spec), var1 = SWAP$pH, 
    var2 = sample(SWAP$pH), method = "randomTF", n = 99, r = 32, fun = WA, col = 1)
    #nonsense second variable

plot(rlgh.js, c("acid", "alkaline"), c("down", "up"))
}

\keyword{multivariate}
\keyword{htest}
\keyword{hplot}

