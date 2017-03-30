\name{obs.cor}
\alias{obs.cor}
\alias{plot.obscor}
\alias{identify.obscor}

\title{
Weighted correlation between weighted averaging optima and constrained ordination species scores.
}
\description{
obs.cor calculates the weighted correlation between the species weighted average optima and the axis one species scores of an ordination constrained by the WA reconstruction. 
}
\usage{
obs.cor(spp, env, fos, ord = rda, n = 99, min.occur = 1)
\method{plot}{obscor}(x, xlab, ylab, f = 5, which = 1, label = "env", abun="abun.calib", p.val=0.95, ...)
\method{identify}{obscor}(x, labels, ...)

}


\arguments{
  \item{spp}{
Data frame of modern training set species data, transformed if required, for example with sqrt  
}
  \item{env}{
Vector of a single environmental variable
}
  \item{fos}{
Data frame of fossil species data. Species codes and tranformations should match those in spp.
}
  \item{ord}{
Constrained ordination. The default is rda. cca will work with the present code. capscale would require the code to be altered.
}
  \item{n}{
Number of random trials
}
  \item{min.occur}{
Minimum number of occurrances of species in the species and fossil data.
}

  \item{x}{
An obscor object.
}
  \item{xlab}{
X-axis label if the default is unsatisfactory.
}
  \item{ylab}{
Y-axis label if the default is unsatisfactory.
}
  \item{f}{
Scale factor for the abundances, the maximum cex of points for the which=1 plot.
}
  \item{which}{
Which type of plot. which = 1 gives a plot of RDA scores against species optima. which = 2 gives a histogram showing the null distribution of correlations between RDA scores and species optima, together with the observed correlation. 
}
  \item{label}{
Label for the observed correlation with which = 2
}
  \item{abun}{
Which species weighting required for plots. See details
}
  \item{p.val}{
P value to draw a line vertical line at (with which=2)  
}
  \item{labels}{
Labels for the points in identify. By default, the species names from intersection of colnames(spp) and colnames(fos) are used.  
}

  \item{\dots}{
Other arguments to plot or identify
}

}
\details{
Obs.cor calculates the weighted correlation between the species WA optima in the calibration set and their ordination axis one scores in the fossil data. Seven different weights for the species are implimented. 
  "abun.fos" - weight by the mean abundance in the fossil data.
  "abun.calib" - weight by the mean abundance in the calibration data
  "abun.joint" - weight by the product of the mean abundance in the fossil and calibration data
  "n2.fos" - weight by the effective number of occurences (Hill's N2) of each species in the fossil data 
  "n2.calib" - weight by the effective number of occurences (Hill's N2) of each species in the calibration data 
  "n2.joint" - weight by the product of n2.calib and n2.fos
  "unweighted" - all species recieve same weight. This is unlikely to be the best option but is included for completeness. 
  It is unclear which of these weights is likely to be best, research is in progress.
A square root transformation of the species data is often useful. n=99 is too small in practice to give a smooth histogram of the null model. n=999 is better.
}
\value{

obs.cor returns an obscor object, which is a list
 \item{ob}{Obsevered correlation. Data.frame with columns Optima, RDA1 and abun containing the species optima, ordination axis 1 scores, and abundance used to weight the species respectively and a vectorcontaining the weighted and unweighted correlations between species optima and ordination scores.}
 \item{sim}{Matrix with the correlation between species weighted average optima and ordination scores in the first column and the weighted correlation in the second column. Each row represents a different random environmental variable.
}
 \item{sigs}{p-value for the observed correlation between species weighted average optima and ordination scores for each of the weights.}


}

\references{
Telford, R. J. and Birks, H. J. B. (2011) A novel method for assessing the statistical significance of quantitative reconstructions inferred from biotic assemblages. \emph{Quaternary Science Reviews} \bold{30}: 1272--1278. DOI: 10.1016/j.quascirev.2011.03.002

}
\author{
Richard Telford \email{richard.telford@bio.uib.no}
}
\note{
The test of the weighted correlation between species optima and ordination axis scores is more powerful, especially with a small number of fossil observations, that the test of variance explained in  \code{\link{randomTF}} but is only applicable to WA and will have a large type II error if there are few species. 
}

\seealso{
\code{\link{randomTF}}, \code{\link{WA}}, \code{\link{MAT}}, \code{\link{WAPLS}}, \code{\link{rda}}, \code{\link{cca}}
}
\examples{
require(rioja)
data(SWAP)
data(RLGH)

rlgh.obs <- obs.cor(spp = sqrt(SWAP$spec), env = SWAP$pH, fos = sqrt(RLGH$spec))
rlgh.obs$sig
plot(rlgh.obs, which=1)
plot(rlgh.obs, which=2)

}

\keyword{multivariate}
\keyword{htest}
\keyword{hplot}