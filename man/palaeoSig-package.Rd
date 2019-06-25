\name{palaeoSig-package}
\alias{palaeoSig-package}
\alias{palaeoSig}
\docType{package}
\title{
Significance Tests of Quantitative Palaeoenvironmental Reconstructions 
}
\description{
Significance tests for quantitative palaeoenvironmental reconstructions derived from transfer functions.
Functions from the \pkg{autocorTF} package are now included in \pkg{palaeoSig}.

}
\details{
This package includes:

 significance tests for quantitative palaeoenvironmental reconstructions (\code{\link{randomTF}}, \code{\link{obs.cor}})  
 
 graphical methods to show autocorrelation in transfer functions (\code{\link{rne}})                                       
 
 null model test of transfer functions performance in a spatially autocorrelated environment - see vignette. 
 
Several functions have from \pkg{autocorTF} version 1.0 and \pkg{palaeoSig} version 1.0 have been rewritten or replaced with more flexible or user friendly functions. See \code{news(package="palaeoSig")} for details.

See also my blog at https://quantpalaeo.wordpress.com/

\tabular{ll}{
Package: \tab palaeoSig\cr
Type: \tab Package\cr
Version: \tab 1.1-2\cr
Date: \tab 2013-04-17\cr
License: \tab GPL - 3\cr
LazyLoad: \tab yes\cr
}
}                                      
\author{
Richard Telford \email{richard.telford@bio.uib.no}
}
\references{
Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. DOI: 10.1016/j.quascirev.2008.12.020  

Telford, R. J. and Birks, H. J. B. (2011) A novel method for assessing the statistical significance of quantitative reconstructions inferred from biotic assemblages. \emph{Quaternary Science Reviews} \bold{30}: 1272--1278. DOI: 10.1016/j.quascirev.2011.03.002
}

\keyword{ package }
\seealso{
\code{\link[rioja:rioja-package]{rioja}}, \code{\link{rda}}
}
