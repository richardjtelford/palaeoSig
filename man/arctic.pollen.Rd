\name{arctic.pollen}
\encoding{latin1}
\alias{arctic.pollen}
\alias{arctic.env}
\docType{data}
\title{ Arctic Pollen and associated environmental data}
\description{
  Arctic pollen percent data and associated environmental data
}
\usage{data(arctic.pollen)}
\format{
  \describe{
    \item{arctic.pollen}{A data frame with 828 observations on the percentage of 39 pollen taxa}
    \item{arctic.env}{Environmental data for the pollen sites}
  }
}
\source{                                                    
  Data extracted from North American Pollen Database and New \emph{et al.} (2002) by Fréchette \emph{et al.} (2008). Following Fréchette (Pers. Comm.), three duplicate sites have been deleted.
}
\references{
 Fréchette, B., de Vernal, A., Guiot, J., Wolfe, A. P., Miller, G. H., Fredskild, B., Kerwin, M. W. and Richard, P. J. H. (2008) Methodological basis for quantitative reconstruction of air temperature and sunshine from pollen assemblages in Arctic Canada and Greenland. \emph{Quaternary Science Reviews} \bold{27}, 1197--1216 doi:10.1016/j.quascirev.2008.02.016
}
\examples{
data(arctic.pollen)
data(arctic.env)
}
\keyword{datasets}