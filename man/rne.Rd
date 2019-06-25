\name{rne}
\alias{rne}
\alias{plot.RNE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random, neighbour, environment deletion analysis for transfer function model }
\description{
  Calculates effect of deleting sites from training set at random, from a geographic neighbourhood, or from an environmental neighbourhood. A simple graphical technique for gauging the effect of spatial autocorrelation on the transfer function model.
}
\usage{
rne(y, env, geodist,fun, neighbours, subsets = c(1, 0.75, 0.5, 0.25, 0.1), ...) 
\method{plot}{RNE}(x, which=1, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Community data, or distance object, or distance matrix}
  \item{env}{Environmental variable}
  \item{geodist}{Matrix of geographical distances between sites}
  \item{fun}{Transfer function}
  \item{neighbours}{ Neighbourhood radii }
  \item{subsets}{ Proportion of sites to retain in random deletion }
  \item{...}{Arguments passed to fun or plot}
  \item{x}{RNE object to be plotted}
  \item{which}{Which column of the results to plot eg if more than one WAPLS component is calculated}
  \item{ylim}{Y-limits of the plot}
  
  
  
}
\details{
  Finds the leave-one-out transfer function performance if sites are deleted at random (repeated 10 times to reduce variance in results), from a neighbourhood zone, or by deleting environmentally close sites.
}
\value{
  Returns an RNE object, list with two components
  
  \item{random }{Performance with random deletion.}
  \item{neighbour }{Performance with deletion by neighourhood, or environment}
  
}
\references{ Telford, R. J. and Birks, H. J. B. (2009) Evaluation of transfer functions in spatially structured environments. \emph{Quaternary Science Reviews} \bold{28}: 1309--1316. DOI: 10.1016/j.quascirev.2008.12.020 }
\author{Richard Telford \email{Richard.Telford@bio.uib.no}  }

\examples{
  require(fields)
  require(rioja)
  data(arctic.env)
  data(arctic.pollen)

  \dontshow{
    #using just the first 20 sites so that code runs in an reasonable time
    arctic.dist <- rdist.earth(cbind(arctic.env$Longitude, arctic.env$Latitude)[1:20,], miles = FALSE)
    
    arctic.rne <- rne(y = arctic.pollen[1:20,], env = arctic.env$tjul[1:20], geodist = arctic.dist, fun = MAT, neighbours = c(0,200), subsets = c(1,.5), k = 5)
    plot(arctic.rne)
  }
  \donttest{
    #using just the first 100 sites so code runs in an reasonable time, but still slow
    arctic.d <- rdist.earth(arctic.env[1:100,c("Longitude","Latitude")], miles = FALSE)
    
    arctic.rne <- rne(y = arctic.pollen[1:100,], env = arctic.env$tjul[1:100],
        geodist = arctic.d, fun = MAT, neighbours = c(0,50,200), 
        subsets = c(1,.5,.1), k = 5)
    plot(arctic.rne)
  }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
