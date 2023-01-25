% --- Source file: man/agelme.Rd ---
\name{agelme}
\alias{agelme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimation of the relationship between Calibrated age and depth }
\description{
  Estimates the relationship of Calibrated age and depth for palaeorecords. The function uses a smooth spline
  of the mgcv library by Simon Wood. It produces predicted confidence interval for the relationship approximating 
  a mixed effect model, as there are two levels of uncertainty, i.e. within dated object and between dated objects.
}
\usage{
agelme(depup, depdo, bpup, bpdo, use, weights=c(1,rep(0,length(depup)-1)), 
    vspan=1, k=length(depup)-1, m=2, diagnostic=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{depup}{ The upper depths of the dated slides }
  \item{depdo}{ The lower depths of the dated slides }
  \item{bpup}{ The younger calibrated ages of the dated slides }
  \item{bpdo}{ The older calibrated ages of the dated slides }
  \item{use}{ Logical vector of dates to include in the model. Default is to use all.}
  \item{weights}{ Weights to be used for the estimation, default is fixed top-layer followed by inverse variance of within dated object }
  \item{vspan}{ The span to be used for the diagnostic plots, default span = 1 }
  \item{k}{ Number of base function to start the shrinkage in the gam estimation procedure }
  \item{m}{The order of penalty for the term, i.e. the degree of continuity at the knots (default, m = 2 gives cubic smooth spline) }
  \item{diagnostic}{Logical, should diagnostic plots be made.}
}
\details{
  Note that the fixation of the top layer is done by a weight = 1, whereas the other weights follows inverse variance within object.
  
  The diagnostic plots is used to check the quality of the estimation and to see if there is a need for an assumption of between
  object variance proportional to mean. The latter however is rarely encountered for palaeodata.
}
\value{
  \item{tdf }{Degrees of freedom used by the cubic smooth spline, a vector with first value for constant variance and second vector for variance equal to mu.}
  \item{weights }{A vector of the weights used by the cubic smooth spline }
%  \item{Constant}{A matrix with the numerical results for the dated points using a constant variance}
%  \item{Muvar}{A matrix with the numerical results for the dated points using variance equal to mu}
  \item{RES}{A vector of the Residual sum of squares}
  \item{Models}{A list with the models from the cubic smooth spline, constant and mu variance, respectively}
  \item{Data}{A data.frame including the data used for the estimation}
}
\references{ Heegaard, E., Birks, HJB. & Telford, RJ. 2005. Relationships between calibrated ages and depth in stratigraphical sequences: an estimation procedure by mixed-effect regression. The Holocene 15: 612-618}
\author{ Einar Heegaard <einar.heegaard@bio.uib.no> }

\examples{
data(STOR)

fit.mod <- with(STOR,agelme(depthup,depthdo,cageup,cagedo))

#Predicting using the constant variance model,
#for each cm between 70 and 400 cm.
fit.pre <- predict(fit.mod,1,70:400)
plot(fit.pre)

}
