\name{predict.agelme}
\alias{predict.agelme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predicts the Calibrated age }
\description{
  This function uses the output from 'agelme' to predict the Calibrated ages for specified depths.
}
\usage{
\method{predict}{agelme}(object, v = 1, depth,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An 'agelme' model }
  \item{v}{Using constant (1) or mu (2) variance }
  \item{depth}{ A vector of the depths to be predicted }
  \item{\dots}{Other arguments, currently unused.}
}

\value{
   A list with three items
   \itemize{
     \item{v} Whether constant variance or mu variance used.
     \item{fit} A data.frame of the predictions including 95\% confidence interval.
     \itemize{
      \item{Depth}{The depths for the predicted ages}
      \item{Estage}{Predicted age}
      \item{Lowlim}{Lower 95\% confidence interval}
      \item{Upplim}{Upper 95\% confidence interval}
      \item{Tsd}{Total standard deviation}
    }
    \item{data} A data.frame containing the age and depth information of the radiocarbon dates.

  }
}

\references{ Heegaard, E., Birks, HJB. & Telford, RJ. 2005. Relationships between calibrated ages and depth in stratigraphical sequences: an estimation procedure by mixed-effect regression. The Holocene 15: 612-618}
\author{ Einar Heegaard <einar.heegaard@bio.bui.no }


\examples{
data(STOR)

fit.mod <- with(STOR,agelme(depthup,depthdo,cageup,cagedo))

#Predicting using the constant variance model,
#for each cm between 70 and 400 cm.
fit.pre <- predict(fit.mod,1,70:400)
plot(fit.pre)
}
