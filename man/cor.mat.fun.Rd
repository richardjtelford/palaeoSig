\name{cor.mat.fun}
\alias{cor.mat.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generates correlation matrix 
%%  ~~function to do ... ~~
}
\description{ Generates a correlation matrix for the environmental variables generated in \code{\link{make.env}} and for correlated species optima in \code{\link{species}}. Only used when correlated environmental variables or optima are generated.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cor.mat.fun(ndim, cors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ndim}{ Number of environmental variables that are subsequently generated with \code{\link{make.env}}.
%%     ~~Describe \code{ndim} here~~
}
\item{cors}{ List of correlations between environmental variables. Each element of the list consists of three numbers, the first two numbers indicate the variables that are correlated, the third number is the correlation coefficient. If correlations between two variables are omitted the correlation remains 0.      
%%     ~~Describe \code{ndim} here~~
}
}
\value{ Correlation matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{ Mathias Trachsel
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{make.env}}, \code{\link{species}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
correlations <- list(c(1, 2, 0.5), c(1, 4, 0.1), c(2, 5, 0.6))
cor.mat <-  cor.mat.fun(5, correlations)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
