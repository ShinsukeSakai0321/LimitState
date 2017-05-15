\name{Lbase}
\alias{Lbase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Lbase
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Base class for the description of limit state function
}
\usage{
Arbitrary limit state functions defined by users are inherited from this base class
}
%- maybe also 'usage' for other objects documented here.

\section{Public function}{
Standard reference class is R6 class \cr
Public function
\tabular{lll}{
\bold{function} \tab \bold{contents} \tab \bold{variable(s)}\cr
initialize(n) \tab initializes the object \tab n:number of variables\cr
GetN() \tab returns the number of variables \tab \cr
GetX() \tab returns the variable list \bold{X} \tab \cr
SetX(x) \tab sets variable list \bold{x} to the variable list \bold{X} \tab x: variable list\cr
GetG() \tab returns calculated g value \tab \cr
SetG(g) \tab sets calculated g-value to the variable \tab g:calculated g \cr
GetdGdX() \tab returns differential variable value list \bold{dGdX} \tab \cr
SetdGdX(x) \tab sets differintial variable list \bold{dGdX} \tab x: \bold{dGdX list}
}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Shinsuke Sakai, Takuyo Kaida
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Example
#"Probability, Reliability and Statistical Methods in
#  Engineering Design"
# Achintya Haldar  & Sankaran Mahadevan
#  P.218 Table 7.5, 7.6
GPRSMED <- R6::R6Class("GPRSMED",
                   inherit=Lbase,
                   public = list(
                     gcalc = function(){
                       X <- super$GetX()
                       As <- X[1]
                       fy <- X[2]
                       fcd <- X[3]
                       b <- X[4]
                       d <- X[5]
                       eta <- X[6]
                       M <- X[7]
                       g <- As*fy*d*(1.0-eta*As*fy/b/d/fcd)-M
                       super$SetG(g)
                     } ,
                     dGdXcalc = function(){
                       X <- super$GetX()
                       dGdX <- super$GetdGdX()
                       As <- X[1]
                       fy <- X[2]
                       fcd <- X[3]
                       b <- X[4]
                       d <- X[5]
                       eta <- X[6]
                       M <- X[7]
                       dGdX[1] <- fy * d * (1 - eta * As * fy/b/d/fcd) - As * fy * d * (eta * fy/b/d/fcd)
                       dGdX[2] <- As * d * (1 - eta * As * fy/b/d/fcd) - As * fy * d * (eta * As/b/d/fcd)
                       dGdX[3] <- As * fy * d * (eta * As * fy/b/d/fcd^2)
                       dGdX[4] <- As * fy * d * (eta * As * fy/b^2/d/fcd)
                       dGdX[5] <- As * fy * (1 - eta * As * fy/b/d/fcd) + As * fy * d * (eta *
                                                                                                   As * fy/b/d^2/fcd)
                       dGdX[6] <- -(As * fy * d * (As * fy/b/d/fcd))
                       dGdX[7] <- -1
                       super$SetdGdX(dGdX)
                     }
                   )

)
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

As the application of the following example of GPRSMED, see example of help command LSFM
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
