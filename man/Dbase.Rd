\name{Dbase}
\alias{Dbase}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
%%  ~~function to do ... ~~
Dbase
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Base class for the expression of distribution. The main function is to evaluate the equivalent normal distribution at specific estimation point.
}
\usage{
Arbitrary distribution is expressed by inheritting this base class
}

\section{Public function}{
%%  ~~ If necessary, more details than the description above ~~
Standard reference class is R6 class \cr
Public function
\tabular{lll}{
\bold{function} \tab \bold{contents} \tab \bold{variable(s)}\cr
initialize(mu,sigmma) \tab initializes the object \tab mu:mean value,sigmma:std. value\cr
Eq(X) \tab returns (mean,std.) list of equivalent norm. distribution \tab X:estimaion point\cr
GetMeanSig() \tab returns (mu,std.) list of the distribution  \tab \cr
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Shinsuke Sakai
}
\note{
%%  ~~further notes~~
Gnormal, GLognormal, GGumbel, GWeibull have been deveoped as the inheritted classes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Definition of class for log normal distribution
GLognormal <- R6Class("GLognormal",
                      inherit=Dbase,
                      public = list(
                        Eq = function(X){
                          zeta <- sqrt(log(1+ (private$sigmmaX / private$muX)**2))
                          lambda <- log(private$muX) - 0.5 * zeta*zeta
                          phi <- dnorm(qnorm(plnorm(X, lambda, zeta)))
                          fXi <- dlnorm(X, lambda, zeta)
                          sigm <- phi / fXi
                          mu <- X - qnorm(plnorm(X, lambda, zeta)) * sigm
                          c(mu, sigm)
                        }
                      )
)
#main program
dist <- GLognormal$new(4,1)
dist$Eq(1.0)
  }


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
