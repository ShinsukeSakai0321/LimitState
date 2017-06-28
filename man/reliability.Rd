\name{reliability}
\alias{reliability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
reliability
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function executes reliability analysis by limit state function method. Returned object is GeneralTreat class which is inheritted from LSFM class.
}
\usage{
reliability(g = "R-S", var = c("R", "S"), dist = c("normal", "normal"), muX = c(200, 100), sigmmaX = c(10, 20))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{ character expression of limit state fuction
%%     ~~Describe \code{g} here~~
}
  \item{var}{ array of variables written in character type
%%     ~~Describe \code{var} here~~
}
  \item{dist}{ array of distribution type written in character type: one of "normal","lognormal","gumbel","weibull"
%%     ~~Describe \code{dist} here~~
}
  \item{muX}{ array of mean value for variables
%%     ~~Describe \code{muX} here~~
}
  \item{sigmmaX}{array of sd. value for variables
%%     ~~Describe \code{sigmmaX} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{GetPOF() }{returns POF value}
\item{GetAlpha() }{returns Alpha vector}
\item{GetBeta()}{returns reliability index}
\item{GetDP()}{returns Design Points}
\item{GetPSF()}{returns PSF vector}
}
%% ...
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Shinsuke Sakai
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
Returned object is GeneralTreat class which is inheritted from LSFM class. As for public function of LSFM class, see LSFM class by help command.
}
\examples{
#Example1
aa <- reliability()
aa$GetPOF()
aa$GetAlpha()
aa$GetBeta()
aa$GetDP()
aa$GetPSF()
# correct answer  3.872108e-06
#
# Example2
#"Probability, Reliability and Statistical Methods in
#  Engineering Design"
# Achintya Haldar  & Sankaran Mahadevan
#  P.218 Table 7.5, 7.6
g<-"As*fy*d*(1.0-eta*As*fy/b/d/fcd)-M"
var <- c("As","fy","fcd","b","d","eta","M")
muX <- c(1.56, 47.7, 3.5, 8.0, 13.2, 0.59, 326.25)
covX <- c(0.036, 0.15, 0.21, 0.045, 0.086, 0.05, 0.17)
sigmmaX <- muX*covX
dist <- c("normal", "normal", "normal", "normal" ,"normal" ,"normal" ,"normal" )
aa <-reliability(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
a1 <- aa$GetBeta()
dist[7] <- "lognormal"
aa <-reliability(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
a2 <- aa$GetBeta()
dist <- c("lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"normal" )
aa <-reliability(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
a3 <- aa$GetBeta()
dist <- c("lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"lognormal" )
aa <-reliability(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
a4 <- aa$GetBeta()
c(a1,a2,a3,a4)
#expected output
# 3.833028 3.761254 4.387684 4.090647
#Example3
# Checking lack of invariance
aa<-reliability()
bb<-reliability(g="1-R/S")
cc<-reliability(g="log(R)-log(S)")
aa$GetBeta()
bb$GetBeta()
cc$GetBeta()
#Example in text p.Basic2.30
aa<-reliability( muX = c(400, 300), sigmmaX = c(20, 30))
aa$GetBeta()
#correct answer:  2.773501
#Example in text p.Basic2.37
g<-"Sb*A-P"
var <- c("Sb","A","P")
muX <- c(400, 100*10*10, 3000e3)
sigmmaX <- c(20, 5*10*10, 300e3)
dist <- c("normal", "normal", "normal" )
aa <-reliability(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
aa$GetBeta()
aa$GetAlpha()
aa$GetDP()
#Correct answer
#> aa$GetBeta()
#[1] 2.458772
#> aa$GetAlpha()
#[1]  0.4696243  0.4696243 -0.7476002
#> aa$GetDP()
#[1]     376.906    9422.650 3551453.579
#Example in text p.Basic2.44
aa<-reliability( muX = c(400, 300), sigmmaX = c(20, 30))
aa$GetPSF()
#[1] 0.9230769 1.2307692
#Example in text Basic3.12
dist <- c("lognormal", "normal" )
aa<-reliability(dist=dist, muX = c(40, 25), sigmmaX = c(4, 5))
aa$GetBeta()
aa$GetPOF()
#> aa$GetBeta()
#[1] 2.377765
#> aa$GetPOF()
#[1] 0.008708954
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
