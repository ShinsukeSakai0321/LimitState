\name{GeneralTreat}
\alias{GeneralTreat}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
GeneralTreat
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{aa<-GeneralTreat$new(g = "R-S", var = c("R", "S"), dist = c("normal", "normal"), muX = c(200, 100), sigmmaX = c(10, 20))}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Calc()}{Execute reliability analysis}
\item{GetPOF() }{returns POF value}
\item{GetAlpha() }{returns Alpha vector}
\item{GetBeta()}{returns reliability index}
\item{GetDP()}{returns Design Points}
\item{GetPSF()}{returns PSF vector}
}
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
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
#Example1
aa <- GeneralTreat$new(g="R-S",var=c("R","S"),dist=c("normal","normal"),
                                                    muX=c(200,100),sigmmaX=c(10,20))
aa$Calc()
aa$GetPOF()
# 3.872108e-06
aa$GetAlpha()
# 0.4472136 -0.8944272
aa$GetBeta()
# 4.472136
aa$GetDP()
# 180 180
aa$GetPSF()
# 0.9 1.8
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
aa <-GeneralTreat$new(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
aa$Calc()
a1 <- aa$GetBeta()
dist[7] <- "lognormal"
aa <-GeneralTreat$new(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
aa$Calc()
a2 <- aa$GetBeta()
dist <- c("lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"normal" )
aa <-GeneralTreat$new(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
aa$Calc()
a3 <- aa$GetBeta()
dist <- c("lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"lognormal" )
aa <-GeneralTreat$new(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
aa$Calc()
a4 <- aa$GetBeta()
c(a1,a2,a3,a4)
#expected output
# 3.833028 3.761254 4.387684 4.090647
##########################
#Example in text p.Basic2.30
aa<-GeneralTreat$new(g="R-S",var=c("R","S"),dist=c("normal","normal"),
                                 muX = c(400, 300), sigmmaX = c(20, 30))
aa$Calc()
aa$GetBeta()
#correct answer:  2.773501
##########################
#Example in text p.Basic2.37
g<-"Sb*A-P"
var <- c("Sb","A","P")
muX <- c(400, 100*10*10, 3000e3)
sigmmaX <- c(20, 5*10*10, 300e3)
dist <- c("normal", "normal", "normal" )
aa <-GeneralTreat$new(g=g,var=var,dist=dist,muX=muX,sigmmaX=sigmmaX)
aa$Calc()
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
############################
#Example in text Basic3.12
dist <- c("lognormal", "normal" )
aa<-GeneralTreat$new(g="R-S",var=c("R","S"),dist=dist, muX = c(40, 25), sigmmaX = c(4, 5))
aa$Calc()
aa$GetBeta()
aa$GetPOF()
#> aa$GetBeta()
#[1] 2.377765
#> aa$GetPOF()
#[1] 0.008708954
}
\keyword{datasets}
