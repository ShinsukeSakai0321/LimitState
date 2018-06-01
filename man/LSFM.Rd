\name{LSFM}
\alias{LSFM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
LSFM
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Base class for the analysis using limit state function method.
}
\usage{
Arbitrary limit state function subjects are inherited from this base class
}

\section{Public function}{
%%  ~~ If necessary, more details than the description above ~~
Standard reference class is R6 class \cr
Public function
\tabular{lll}{
\bold{function} \tab \bold{contents} \tab \bold{variable(s)}\cr
initialize(name,n,Mu,sigmmaX,dist) \tab initializes the object \tab name: name of the subject,n:number of variables, Mu:list of mean values of variables, sigmmaX: list of std. values of variables,dist: list of distribution names(characters)\cr
GetN() \tab returns the number of variables \tab \cr
DefineG(glim) \tab difines control class for G-calculation \tab glim: object of g-definition class which is inheritted class of Lbase \cr
GetMu() \tab returns  list of mean values of variables \tab \cr
GetSigm() \ tab returns list of std. values of variables \tab \cr
GetBeta() \tab retuns reliability index \tab \cr
GetAlpha() \tab returns the list of differential vector \tab \cr
GetPOF() \tab returns value of probability of failure \tab \cr
GetDP() \tab returns list of DP vector \tab \cr
GetConv() \tab returns number of convergence \tab \cr
GetG() \tab returns g value at the mean vector of the variables \tab \cr
SetMu(Xm) \tab sets mean values of variables from list Xm  \tab Xm:list of mean values of variables\cr
RF() \tab searches the design point usein modified Rackwitz-Fiessler algorithm \tab  \cr

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
Shinsuke Sakai, Takuyo Kaida
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Most simple example
# Definition of G and dGdX
GRS <- R6::R6Class("GRS",
                       inherit=Lbase,
                       public=list(
                         gcalc = function(){
                           X <- super$GetX()
                           R <- X[1]
                           S <- X[2]
                           g<- R-S
                           super$SetG(g)
                         },
                         dGdXcalc = function(){
                           X <- super$GetX()
                           dGdX <- super$GetdGdX()
                           dGdX[1] <- 1
                           dGdX[2] <- -1
                           super$SetdGdX(dGdX)
                         }
                       )
                       )
# Reliability Analysis
Sexample <- R6::R6Class("Sexample",
                        inherit=LSFM,
                        public= list(
                          Calc=function(){
                            super$DefineG(GRS$new(super$GetN()))
                            super$RF()
                          }
                        )
                        )
####main###
muX<-c(200,100)
sigmmaX<-c(10,20)
dist <- c("normal", "normal")
exam <- Sexample$new("Sexample",2,muX,sigmmaX,dist)
exam$Calc()
exam$GetPOF()
# correct answer  3.872108e-06


#Example of Lbase must be executed before this example.
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
                           dGdX[1] <- eval(D(dd,"As"))
                           dGdX[2] <- eval(D(dd,"fy"))
                           dGdX[3] <- eval(D(dd,"fcd"))
                           dGdX[4] <- eval(D(dd,"b"))
                           dGdX[5] <- eval(D(dd,"d"))
                           dGdX[6] <- eval(D(dd,"eta"))
                           dGdX[7] <- eval(D(dd,"M"))
                           super$SetdGdX(dGdX)
                         }
                       )

)
PRSMED <- R6::R6Class("PRSMED",
                     inherit=LSFM,
                     public = list(
                       Calc = function(){
                         super$DefineG(GPRSMED$new(super$GetN()))
                         super$RF()
                       }
                     )
)
############### Main Program #######################
##################### Input Data (mean value and COV) ##################
Ay <- 1.56 ; Aycov <- 0.036
fy <- 47.7; fycov <- 0.15
fcd <- 3.5; fcdcov <- 0.21
b <- 8.0; bcov <- 0.045
d <- 13.2; dcov <- 0.086
eta <- 0.59; etacov <- 0.05
M <- 326.25; Mcov <- 0.17

##list of mean values
muX <-  c(Ay,fy,fcd,b,d,eta,M)
##list of COV values
COVX <- c(Aycov,fycov,fcdcov,bcov,dcov,etacov,Mcov)
sigmmaX <- muX * COVX
##sets names of distribution
dist <- c("normal", "normal", "normal", "normal" ,"normal" ,"normal" ,"normal" )

robj <- PRSMED$new("PRSMDE",7,muX,sigmmaX,dist)
robj$Calc()
cat('Reliability index beta = ', robj$GetBeta(), '\n')
dist[7] <- "lognormal"
robj <- PRSMED$new("PRSMDE",7,muX,sigmmaX,dist)
robj$Calc()
cat('Reliability index beta = ', robj$GetBeta(), '\n')
dist <- c("lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"normal" )
robj <- PRSMED$new("PRSMDE",7,muX,sigmmaX,dist)
robj$Calc()
cat('Reliability index beta = ', robj$GetBeta(), '\n')
dist <- c("lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"lognormal" )
robj <- PRSMED$new("PRSMDE",7,muX,sigmmaX,dist)
robj$Calc()
cat('Reliability index beta = ', robj$GetBeta(), '\n')

#######################################################
#### expected results of the example program ####
#Reliability index beta =  3.833028
#Reliability index beta =  3.761254
#Reliability index beta =  4.387684
#Reliability index beta =  4.090647
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
