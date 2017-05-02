#Description of distribution
# Eq(x)      : Equivalent mean and stand.dev. for normal distribution at value of x
library(evd)
library(R6)
Dbase <- R6::R6Class("Dbase",
                 public = list(
                   initialize = function(mu,sigmma){
                     private$muX<-mu
                     private$sigmmaX<-sigmma
                   },
                   GetMeanSig = function(){c(private$muX,private$sigmmaX)},
                   DushCalc = function(X){
                     #under construction
                     private$sigmmaDush <- dnorm(qnorm(self$CDF(X)))/self$PDF(X)
                     private$muDush <- X - (qnorm(self$CDF(X)))*private$sigmmaDush
                   },
                   # Eq(x)      : Equivalent mean and stand.dev. for normal distribution at value of x
                   Eq=function(x)print("Specify Eq()")
                 ),
                 private = list(
                   muX = 0.0,
                   sigmmaX = 0.0,
                   sigmmaDush=0.0,
                   muDush = 0.0
                 )
)
GNormal <- R6::R6Class("GNormal",
                   inherit=Dbase,
                   public = list(
                     Eq = function(X){c(private$muX,private$sigmmaX)}
                   )
)
GLognormal <- R6::R6Class("GLognormal",
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
GGumbel <- R6::R6Class("GGumbel",
                   inherit=Dbase,
                   public = list(
                     Eq = function(X){
                       eta <- sqrt(6.) * private$sigmmaX /pi #shape
                       mu <- private$muX - 0.57722 * eta #loc
                       phi <- dnorm(qnorm(evd::pgumbel(X, loc = mu, scale = eta)))
                       fXi <- evd::dgumbel(X, loc = mu, scale = eta)
                       sigm <- phi / fXi
                       mu <- X - qnorm(evd::pgumbel(X, loc = mu, scale = eta)) * sigm
                       c(mu, sigm)
                     },
                     Param = function(){
                       private$eta <- sqrt(6.) * private$sigmmaX /pi #scale
                       private$mu <- private$muX - 0.57722 * private$eta #loc
                       c(private$mu,private$eta)
                     },
                     SetParam = function(mu,eta){
                       private$sigmmaX <- pi/sqrt(6)*eta
                       private$muX <- mu + 0.57722*eta
                       c(private$muX,private$sigmmaX)
                     },
                     MuSigm = function(mu,eta){
                       c(mu+0.57722*eta,pi*eta/sqrt(6))
                     }
                   ),
                   private = list(
                     eta=0.0,
                     mu=0.0
                   )
)
GWeibull <- R6::R6Class("GWeibull",
                    inherit=Dbase,
                    public = list(
                      Eq = function(X){
                        self$FactorCalc()
                        super$DushCalc(X)
                        c(private$muDush,private$sigmmaDush)
                      },
                      FactorCalc = function(){
                        cov <- private$muX / private$sigmmaX
                        private$ar <- cov ^ (-1.08) #Prof.Ichikawa's theory
                        private$beta <- private$muX/gamma(1.0+1/private$ar)
                      },
                      PDF = function(X){
                        dweibull(X,private$ar,private$beta)
                      },
                      CDF = function(X){
                        pweibull(X,private$ar,private$beta)
                      }
                    ),
                    private = list(
                      ar=0.0,
                      beta=0.0
                    )
)
bb <- GWeibull$new(10,1)
bb$Eq(1)

