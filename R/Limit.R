# Calss for Limit State Function Method
library(R6)
# Base class for the description of limit state function
Lbase <- R6::R6Class("Lbase",
                 public = list(
                   initialize = function(n){
                     private$n <- n
                     private$X <- numeric(n)
                     private$dGdX <- numeric(n)
                   },
                   GetN = function(){private$n},
                   GetX = function(){private$X},
                   SetX = function(X){
                     private$X <- X
                   },
                   GetG = function(){private$g},
                   SetG = function(g){private$g <- g},
                   GetdGdX = function(){private$dGdX},
                   SetdGdX = function(x){
                     private$dGdX <- x
                   }
                 ),
                 private = list(
                   n = 0,   # Number of variables
                   X = numeric(0),
                   dGdX = numeric(0),
                   g = 0.0
                 )
)

LSFM <- R6::R6Class("LSFM",
                public = list(
                  name=NA,
                  initialize = function(name,n,Mu,sigmmaX,dist){
                    if(!missing(name))self$name <- name
                    private$n <- n  #Number of parameters
                    private$muX <- Mu
                    private$sigmmaX <- sigmmaX
                    private$dist <- dist
                    private$alphai <- numeric(n)
                    if(private$dist[1]=="normal")private$Distr=c(GNormal$new(private$muX[1],private$sigmmaX[1]))
                    if(private$dist[1]=="lognormal")private$Distr=c(GLognormal$new(private$muX[1],private$sigmmaX[1]))
                    if(private$dist[1]=="gumbel")private$Distr=c(GGumbel$new(private$muX[1],private$sigmmaX[1]))
                    if(private$dist[1]=="weibull")private$Distr=c(GWeibull$new(private$muX[1],private$sigmmaX[1]))
                    for (i in 2:n){
                      if(private$dist[i]=="normal")private$Distr=c(private$Distr,GNormal$new(private$muX[i],private$sigmmaX[i]))
                      if(private$dist[i]=="lognormal")private$Distr=c(private$Distr,GLognormal$new(private$muX[i],private$sigmmaX[i]))
                      if(private$dist[i]=="gumbel")private$Distr=c(private$Distr,GGumbel$new(private$muX[i],private$sigmmaX[i]))
                      if(private$dist[i]=="weibull")private$Distr=c(private$Distr,GWeibull$new(private$muX[i],private$sigmmaX[i]))
                    }
                  },
                  GetN = function()private$n,
                  GetMu = function()private$muX,
                  SetMu = function(aa){private$muX <- aa},
                  RF = function(){
                    betaold = 40
                    delta = 1e-6
                    munormX <- numeric(private$n)
                    sigmmanormX <- numeric(private$n)
                    # step 1 : Define the appropriate performance function
                    # step 2 : Assume initial values of the design point
                    X <- private$muX ## * rnum
                    # step 3 : Compute the mean and standard deviation at the
                    #          design point of the equivalent normal distribution
                    for (i in 1:100){
                      for (j in 1:private$n){
                        Valu <- private$Distr[[j]]$Eq(X[j])
                        munormX[j] <- Valu[1]
                        sigmmanormX[j] <- Valu[2]
                      }
                      Xdush <- (X - munormX) / sigmmanormX
                      Xdush.old <- Xdush
                      # step 4 : Compute the limit state function g and partial derivative
                      #          dg/dXi evaluated at design point
                      private$lim$SetX(X)
                      private$lim$gcalc()
                      private$lim$dGdXcalc()
                      g <- private$lim$GetG()
                      dgdX <- private$lim$GetdGdX()
                      # step 6 : Compute the new values for the design point
                      #          in the equivalent standard normal space
                      dgdXdush <- dgdX * sigmmanormX
                      A <- 1 / sum(dgdXdush * dgdXdush) * (sum(dgdXdush * Xdush) - g)
                      Xdush <- A * dgdXdush
                      private$alphai <- dgdXdush / sqrt(sum(dgdXdush * dgdXdush))
                      Xdush.new <- Xdush
                      # step 7 : Compute the distance to this new design point
                      betanew <- sqrt(sum(Xdush * Xdush))
                      hantei <- is.nan(betanew)
                      # step 8 : Compute the new values for the design point
                      #          in the priginal space
                      if(is.nan(betanew)) betaold <- betaold
                      else {
                        X <- munormX + sigmmanormX * Xdush
                        if(abs(betaold - betanew) < delta) break
                        betaold <- betanew
                      }
                      ## for moddegied RF
                      X.t <- munormX + sigmmanormX * Xdush.old
                      private$lim$SetX(X.t)
                      private$lim$gcalc()
                      g.hantei <- private$lim$GetG()
                      #if(is.na(g.hantei)){
                      #cat("X=",X.t,"\r\n")
                      #cat("munormX=",munormX,"\r\n")
                      #cat("sigm=",sigmmanormX,"\r\n")
                      #cat("xdush=",Xdush.old,"\r\n")
                      #}
                      deltan <- 50
                      dXdush <- (Xdush.new - Xdush.old) / deltan
                      for(i1 in 1:deltan){
                        Xdush.t <- Xdush.old + i1 * dXdush
                        X.t <- munormX + sigmmanormX * Xdush.t
                        private$lim$SetX(X.t)
                        private$lim$gcalc()
                        g.hantein <- private$lim$GetG()
                        if(is.na(g.hantei))return()   #added by S.sakai 2017.4.13
                        if(g.hantei^2 > g.hantein^2){
                          Xdush <- Xdush.t
                          g.hantei <- g.hantein
                        }
                        else{
                          g.hantei <- g.hantei
                        }
                      }
                      X <- munormX + sigmmanormX * Xdush
                    }
                    private$beta <- betanew
                    private$POF <- pnorm(betanew, lower.tail = FALSE)
                    private$DP <- X
                    private$ncon <- i
                  },
                  GetSigm = function()private$sigmmaX,
                  GetBeta = function()private$beta,
                  GetAlpha = function()private$alphai,
                  GetPOF = function()private$POF,
                  GetDP = function()private$DP,
                  GetConv = function()private$ncon,
                  GetPSF = function(){
                    private$DP/private$muX
                  },
                  DefineG = function(glim){private$lim <- glim},
                  GetG = function(){
                    #added by S.sakai 2017.4.13
                    #平均値位置でのg()の計算値を返す，負のとき破損領域にある
                    X <- private$muX
                    private$lim$SetX(X)
                    private$lim$gcalc()
                    private$lim$GetG()
                  }
                ),
                private = list(
                  n = 0,   # Number of variables
                  ncon = 0, #number of convergence
                  muX = c(0),
                  sigmmaX = c(0),
                  alphai = c(0),
                  DP = c(0),
                  dist = c(" "),
                  Distr = c(Dbase),
                  lim = Lbase,   #Object for limit state function
                  beta = 0.0,
                  POF = 0.0
                )
)


