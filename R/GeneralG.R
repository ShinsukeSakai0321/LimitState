GeneralG <- R6::R6Class("GeneralG",
                        inherit=Lbase,
                        public = list(
                          initialize = function(gg,var){
                            private$n <- length(var)
                            super$initialize(private$n)
                            private$gg <- gg
                            private$var <- var
                          },
                          gcalc = function(){
                            expr <- parse(text=private$gg)
                            X <- super$GetX()
                            for(i in 1:private$n){
                              str1 <- paste(private$var[i],"<- X[",sep="")
                              str1 <- paste(str1,as.character(i),sep="")
                              str1 <- paste(str1,"]",sep="")
                              exp_sub <- parse(text=str1)
                              eval(exp_sub)
                            }
                            super$SetG(eval(expr))
                          },
                          dGdXcalc = function(){
                            X <- super$GetX()
                            dGdX <- super$GetdGdX()
                            for(i in 1:private$n){
                              str1 <- paste(private$var[i],"<- X[",sep="")
                              str1 <- paste(str1,as.character(i),sep="")
                              str1 <- paste(str1,"]",sep="")
                              exp_sub <- parse(text=str1)
                              eval(exp_sub)
                            }
                            expr <- parse(text=private$gg)
                            for(i in 1:private$n){
                              dstr1 <- D(expr,private$var[i])
                              dGdX[i] <- eval(dstr1)
                            }
                            super$SetdGdX(dGdX)
                          }
                        ),
                        private = list(
                          gg = "",
                          var = c(""),
                          n = 0
                        )
                        )
GeneralTreat <- R6::R6Class("GeneralTreat",
                            inherit=LSFM,
                            public = list(
                              initialize = function(g,var,dist,Mu,sigmmaX){
                                n <- length(var)
                                private$g <- g
                                private$var <- var
                                super$initialize("GeneralTreat",n,Mu,sigmmaX,dist)
                              },
                              Calc = function(){
                                gg <- LimitState::GeneralG$new(private$g,private$var)
                                super$DefineG(gg)
                                super$RF()
                              }
                            ),
                            private = list(
                              g = "",
                              var = c("")
                            )
)
reliability <- function(g="R-S",var=c("R","S"),
                        dist=c("normal","normal"),muX=c(200,100),sigmmaX=c(10,20)){
  aa <- LimitState::GeneralTreat$new(g,var,muX,sigmmaX,dist)
  aa$Calc()
  aa
}


