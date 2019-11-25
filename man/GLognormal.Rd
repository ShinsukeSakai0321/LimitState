\name{GLognormal}
\alias{GLognormal}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Class definition to control lognormal distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{lnorm<-GLognormal$new(mu=350.25,sigmma=105.075)}
\format{
SetParam(lambda=lambda,zeta=zeta)
    Setting lognormal parameters and calculating mu and sigm. lambda: mean of lnorm, zeta: std. of lnorm
    Param(): Calculating lambda and mu
    Eq(x): Calculating equivalent normal distribution for adujusting point x
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
#fundamental usage
lnorm<-GLognormal$new(mu=350.25,sigmma=105.075)
res<-lnorm$Param()
cat("lambda=",res[2],",zeta=",res[1])
# lambda= 5.815558 ,zeta= 0.2935604
#example: giving lambda and eta parameters and returning mu,sigm
lambda<- 5.815558 ;zeta<- 0.2935604
res2<-lnorm$SetParam(lambda=lambda,zeta=zeta)
cat("mu=",res2[1],",sigm=",res2[2])
# mu= 350.2499 ,sigm= 105.075
#Calculating equivalent norm distribution
xval<-0.5
res3<-lnorm$Eq(X=xval)
cat("mu_eq=",res3[1],",sig_eq=",res3[2])
# mu_eq= 3.754353 ,sig_eq= 0.1467802
}
\keyword{datasets}
