\name{GGumbel}
\alias{GGumbel}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Class definition to control gumbel distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("GGumbel")}
\format{
  The format is:
Class 'R6ClassGenerator' <GGumbel> object generator
  Inherits from: <Dbase>
  Public:
    Eq: function (X)
    Param: function ()
    SetParam: function (mu, eta)
    MuSigm: function (mu, eta)
    clone: function (deep = FALSE)
  Private:
    eta: 0
    mu: 0
  Parent env: <environment: namespace:LimitState>
  Locked objects: TRUE
  Locked class: FALSE
  Portable: TRUE
 - attr(*, "name")= chr "GGumbel_generator"
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
gum<-GGumbel$new(mu=1.0,sigmma=0.1)
res<-gum$Param()
cat("Shape=",res[1],",Location=",res[2])
#example: giving loc and scale parameters and returning mu,sigm
loc<-0.05; shape<-0.8
res2<-gum$MuSigm(mu=loc,eta=shape)
cat("mu=",res2[1],",sigm=",res2[2])
#Calculating equivalent norm distribution
gum<-GGumbel$new(mu=1.0,sigmma=0.1)
xval<-0.5
res3<-gum$Eq(X=xval)
cat("mu_eq=",res3[1],",sig_eq=",res3[2])
}
\keyword{datasets}
