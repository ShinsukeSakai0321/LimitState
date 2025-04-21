
# LimitState

Software for Advanced First Order Second Moment Reliability Method

# Author
Shinsuke Sakai   
Yokohama National University

## Features



## Installation

You can install the package via pip:

```bash
pip install LimitState
```
## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Usage
### When defining and using a limit state function in a class
Example for RS model
```python
from LimitState import LimitState as ls
class GRS(ls.LSFM):  #Inherit from the LSFM class
    def __init__(self,n,muX,sigmmaX,dist):
        super().__init__(n,muX,sigmmaX,dist)
    def gcalc(self,X): #Virtural function for calculating g value
        R=X[0]
        S=X[1]
        g=R-S  #Definition of limit statue function
        return g  #Make sure to  return the g value.
    def dGdXcalc(self,X): #Virtural function for calculating the gradient vector of g
        R=X[0]
        S=X[1]
        dGdX=[0]*len(X) #Create a list of the same size as X.
        dGdX[0] =1.0
        dGdX[1] =-1.0
        return dGdX #Make sure to  return the gradient vector of g.
dist=['normal','normal'] #分布形はweibull,normal,lognormal,gumbel,uniformのいずれか
muX=[200,100]
sigmmaX=[10,20]
n=2
rs=GRS(n,muX,sigmmaX,dist) #インスタンスの生成
rs.RFn()
rs.GetBeta()
# 4.47213595499958
```

### In the RS model, when g is defined as a string
```python
g="r-s"#Definition of limit state function
var =['r','s']#Definition of random variables
muX =[200,100]#Mean value of random variables
sigmmaX=[10,20]#Standard deviation of random variables
dist =["normal", "normal"]#Definition of distribution type of random variables
aa =ls.GeneralG(g,var,len(var),muX,sigmmaX,dist)
aa.RFn()
aa.GetBeta()
# 4.47213595499958
```

### Another example for defining g as a string
Verify using the example problem described in the following reference.   
**Reference**  "Probability, Reliability and Statistical Methods
in　Engineering Design" Achintya Haldar & Sankaran Mahadevan
P.218 Table 7.5, 7.6
```python
import numpy as np
g="As*fy*d*(1.0-eta*As*fy/b/d/fcd)-M"#Definition of limit state function
var =["As","fy","fcd","b","d","eta","M"]
muX =[1.56, 47.7, 3.5, 8.0, 13.2, 0.59, 326.25]
covX =[0.036, 0.15, 0.21, 0.045, 0.086, 0.05, 0.17]
sigmmaX = list(np.array(covX)*np.array(muX))
dist =["normal", "normal", "normal", "normal" ,"normal" ,"normal" ,"normal" ]
aa =ls.GeneralG(g,var,len(var),muX,sigmmaX,dist)
aa.RFn()
a1=aa.GetBeta()
dist[6]="lognormal"
aa =ls.GeneralG(g,var,len(var),muX,sigmmaX,dist)
aa.RFn()
a2=aa.GetBeta()
dist=["lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"normal"]
aa =ls.GeneralG(g,var,len(var),muX,sigmmaX,dist)
aa.RFn()
a3=aa.GetBeta()
dist=["lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"lognormal"]
aa =ls.GeneralG(g,var,len(var),muX,sigmmaX,dist)
aa.RFn()
a4=aa.GetBeta()
[a1,a2,a3,a4]
#[3.8330281658128675, 3.7612536656241864, 4.387684391077286, 4.090647398473586]
```

### When defining the limit state function as a class and evaluating the derivative of g using symbolic processing
```python
from sympy import *
import numpy as np
class Gmulti(ls.LSFM):
    def __init__(self,n,muX,sigmmaX,dist):
        As,fy,fcd,b,d,eta,M = symbols('As fy fcd b d eta M')
        g=As*fy*d*(1.0-eta*As*fy/b/d/fcd)-M
        #Precompute the derivative expression using symbolic computation.
        self.gg=str(g)
        #Initialize the parent class after defining self.gg
        super().__init__(n,muX,sigmmaX,dist)
        self.d0=str(diff(g,As))
        self.d1=str(diff(g,fy))
        self.d2=str(diff(g,fcd))
        self.d3=str(diff(g,b))
        self.d4=str(diff(g,d))
        self.d5=str(diff(g,eta))
        self.d6=str(diff(g,M))
    def gcalc(self,X):
        As=X[0]
        fy=X[1]
        fcd=X[2]
        b=X[3]
        d=X[4]
        eta=X[5]
        M=X[6]
        g=eval(self.gg)
        return g
    def dGdXcalc(self,X):
        As=X[0]
        fy=X[1]
        fcd=X[2]
        b=X[3]
        d=X[4]
        eta=X[5]
        M=X[6]
        dGdX=[0]*len(X)
        #Evaluation of derivative vector
        dGdX[0]=eval(self.d0)
        dGdX[1]=eval(self.d1)
        dGdX[2]=eval(self.d2)
        dGdX[3]=eval(self.d3)
        dGdX[4]=eval(self.d4)
        dGdX[5]=eval(self.d5)
        dGdX[6]=eval(self.d6)
        return dGdX
dist =["normal", "normal", "normal", "normal" ,"normal" ,"normal" ,"normal" ]
muX =[1.56, 47.7, 3.5, 8.0, 13.2, 0.59, 326.25]
covX =[0.036, 0.15, 0.21, 0.045, 0.086, 0.05, 0.17]
sigmmaX = list(np.array(covX)*np.array(muX))
n=7
rs =Gmulti(n,muX,sigmmaX,dist)
rs.RFn()
a1=rs.GetBeta()
dist[6]="lognormal"
rs=Gmulti(n,muX,sigmmaX,dist)
rs.RFn()
a2=rs.GetBeta()
dist=["lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"normal"]
rs=Gmulti(n,muX,sigmmaX,dist)
rs.RFn()
a3=rs.GetBeta()
dist=["lognormal", "lognormal", "lognormal", "lognormal" ,"lognormal" ,"lognormal" ,"lognormal"]
rs=Gmulti(n,muX,sigmmaX,dist)
rs.RFn()
a4=rs.GetBeta()
[a1,a2,a3,a4] 
#[3.8330281658128675, 3.761253665624187, 4.387684391077285, 4.090647398473588]
```

### Define a class for the limit state function, with input data provided in dictionary format.
```python
class GRS(ls.LSFM):  #Inherit from the LSFM class
    def __init__(self,n,muX,sigmmaX,dist):
        super().__init__(n,muX,sigmmaX,dist)
    def gcalc(self,X): #Virtural function for calculating g value
        R=X[0]
        S=X[1]
        g=R-S  #Definition of limit statue function
        return g  #Make sure to  return the g value.
    def dGdXcalc(self,X): #Virtural function for calculating the gradient vector of g
        R=X[0]
        S=X[1]
        dGdX=[0]*len(X) #Create a list of the same size as X.
        dGdX[0] =1.0
        dGdX[1] =-1.0
        return dGdX #Make sure to  return the gradient vector of g.
# main program #
#Representation of input data in dictionary format
data={"r":{"mean":200,"cov":0.05,"dist":"normal"},
"s":{"mean":100,"cov":0.2,"dist":"normal"}            }
n,muX,sigmmaX,dist=ls.dict2data(data)#Conversion of dictionary data into input data
rs=GRS(n,muX,sigmmaX,dist)
print('val=',rs.gcalc(muX))#Limit state function value at the initial point
rs.RFn()
print('beta=',rs.GetBeta())#Reliability index
print('alpha=',rs.GetAlpha())#Sensitivity
print('POF=',rs.GetPOF())#Probability of failure
#val= 100
#beta= 4.47213595499958
#alpha= [ 0.4472136  -0.89442719]
#POF= 3.872108215522035e-06 
```

### Vary the initial values over a wide range to search for the minimum reliability index.
```python
class GRS(ls.LSFM):  #Inherit from the LSFM class
    def __init__(self,n,muX,sigmmaX,dist):
        super().__init__(n,muX,sigmmaX,dist)
    def gcalc(self,X): #Virtural function for calculating g value
        R=X[0]
        S=X[1]
        g=R-S  #Definition of limit statue function
        return g  #Make sure to  return the g value.
    def dGdXcalc(self,X): #Virtural function for calculating the gradient vector of g
        R=X[0]
        S=X[1]
        dGdX=[0]*len(X) #Create a list of the same size as X.
        dGdX[0] =1.0
        dGdX[1] =-1.0
        return dGdX #Make sure to  return the gradient vector of g.
# main program #
#Representation of input data in dictionary format
data={"r":{"mean":200,"cov":0.05,"dist":"normal"},
"s":{"mean":100,"cov":0.2,"dist":"normal"}            }
n,muX,sigmmaX,dist=ls.dict2data(data)#辞書型データの入力データへの変換
import numpy as np
import optuna
    
class Objective:
    def __init__(self,data):
        # Initialization of variables
        self.data=data
    def __call__(self, trial):
        n,muX,sigmmaX,dist=ls.dict2data(self.data)
        rs=GRS(n,muX,sigmmaX,dist)
        ### "Specification of the search range for the starting point
        r = trial.suggest_uniform('r', 50,300)
        s = trial.suggest_uniform('s', 50,300)
        start=[r,s]
        rs.RFn(start='Coordinate',Xstart=start)
        return rs.GetBeta()  #Return the reliability index
# main program #
#Representation of input data in dictionary format
data={"r":{"mean":200,"cov":0.05,"dist":"normal"},
"s":{"mean":100,"cov":0.2,"dist":"normal"}            }
optuna.logging.disable_default_handler() #Suppress the output of Optuna
objective = Objective(data)
study = optuna.create_study(direction='minimize')
study.optimize(objective, timeout=60)
print('params:', study.best_params,',best_value:',study.best_value) 
#params: {'r': 296.52167806771865, 's': 124.24549136860253} ,best_value: 4.47213595499958
```



## methods of LimitState.LimitState
|method|contents|
|------|---------|
|RFn()|Design point search using the Rackwitz-Fiessler method|
|GetBeta()|Acquisition of the reliability index|
|GetAlpha()|Acquisition of sensitivity vector|
|GetPOF()|Acquisition of probability of failure value|
|GetDP()|Acquisition of design point|
|GetConv()|Acquisition of the number of conversion|
|GetPSF()|Acquisition of partial safety factor|
