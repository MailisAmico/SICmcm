## Maïlis Amico - KU Leuven - mailis.amico@kuleuven.be
## (last update : January 4, 2018) 
## The single-index/Cox mixture cure model examples

source('Functions-SIC.R')
library(smcure)
library(np)
data('e1684')
e1684 <- subset(e1684,e1684$AGE!='') # remove observation for which there is missing data (only one)
attach(e1684)

Z <- cbind(AGE,TRT,SEX)
LCmm.estim <- smcure(Surv(FAILTIME,FAILCENS)~AGE+TRT+SEX,cureform=~AGE+TRT+SEX,model="ph",data=e1684,Var=FALSE)
gamma.init <- glm(FAILCENS~Z,family=binomial(link=logit),data=e1684)$coefficients
beta.init <- coxph(Surv(FAILTIME,FAILCENS)~Z,subset = FAILCENS!=0,data = e1684,method = "breslow")$coef 

# EM algorithm
SICmm.estim1 <- SIC(Time=FAILTIME,Status=FAILCENS,gamma.init,beta.init,Z,Z,LCMM=LCmm.estim,eps=1e-5,emmax=20,dataset=e1684,rescale=TRUE)
  
# Estimation + bootstrap
SICmm.estim2 <- SICMM(Time=FAILTIME,Status=FAILCENS,Z,Z,LCMM=LCmm.estim,eps=1e-5,emmax=20,dataset=e1684,rescale=TRUE,bootstrap=TRUE,nboot=250)
