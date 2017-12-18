## (Author : Maïlis Amico - KU Leuven - mailis.amico@kuleuven.be) ##
## (last update : December 18, 2017) ##
## Estimation of a mixture cure model with a single-index model for the incidence ##

library(np) 
#library(nloptr)
setwd('/Volumes/Mailis/WORK/SuperComputer/Simulations-SIC') # A changer en fonction de l'emplacement de mes codes
#source("npindexbw_code.R",local=TRUE)
source("np_util.R",local=TRUE)

smsurv <-
  function(Time,Status,X,beta,w){    
    death_point <- sort(unique(subset(Time, Status==1)))
    coxexp <- exp((beta)%*%t(X[,-1]))  
    lambda <- numeric()
    event <- numeric()
    for(i in 1: length(death_point)){
      event[i] <- sum(Status*as.numeric(Time==death_point[i]))
      temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
      temp1 <- event[i]
      lambda[i] <- temp1/temp
    }
    HHazard <- numeric()
    for(i in 1:length(Time)){
      HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
      if(Time[i]>max(death_point))HHazard[i] <- Inf
      if(Time[i]<min(death_point))HHazard[i] <- 0
    }
    survival <- exp(-HHazard)
    list(survival=survival)
  }

#######
SIC <- function(n,Y,Delta,beta.init,gamma.init,X,Z,rescale,LCMM,eps,emmax,dataset){ 
  # Init values
  iteration = NULL
  m <- 0
  difference <- 1000
  k <- 1
  beta.hat <- coxph(Surv(Y,Delta)~Z,subset = Delta!=0,data = dataset,method = "breslow")$coef 
  B <- event
  gamma.hat <- glm(Delta~X,family=binomial(link=logit),data=dataset)$coefficients
  
  if(LCMM$b[2]>0){SIC1 <- 1}else{SIC1 <- -1}
  T_max <- max(Y[Delta==1])
  Z.fit <- cbind(rep(1,n),Z)
  s <- smsurv(Y,Delta,Z.fit,beta.hat,Delta)$survival # Baseline conditional survival function
  S_u <- drop(s^(exp((beta.init)%*%t(Z.fit[,-1]))))
  X.fit <- cbind(rep(1,n),X)
  p_X <- matrix(exp((gamma.init)%*%t(X.fit))/(1+exp((gamma.init)%*%t(X.fit))),ncol=1)
  
  gamma.hat <- gamma.init[-1]
  beta.hat <- beta.init
  
  while (difference>eps & k<emmax){
    # E-step
    w <- rep(NA,n)
    for(i in 1:n){
      if(Delta[i]==1) w[i] <- Delta[i] 
      if(Delta[i]==0) w[i] <- ((p_X[i]*S_u[i])/((1-p_X[i]) + p_X[i]*S_u[i]))
    }
    w[Delta==0 & Y>T_max] <- 0 # zero-tail constraint
    
    # M-step of the algorithm 
    # INCIDENCE 
    B <- w
    # Bandwidth 
    CVh <- function(param){                                   
      h <- param[1]
      h.comp <- function(h){
        index <- X%*%(gamma.hat)
        W <- as.matrix(data.frame(w,1))
        K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=h,ckertype="epanechnikov",leave.one.out=T)$ksum
        p2 <- K.sum[1,2,]/NZD(K.sum[2,2,])
        
        Floor <- sqrt(.Machine$double.eps) 
        p2[which(p2<Floor)] <- Floor 
        p2[which(p2>1-Floor)] <- 1 - Floor 
        contrib <- rep(NA,n)
        for(i in 1:dim(X)[1]){
          if(B[i]==0) contrib[i] <- log(1-p2[i])
          if(B[i]>0 & B[i]<1) contrib[i] <- B[i]*log(p2[i])+(1-B[i])*log(1-p2[i])
          if(B[i]==1) contrib[i] <- log(p2[i])
        } 
        CV <- -sum(contrib)
        
        return(CV)  
      }
      if(h>0){return(h.comp(h))}else{return(sqrt(.Machine$double.xmax))}
    }  
    
    if(m==0){h.init <- 0.5}else{h.init <- h}
    h.CV <- optim(par=h.init,CVh,method='Brent',lower=0.4,upper=1,control=list(fnscale=1)) 
    h <- h.CV$par
    bdwth <- h 
    
    l_I <- function(gamma2){    
      Floor <- sqrt(.Machine$double.eps)
      index <- X%*%c(SIC1,gamma2)
      W <- as.matrix(data.frame(B,1))
      K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=bdwth,ckertype="epanechnikov",leave.one.out=T)$ksum
      p2 <- K.sum[1,2,]/NZD(K.sum[2,2,])
      p2[which(p2<Floor)] <- Floor
      p2[which(p2>1-Floor)] <- 1 - Floor
      contrib <- rep(NA,n)
      for(i in 1:n){
        if(B[i]==0) contrib[i] <- log(1-p2[i])
        if(B[i]>0 & B[i]<1) contrib[i] <- B[i]*log(p2[i])+(1-B[i])*log(1-p2[i])
        if(B[i]==1) contrib[i] <- log(p2[i])
      }  
      likelihood <- -mean(contrib)     
    }
    
    starting <- gamma.hat[2:length(gamma.hat)] 
    incidence_est <- optim(par=starting,l_I,method='BFGS',control=list(fnscale=1)) 
    
    # P_X computation
    index <- X%*%c(SIC1,incidence_est$par)
    W <- as.matrix(data.frame(B,1))
    K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=bdwth,ckertype="epanechnikov",leave.one.out=F)$ksum
    p_X <- K.sum[1,2,]/NZD(K.sum[2,2,])
    
    
    # LATENCY 
    beta.update <- coxph(Surv(Y, Delta)~Z+offset(log(w)), subset=w!=0, method="breslow")$coef
    s.update <- smsurv(Y,Delta,Z.fit,beta.update,w)$survival 
    S_u <- drop(s.update^(exp((beta.update)%*%t(Z.fit[,-1]))))
    
    difference <- sum(c(incidence_est$par-gamma.hat[2:length(gamma.hat)],beta.update-beta.hat)^2)+sum((s-s.update)^2) 
    s <- s.update
    beta.hat <- beta.update
    gamma.hat[1] <- SIC1 
    gamma.hat[2:length(gamma.hat)] <- incidence_est$par
    iteration <- c(m,gamma.hat,beta.hat,bdwth,difference)
    m = m+1 
    if(m>emmax) break()
  }   
  
  # RESCALING
  if(rescale==TRUE){ 
    euclnorm.si <- sqrt(sum(gamma.hat^2))  
    gamma.new <- gamma.hat*1/euclnorm.si 
    h.new <- bdwth*(1/euclnorm.si) 
    # Predictions
    index <- X%*%gamma.new
    W <- as.matrix(data.frame(B,1))
    K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=h.new,ckertype="epanechnikov",leave.one.out=F)$ksum
    p_X <- K.sum[1,2,]/NZD(K.sum[2,2,])
  }else{gamma.new <- gamma.hat
  h.new <- bdwth 
  # Predictions
  index <- index
  p_X <- p_X
  } 
  
  cat("Cure probability model: \n")
  cat('        Estimate \n')
  print(rbind(as.data.frame(gamma.new),h.new))
  cat("\n")
  
  cat("Failure time distribution: \n")
  cat('        Estimate \n')
  print(as.data.frame(beta.hat))
  
  return(list(b=gamma.new,h=h.new,beta=beta.hat,interation=iteration,survival=s,cured=w,pred=cbind(index,p_X),h.conv=h.CV$convergence,conv=incidence_est$convergence)) 
  
}