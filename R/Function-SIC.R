## (Author : Maïlis Amico - KU Leuven - mailis.amico@kuleuven.be) 
## (last update : December 18, 2017) 
## SIC cure model functions

# smsurv : Estimation of the baseline conditional survival function (from smcure package)
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
# EM algorithm to estimate the single-index/Cox mixture cure model
SIC <- function(Y,Delta,beta.init,gamma.init,X,Z,LCMM,eps,emmax,dataset,rescale){ 
  # Init values
  iteration = NULL
  m <- 0
  difference <- 1000
  k <- 1
  n <- dim(dataset)[1]
  
  if(LCMM$b[2]>0){SIC1 <- 1}else{SIC1 <- -1}
  T_max <- max(Y[Delta==1])
  Z.fit <- cbind(rep(1,n),Z)
  s <- smsurv(Y,Delta,Z.fit,beta.init,Delta)$survival # Baseline conditional survival function
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
  
  # NORMALIZATION
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
  p_X <- p_X}
    
  return(list(b=gamma.new,h=h.new,beta=beta.hat,interation=iteration,survival=s,cured=w,pred=cbind(index,p_X),h.conv=h.CV$convergence,conv=incidence_est$convergence)) 
  
}


#######
SICMM <- function(Y,Delta,X,Z,LCMM,eps,emmax,dataset,rescale,bootstrap,nboot){
  n <- dim(dataset)[1]
  beta.init <- coxph(Surv(Y,Delta)~Z,subset = Delta!=0,data =dataset,method = "breslow")$coef 
  gamma.init <- glm(Delta~X,family=binomial(link=logit),data=dataset)$coefficients
  SIC.model <- SIC(Y,Delta,beta.init,gamma.init,X,Z,LCMM,eps,emmax,dataset,rescale)
  
  if(bootstrap==TRUE){
    bootstrap_gamma <- NULL
    bootstrap_beta <- NULL
    ngamma <- ncol(X)
    gammanm <- colnames(X)
    nbeta <- ncol(as.matrix(Z))
    betanm <- colnames(Z)
    gamma_boot <- matrix(rep(0,nboot*(ngamma)),nrow=nboot)
    beta_boot <- matrix(rep(0,nboot*(nbeta)),nrow=nboot)
    iter <- matrix(rep(0,nboot),ncol=1)
    tempdata <- cbind(time,event,X,Z)
    inb <- 1
    while(inb <= nboot){
      set.seed(inb)
      id <- sample(1:n,n,replace=TRUE) # resampling from dataset
      bootdata <- tempdata[id,] 
      bootX <- bootdata[,gammanm] # covariates for the incidence for bootstrap resample
      colnames(bootX) <- c("bootX1","bootX2","bootX3","bootX4")
      bootZ <- bootdata[,betanm] # covariates for the latency for bootstrap 
      
      bootdata <- data.frame(bootdata[,1],bootdata[,2],bootX,bootZ)
      bootfit <- SIC(Y=bootdata[,1],Delta=bootdata[,2],beta.init=beta.init,gamma.init=gamma.init,X=bootX,Z=bootZ,LCMM,eps,emmax,bootdata,rescale)
      gamma_boot[inb,] <- bootfit$b
      beta_boot[inb,] <- bootfit$beta
      
      inb <- inb+1
      print(inb)
      
    }
    gamma_se <- sqrt(apply(gamma_boot,2,var))
    gamma_zvalue <- SIC.model$b/gamma_se
    gamma_pvalue <- (1-pnorm(abs(gamma_zvalue)))*2
    beta_se <- sqrt(apply(beta_boot,2,var))
    beta_zvalue <- SIC.model$beta/beta_sd
    beta_pvalue <- (1-pnorm(abs(beta_zvalue)))*2
    
    return(list(b=SIC.model$b,h=SIC.model$h,beta=SIC.model$beta,survival=SIC.model$survival,pred=SIC.model$pred,h.conv=SIC.model$h.conv,conv=SIC.model$conv))
    
  }else{
    return(list(b=SIC.model$b,h=SIC.model$h,beta=SIC.model$beta,survival=SIC.model$survival,pred=SIC.model$pred,h.conv=SIC.model$h.conv,conv=SIC.model$conv))
  }
}



