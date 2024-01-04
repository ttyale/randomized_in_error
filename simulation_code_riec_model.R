require("lme4")
require("LaplacesDemon")
require("mvtnorm")
require("truncnorm")

### sample from the Wishart distribution
# rwish<-function(n,nu0,S0)
# {
#   sS0 <- chol(S0)
#   S<-array( dim=c( dim(S0),n ) )
#   for(i in 1:n)
#   {
#     Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
#     S[,,i]<- t(Z)%*%Z
#   }
#   S[,,1:n]
# }

gendata<-function(N=500,beta0=3.1){
  ##study parameters
  #N=500   #sample size  
  #err=0.1 #randomization error rate
  
  trt<-rep(c(1,0),each=N/2)

  ##assume two covariates
  x1<-rnorm(N,0,2) #normal continuous 
  x2<-runif(N,-5,5)  #uniform continuous 
  X<-as.matrix(data.frame(rep(1,N),x1,x2))   #covariate matrix
  
  ##strata model parameters
  #phi2<-1 #cluster variance
  #generate latent variable for G group status
  beta=c(beta0,1.2,0.5)  #beta0=3.1 is 15% 
  S = rbinom(N, 1, pnorm(X%*%beta))
  
  #outcome model parameters
  alpha0<-c(-1.5,0.9,0.5) 
  alpha1<-c(1.5,0.5,0.8)   
  sigma2<-5 #outcome variance
  resid<-rnorm(N,0,sqrt(sigma2))
  
  #unobserved model for outcome of patients randomized in error
  gamma0 = c(-5,-3, -4)
  gamma1 = c(+5, 2, 3)
  tau2<-10
  resid2<-rnorm(N,0,sqrt(tau2))
  
  
  Y1<-rep(NA,N) #outcome 
  Y1[intersect(which(S==1),which(trt==1))]<-X[intersect(which(S==1),which(trt==1)),]%*%alpha1+resid[intersect(which(S==1),which(trt==1))]
  Y1[intersect(which(S==1),which(trt==0))]<-X[intersect(which(S==1),which(trt==0)),]%*%alpha0+resid[intersect(which(S==1),which(trt==0))]
  Y1[intersect(which(S==0),which(trt==0))]<-X[intersect(which(S==0),which(trt==0)),]%*%gamma0+resid2[intersect(which(S==0),which(trt==0))]
  

  Y2<-rep(NA,N) #outcome with rie Y's measured in treatment
  Y2[intersect(which(S==1),which(trt==1))]<-X[intersect(which(S==1),which(trt==1)),]%*%alpha1+resid[intersect(which(S==1),which(trt==1))]
  Y2[intersect(which(S==1),which(trt==0))]<-X[intersect(which(S==1),which(trt==0)),]%*%alpha0+resid[intersect(which(S==1),which(trt==0))]
  Y2[intersect(which(S==0),which(trt==0))]<-X[intersect(which(S==0),which(trt==0)),]%*%gamma0+resid2[intersect(which(S==0),which(trt==0))]
  Y2[intersect(which(S==0),which(trt==1))]<-X[intersect(which(S==0),which(trt==1)),]%*%gamma1+resid2[intersect(which(S==0),which(trt==1))]
  
  
  
  #true TE 
  TE<-colMeans(X[which(S==1),]%*%alpha1-X[which(S==1),]%*%alpha0)
  
  #set S in control as unknown
  S0<-S
  S[trt==0]<-NA
  
  #For G, the true status of 0, 1 or 2 will not be known, and death in treatment group and survival in control group have known G status   
  return(data.frame(cbind(Y1,Y2,x1,x2,trt,S,S0,TE)))
  
}



#error rate
set.seed(245798)
err<-te<-lmest<-rep(NA, 1000)

for (i in 1:1000){
 
  df<-gendata(N=100, beta0=3.1)
  err[i]<-sum(df$S0)/length(df$S0)
   te[i]<-df$TE[1]
   dfX<-df[-which(df$S==0),]
   lmest[i]<-summary(lm(Y1~trt, data=dfX))$coef[2,1]
}
mean(err) 
mean(te)
mean(lmest)

#set.seed(87345313)
#err<-te<-rep(NA, 1000)
#for (i in 1:1000){
#  df<-gendata(N=500, beta0=1.6)
#  df<-gendata(N=100, beta0=1.6)
#  err[i]<-sum(df$S0)/length(df$S0)
#   te[i]<-df$TE[1]
#}
#mean(err)
#mean(te)

#df<-gendata(N=500)

#Gibbs sampler with probit model when outcome in trt among rie is not measured 
rie_sampler<-function(df=df, ITER=10000, beta0=3.1){
  
  #df<-gendata(N=1500,s=60)
  #model information
  N=nrow(df)
  p=3
  X<-as.matrix(cbind(rep(1,N),df$x1,df$x2))
  
  #if (ymeasure == F) 
  Y= df$Y1
  
  D<-df$trt
  
  #outcome model parameters
 
  alpha1<-c(1.5,0.5,0.8)   
  alpha0<-c(-1.5,0.9,0.5) 
  gamma0<-c(-5,-3, -4)
  
  sigma2<-5
  tau2<-10
  #initial value for missing S
  S=df$S
  S[D==0]<-rbinom(length(which(D==0)), 1, 0.9)
  beta=matrix(ncol=3, nrow=1,c(beta0,1.2,0.5)) #group00 vs group11
  #beta=matrix(ncol=3, nrow=1,c(1.6,1.2,0.5)) #group00 vs group11
  
  ##priors 
  #outcome model
  #regression coefficients
  #a111<-a101<-a110<-rep(0,p) set all prior means of regression parameters as 0
  Sigma1<-Sigma0<-Gamma0<-Lambda<-diag(p)*1000

  cc=dd=ee=ff=0.001

  #strata model
  #regression coefficients
  beta0 = rep(0,p)
  
  ##store chain S
  ALPHA1<-ALPHA0<-GAMMA0<-matrix(NA,p,ITER)
  SIGMA<-TE<-TAU<-rep(0,ITER) #sample TE by estimating the potential outcomes of G=2.
  
  BETA<-matrix(NA,p,ITER)
  GV<-matrix(NA,N,ITER)

  #store the random intercept of latent strata model
  
  for(i in 1:ITER){
    
    #update alpha
    tgt<-intersect(which(D==1),which(S==1)) 
    v=solve(crossprod(X[tgt,])/sigma2+solve(Sigma1))
    m=v%*%(t(X[tgt,])%*%Y[tgt]/sigma2)
    alpha1=rmvnorm(1,m,v)
    ALPHA1[,i]<-alpha1
    
    tgt<-intersect(which(D==0),which(S==1)) 
    v=solve(crossprod(X[tgt,])/sigma2+solve(Sigma0))
    m=v%*%(t(X[tgt,])%*%Y[tgt]/sigma2)
    alpha0=rmvnorm(1,m,v)
    ALPHA0[,i]<-alpha0
    
    tgt<-intersect(which(D==0),which(S==0)) 
    if(length(tgt)==1){
      v=solve(tcrossprod(X[tgt,])/tau2+solve(Gamma0))
      m=v%*%((X[tgt,])*(Y[tgt])/tau2)
    }else{
      v=solve(crossprod(X[tgt,])/tau2+solve(Gamma0))
      m=v%*%(t(X[tgt,])%*%Y[tgt]/tau2)
    }
    gamma0=rmvnorm(1,m,v)
    GAMMA0[,i]<-gamma0
    
    #update sigma2
    tgt1<-intersect(which(S==1),which(D==1))
    tgt2<-intersect(which(S==1),which(D==0))
  
    rate<-dd+0.5*(sum((Y[tgt1]-X[tgt1,]%*%t(alpha1))^2))+
             0.5*(sum((Y[tgt2]-X[tgt2,]%*%t(alpha0))^2))
    
    sigma2<-rgamma(1,shape=(cc+sum(S,na.rm=T)/2),rate=rate)^(-1)
    SIGMA[i]<-sigma2
    
    #update tau2
    tgt1<-intersect(which(S==0),which(D==0))
    #tgt2<-intersect(which(S==1),which(D==0))
    
    rate<-ff+0.5*(sum((Y[tgt1]-X[tgt1,]%*%t(gamma0))^2))
    tau2<-rgamma(1,shape=(ee+sum((1-S)*(1-D),na.rm=T)/2),rate=rate)^(-1)
    TAU[i]<-tau2
    
    
    #potential outcomes
    tgt1<-intersect(which(S==1),which(D==1))
    tgt2<-intersect(which(S==1),which(D==0))
    dff<-X[c(tgt1,tgt2),]
    #X1<-as.matrix(cbind(rep(1,nrow(df1)),df1$x1,df1$x2))
    
    Y1<-dff%*%t(alpha1)
    Y0<-dff%*%t(alpha0)
    te<-mean(Y1-Y0)
    TE[i]<-te 
    
    ######################################
    #### mixed membership model
 
    Z<-rep(NA,N)
    zmean<-X%*%t(beta)
    
    for (j in 1:N){
        if (S[j] == 1) Z[j] = rtruncnorm(1,a=0,mean=as.numeric(zmean[j]), sd =1) else Z[j] = rtruncnorm(1,b=0,mean=as.numeric(zmean[j]), sd =1)
    }
      
    v=solve(crossprod(X)+solve(Lambda))
    m=v%*%(crossprod(X,Z)+solve(Lambda)%*%beta0)
    beta=rmvnorm(1,m,v)
    BETA[,i]<-beta
    
    #Update membership
    ################################# this part needs a fix. other parts look all correct#
    #membership for RIE in control
    #P01<-rep(NA,length(which(D==0)))  
    #P00<-rep(NA,length(which(D==0)))  
    
    P01<-dnorm(as.numeric(Y[D==0]-X[D==0,]%*%t(alpha0))/sqrt(sigma2))
    P00<-dnorm(as.numeric(Y[D==0]-X[D==0,]%*%t(gamma0))/sqrt(tau2))
   
    Pstrata<-pnorm(X[D==0,]%*%t(beta))
    P=Pstrata*P01/(Pstrata*P01+(1-Pstrata)*P00)     
    S[D==0]<-rbinom(length(which(D==0)),1,P)
    
    #store G
    GV[,i]<-S
    if (i%%500==0) print(i)
    #print(i)
    #GT[,i]<-G
    
  }  
  return(list(ALPHA1=ALPHA1,ALPHA0=ALPHA0, GAMMA0=GAMMA0,BETA=BETA,SIGMA=SIGMA, TAU=TAU, GV=GV,TE=TE))
}

#sima<-rie_sampler1(df=df, ITER=1000)



sim_function<-function(S=10000, nsim=500, N=100, beta0=3.1){
  p=3
  true_TE<-array(NA,c(nsim,1))
  true_GV<-array(NA,c(nsim,2))
  ALPHA1<-ALPHA0<-GAMMA0<-array(NA,c(nsim,p,S))
 
  GV<-array(NA,c(nsim,N,S))
  
  BETA<-array(NA,c(nsim,3,S))

  SIGMA<-TAU<-TE<-matrix(NA,nrow=nsim,ncol=S)
  
  TE0<-TE1<-TE2<-matrix(NA,nsim,3)
  
  for(i in 1:nsim){
    #set.seed(i*10)
    print(paste("replicate:",i))
    df<-gendata(N=N, beta0=beta0)
    true_TE[i,]<-c(df$TE[1])
    true_GV[i,1]<-sum(df$S0)/nrow(df)
    true_GV[i,2]<-1-sum(df$S0)/nrow(df)
    ########################################################################
    #otp1<-try(rie_sampler(df=df,ITER=S)) #_NP is the function for known strata 
    #if (class(otp1) =="try-error"){
    #  next
    #}
    
    otp1<-rie_sampler(df=df,ITER=S,beta0=beta0)
    
    #otp<-combined.Gibbs(df=df,S=S) 
    ALPHA1[i,,]<-otp1$ALPHA1
    ALPHA0[i,,]<-otp1$ALPHA0
    GAMMA0[i,,]<-otp1$GAMMA0
    
    BETA[i,,]<-otp1$BETA
    SIGMA[i,]<-otp1$SIGMA
    TAU[i,]<-otp1$TAU
    TE[i,]<-otp1$TE
    GV[i,,]<-otp1$GV
    
    #########################################################################  
    #frequentist
    dfX<-df[-which(df$S==0),]
    otp0<-try(lm(Y1~trt, data=dfX)) #_NP is the function for known strata 
    if (class(otp0) =="try-error"){
      next
    }
    
    otp1<-try(lm(Y1~trt+x1+x2, data=dfX))
    if (class(otp1) =="try-error"){
      next
    }
    
    df1<-df[complete.cases(df[,c("Y1","x1","x2","trt")]),]
    m<-glm(trt~x1+x2,data=df1,family=binomial(link=logit))
    p<-predict(m,type="response")
    w<-1/(p*I(df1$trt==1)+(1-p)*I(df1$trt==0))
    
    otp2<-try(lm(Y1~trt,data=df1, weights = w))
    if (class(otp1) =="try-error"){
      next
    }
    
    TE0[i,]<-c(coef(otp0)[2], coef(otp0)[2]-1.96*summary(otp0)$coef[2,2] , coef(otp0)[2]+1.96*summary(otp0)$coef[2,2])
    
    TE1[i,]<-c(coef(otp1)[2], coef(otp1)[2]-1.96*summary(otp1)$coef[2,2] , coef(otp1)[2]+1.96*summary(otp1)$coef[2,2])
    
    TE2[i,]<-c(coef(otp2)[2], coef(otp2)[2]-1.96*summary(otp2)$coef[2,2] , coef(otp2)[2]+1.96*summary(otp2)$coef[2,2])
  
  }
  
  output<-list(ALPHA1=ALPHA1,ALPHA0=ALPHA0, GAMMA0=GAMMA0,
               BETA=BETA,
               GV = GV,
               SIGMA=SIGMA,
               TAU=TAU,
               TE=TE,
               TE0=TE0,
               TE1=TE1,
               TE2=TE2,
               true_TE=true_TE,
               true_GV=true_GV)
  
  return(output)
}





#1
#N=250; error rate = 30%
set.seed(8695863)
sim1<-sim_function(S=10000, nsim=500, N=250, beta0=1.6)
saveRDS(sim1, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n250_erate_30_sim_001-500.RDS" )

set.seed(10137341)
sim2<-sim_function(S=10000, nsim=500, N=250, beta0=1.6)
saveRDS(sim2, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n250_erate_30_sim_501-1000.RDS" )

#2
#N=250; error rate = 15%
set.seed(57928354)
sim3<-sim_function(S=10000, nsim=500, N=250, beta0=3.1)
saveRDS(sim3, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n250_erate_15_sim_001-500.RDS" )

set.seed(8623342)
sim4<-sim_function(S=10000, nsim=500, N=250, beta0=3.1)
saveRDS(sim4, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n250_erate_15_sim_501-1000.RDS" )


#3
#N=500; error rate = 30%
set.seed(76839235)
sim5<-sim_function(S=10000, nsim=500, N=500, beta0=1.6)
saveRDS(sim5, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n500_erate_30_sim_001-500.RDS" )

set.seed(58173643)
sim6<-sim_function(S=10000, nsim=500, N=500, beta0=1.6)
saveRDS(sim6, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n500_erate_30_sim_501-1000.RDS" )

#4
#N=500; error rate = 15%
set.seed(768392234)
sim7<-sim_function(S=10000, nsim=500, N=500, beta0=3.1)
saveRDS(sim7, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n500_erate_15_sim_001-500.RDS" )

set.seed(58173589)
sim8<-sim_function(S=10000, nsim=500, N=500, beta0=3.1)
saveRDS(sim8, "/home/cl2492/palmer_scratch/SACE/RIE/RIE_miss/rie_simulation_n500_erate_15_sim_501-1000.RDS" )






################################################################################################################################
#simple linear regression


sim_function_lm<-function(nsim=500, N=100, beta0=3.1, seed=1234567){
  p=3
  set.seed(seed)

  TE0<-TE1<-TE2<-matrix(NA,nsim,3)

  for(i in 1:nsim){
    print(paste("replicate:",i))
    df<-gendata(N=N, beta0=beta0)
    
    ########################################################################
    otp0<-try(lm(Y1~trt, data=df)) #_NP is the function for known strata 
    if (class(otp0) =="try-error"){
      next
    }
    
    otp1<-try(lm(Y1~trt+x1+x2, data=df))
    if (class(otp1) =="try-error"){
      next
    }
  
    df1<-df[complete.cases(df[,c("Y1","x1","x2","trt")]),]
    m<-glm(trt~x1+x2,data=df1,family=binomial(link=logit))
    p<-predict(m,type="response")
    w<-1/(p*I(df1$trt==1)+(1-p)*I(df1$trt==0))
  
    otp2<-try(lm(Y1~trt,data=df1, weights = w))
    if (class(otp1) =="try-error"){
      next
    }
    
    TE0[i,]<-c(coef(otp0)[2], coef(otp0)[2]-1.96*summary(otp0)$coef[2,2] , coef(otp0)[2]+1.96*summary(otp0)$coef[2,2])
    
    TE1[i,]<-c(coef(otp1)[2], coef(otp1)[2]-1.96*summary(otp1)$coef[2,2] , coef(otp1)[2]+1.96*summary(otp1)$coef[2,2])
    
    TE2[i,]<-c(coef(otp2)[2], coef(otp2)[2]-1.96*summary(otp2)$coef[2,2] , coef(otp2)[2]+1.96*summary(otp2)$coef[2,2])
    
  }
  
  output<-list(TE0=TE0,
               TE1=TE1,
               TE2=TE2)
  
  return(output)
}


lmsim<-sim_function_lm(nsim=500, N=100, beta0=3.1, seed=20230409)

colMeans(lmsim$TE0)
colMeans(lmsim$TE1)
colMeans(lmsim$TE2)


lmsim<-sim_function_lm(nsim=500, N=500, beta0=3.1, seed=20230410)

colMeans(lmsim$TE0)
colMeans(lmsim$TE1)
colMeans(lmsim$TE2)


lmsim<-sim_function_lm(nsim=500, N=100, beta0=1.6, seed=20230411)

colMeans(lmsim$TE0)
colMeans(lmsim$TE1)
colMeans(lmsim$TE2)


lmsim<-sim_function_lm(nsim=500, N=500, beta0=1.6, seed=20230412)

colMeans(lmsim$TE0)
colMeans(lmsim$TE1)
colMeans(lmsim$TE2)

