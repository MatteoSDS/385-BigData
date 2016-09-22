negloglikelihood<-function(m,y,X,beta){
  total<-0
  N<-length(y)
  for(i in 1:N) total<-total+(m[i]-y[i])*t(X[i,])%*%beta+m[i]*log(1+exp(-t(X[i,])%*%beta))
  return(total)
}

unit_negloglikelihood<-function(m,y,X,beta,unit){
  result<-(m[unit]-y[unit])*t(X[unit,])%*%beta+m[unit]*log(1+exp(-t(X[unit,])%*%beta))
  return(result)
}

unit_gradient_negloglik<-function(m,y,X,beta,unit){
  unitX<-X[unit,]
  unitY<-y[unit]
  unitM<-m[unit]
  unitW<-1/(1+exp(-crossprod(unitX,beta)))
  S<-unitM*unitW-unitY
  grad<-unitX*S
  return(grad)
}

stocastic_gradient_descent<-function(m,y,X,beta0,stepsize,maxstepnumber){
  n<-dim(X)[1]
  unit<-sample(seq(1,n),1)
  betas<-matrix(NA,nrow=maxstepnumber,ncol=dim(X)[2])
  betas[1,]<-beta0
  negloglik<-numeric(maxstepnumber)
  negloglik[1]<-negloglikelihood(m,y,X,beta0)
  unit_negloglik<-numeric(maxstepnumber)
  unit_negloglik[1]<-negloglik[1]/n
  gradient<-unit_gradient_negloglik(m,y,X,beta0,unit)
  i<-1
  while(!(i==maxstepnumber)){
    i<-i+1
    beta1<-beta0-stepsize*gradient
    beta0<-beta1
    negloglik[i]<-negloglikelihood(m,y,X,beta0)
    unit_negloglik[i]<-unit_negloglikelihood(m,y,X,beta0,unit)
    betas[i,]<-beta0
    unit<-sample(seq(1,n),1)
    gradient<-unit_gradient_negloglik(m,y,X,beta0,unit)
  }
  return(list(betahat=beta0,negloglik=negloglik,
              unit_negloglik=unit_negloglik,betas=betas,step=i))
}

data_wdbc<-read.csv("./wdbc.csv", header=FALSE)
X<-as.matrix(cbind(rep(1,569),scale(data_wdbc[,3:12])))
y<-data_wdbc[,2]
y<-as.numeric(y=="M")
m<-rep(1,569)

beta0<-rep(0,11)
stepsize<-0.02
maxstepnumber<-100000
res_stoc_grad_desc<-stocastic_gradient_descent(m,y,X,beta0,stepsize,maxstepnumber)

plot(res_stoc_grad_desc$negloglik[1:res_stoc_grad_desc$step],
     main = "Negative Log-likelihood at each iteration - Fixed step",
     xlab="Iteration",ylab="negloglik",type="l",log="x")

running_avg_fix<-numeric(maxstepnumber)
running_avg_fix[1]<-res_stoc_grad_desc$unit_negloglik[1]
for(i in 2:maxstepnumber){
  running_avg_fix[i]<-(running_avg_fix[i-1]*(i-1)+res_stoc_grad_desc$unit_negloglik[i])/(i)}
plot(569*running_avg_fix[1:res_stoc_grad_desc$step],
     main = "Running average of the Unit Negative Log-likelihood at each iteration - Fixed step",
     xlab="Iteration",ylab="Avg unit negloglik",type="l",log="x")

par(mfrow=c(3,4))
estimates<-glm(y~0+X,family='binomial')$coefficients
for(i in 1:11){
plot(res_stoc_grad_desc$betas[1:res_stoc_grad_desc$step,i],
     main = paste("Beta",i,"by iteration - Fixed step"),
     xlab="Iteration",ylab="negloglik",type="l")
abline(h=estimates[i])}
par(mfrow=c(1,1))

robbins_monro_step<-function(iter,C,a){
  weight<-C*(iter+1)^-a
  return(weight)
}

sgd_decaying_step<-function(m,y,X,beta0,C,a,maxstepnumber){
  n<-dim(X)[1]
  unit<-sample(seq(1,n),1)
  betas<-matrix(NA,nrow=maxstepnumber,ncol=dim(X)[2])
  betas[1,]<-beta0
  negloglik<-numeric(maxstepnumber)
  negloglik[1]<-negloglikelihood(m,y,X,beta0)
  unit_negloglik<-numeric(maxstepnumber)
  unit_negloglik[1]<-negloglik[1]/n
  gradient<-unit_gradient_negloglik(m,y,X,beta0,unit)
  i<-1
  while(!(i==maxstepnumber)){
    i<-i+1
    beta1<-beta0-robbins_monro_step(i,C,a)*gradient
    beta0<-beta1
    negloglik[i]<-negloglikelihood(m,y,X,beta0)
    unit_negloglik[i]<-unit_negloglikelihood(m,y,X,beta0,unit)
    betas[i,]<-beta0
    unit<-sample(seq(1,n),1)
    gradient<-unit_gradient_negloglik(m,y,X,beta0,unit)
  }
  return(list(betahat=beta0,negloglik=negloglik,
              unit_negloglik=unit_negloglik,betas=betas,step=i))
}

C<-5
a<-0.8
beta0<-rep(0,11)
maxstepnumber<-100000
res_sgd_decaying_step<-sgd_decaying_step(m,y,X,beta0,C,a,maxstepnumber)

plot(res_sgd_decaying_step$negloglik[1:res_sgd_decaying_step$step],
     main = "Negative Log-likelihood at each iteration - Decaying step",
     xlab="Iteration",ylab="negloglik",type="l",log="x")

running_avg_decay<-numeric(maxstepnumber)
running_avg_decay[1]<-res_sgd_decaying_step$unit_negloglik[1]
for(i in 2:maxstepnumber){
  running_avg_decay[i]<-(running_avg_decay[i-1]*(i-1)+res_sgd_decaying_step$unit_negloglik[i])/(i)}
plot(569*running_avg_decay[1:res_sgd_decaying_step$step],
     main = "Running average of the Unit Negative Log-likelihood at each iteration - Decaying step",
     xlab="Iteration",ylab="Avg unit negloglik",type="l",log="x")

par(mfrow=c(3,4))
estimates<-glm(y~0+X,family='binomial')$coefficients
for(i in 1:11){
  plot(res_sgd_decaying_step$betas[1:res_sgd_decaying_step$step,i],
       main = paste("Beta",i,"by iteration - Decaying step"),
       xlab="Iteration",ylab="negloglik",type="l")
  abline(h=estimates[i])}
par(mfrow=c(1,1))

running_avg_beta_fix<-matrix(NA,nrow=maxstepnumber,ncol=11)
running_avg_beta_decay<-matrix(NA,nrow=maxstepnumber,ncol=11)
running_avg_beta_fix[1,]<-res_stoc_grad_desc$betas[1,]
running_avg_beta_decay[1,]<-res_sgd_decaying_step$betas[1,]
estimates<-glm(y~0+X,family='binomial')$coefficients
for(i in 2:maxstepnumber){
  running_avg_beta_fix[i,]<-(running_avg_beta_fix[i-1,]*(i-1)+res_stoc_grad_desc$betas[i,])/(i)
  running_avg_beta_decay[i,]<-(running_avg_beta_decay[i-1,]*(i-1)+res_sgd_decaying_step$betas[i,])/(i)}
running_negloglik_fix<-numeric(maxstepnumber)
running_negloglik_decay<-numeric(maxstepnumber)
running_negloglik_fix[1]<-negloglikelihood(m,y,X,running_avg_beta_fix[1,])
running_negloglik_decay[1]<-negloglikelihood(m,y,X,running_avg_beta_decay[1,])
for(i in 2:maxstepnumber){
  running_negloglik_fix[i]<-(running_negloglik_fix[i-1]*(i-1)+negloglikelihood(m,y,X,running_avg_beta_fix[i,]))/(i)
  running_negloglik_decay[i]<-(running_negloglik_decay[i-1]*(i-1)+negloglikelihood(m,y,X,running_avg_beta_decay[i,]))/(i)}
par(mfrow=c(1,2))
plot(running_negloglik_fix[1:res_stoc_grad_desc$step],
     main = "NegLogLik for averaged Betas - Fixed step",
     xlab="Iteration",ylab="negloglik",type="l",log="x")
plot(running_negloglik_decay[1:res_sgd_decaying_step$step],
     main = "NegLogLik for averaged Betas - Decaying step",
     xlab="Iteration",ylab="negloglik",type="l",log="x")
par(mfrow=c(1,1))