require(Matrix)
require(readr)

#calculate matrix L
fill_mL<-function(rho,propJ,propK,size){
  cJ<-dim(propJ)[2]
  cK<-dim(propK)[2]
  cN<-length(size)
  #our matrix for the lagrangian is going to contain
  #cJ*cK parameters for each cN precint
  #cJ+cK multipliers for linear constraints for each cN precint
  #cJ*cK parameters for the mean point
  #cJ*cK multipliers for linear constraints for the mean point
  cJK<-cJ*cK
  mL_dim<-(cJK+cJ+cK)*cN+2*cJK
  mA_dim<-cJK*cN+cJK
  mL<-Matrix(0,nrow=mL_dim,ncol=mL_dim)
  #filling matrix A
  sumSize<-sum(size)
  for(i in 1:cN){
    mL[1:cJK+cJK*(i-1),1:cJK+cJK*(i-1)]<- size[i]*diag(cJK)
    mL[1:cJK+cJK*(i-1),1:cJK+cJK*cN]<- -size[i]*diag(cJK)
    mL[1:cJK+cJK*cN,1:cJK+cJK*(i-1)]<- -size[i]*diag(cJK)
  }
  mL[1:cJK+cJK*cN,1:cJK+cJK*cN]<- sumSize*diag(cJK)
  #filling linear constrains
  #sum b=1
  mTemp<-Matrix(0,nrow=cK,ncol=cJK)
  for(k in 1:cK){
    mTemp[k,]<-c(rep(0,k-1),rep(c(1,rep(0,cK-1)),cJ-1),1,rep(0,cK-k))}
  mTemp_t<-t(mTemp)
  for(i in 1:cN){
    mL[1:cK+mA_dim+cK*(i-1),1:cJK+cJK*(i-1)]<-mTemp
    mL[1:cJK+cJK*(i-1),1:cK+mA_dim+cK*(i-1)]<-mTemp_t
  }
  #y=sum b*x
  for(i in 1:cN){
    mTemp<-Matrix(0,nrow=cJ,ncol=cJK)
    for(j in 1:cJ){
      mTemp[j,]<-c(rep(0,cK*(j-1)),propK[i,],rep(0,cK*(cJ-j)))}
    mL[1:cJ+mA_dim+cK*cN+cK*(i-1),1:cJK+cJK*(i-1)]<-mTemp
    mL[1:cJK+cJK*(i-1),1:cJ+mA_dim+cK*cN+cK*(i-1)]<-t(mTemp)
  }
  #mean b
  mTemp<-Matrix(0,nrow=cJK,ncol=cJK*cN+cJK)
  for(i in 1:cN){
    mTemp[1:cJK,1:cJK+cJK*(i-1)]<-size[i]*diag(cJK)}
  mTemp[1:cJK,1:cJK+cJK*cN]<- -sumSize*diag(cJK)
  mL[1:cJK+mA_dim+cK*cN+cJ*cN,1:(cJK*cN+cJK)]<-mTemp
  mL[1:(cJK*cN+cJK),1:cJK+mA_dim+cK*cN+cJ*cN]<-t(mTemp)
  #extract D and sum it squared
  mD<-mL[(mA_dim+1):mL_dim,1:mA_dim]
  mL[1:mA_dim,1:mA_dim]<-mL[1:mA_dim,1:mA_dim]+rho*crossprod(mD)
  #completed matrix
  return(mL)
}

#calculate vector E
fill_vE<-function(rho,mL,propJ,propK,size){
  cJ<-dim(propJ)[2]
  cK<-dim(propK)[2]
  cN<-length(size)
  cJK<-cJ*cK
  mL_dim<-(cJK+cJ+cK)*cN+2*cJK
  mA_dim<-cJK*cN+cJK
  vC<-c(rep(1,cK*cN),as.vector(t(propJ)),rep(0,cJK))
  t_mD<-mL[1:mA_dim,(mA_dim+1):mL_dim]
  t_mD_vC<-as.vector(t_mD%*%vC)
  vE<-c(-rho*t_mD_vC,-vC)
  return(vE)
}

#just to add some verbose to the optimization
obj_fun<-function(vBetas,cJK){
  nBetas<-length(vBetas)
  bar_betas<-vBetas[(1+nBetas-cJK):nBetas]
  betas<-vBetas[1:(nBetas-cJK)]
  totale<-0
  for(i in 1:cJK) totale<-totale+sum((betas[seq(i,nBetas-cJK,by=cJK)]-bar_betas[i])^2)
  return(totale)
}

#main function for the projected gradient descent
proj_gradient_descent<-function(cJ,cK,cN,mL,vE,mu0,stepsize,maxstepnumber,eps){
  cJK<-cJ*cK
  mA_dim<-cJK*cN+cJK
  gradient<-mL%*%mu0+vE
  diff<-eps+1
  i<-1
  while(!(i==maxstepnumber)&&(eps<diff)){
    i<-i+1
    mu1<-mu0-stepsize*gradient
    mu1[mu1[1:(mA_dim)]<0]<-0
    mu1[mu1[1:(mA_dim)]>1]<-1
    diff<-sum(abs(mu0-mu1))
    mu0<-mu1
    gradient<-mL%*%mu0+vE
    if(i%%1000==0){
      print(c(i,obj_fun(mu0[1:mA_dim],cJK),sum(abs(gradient[1:mA_dim]))))
      print(base_betas)
      print(round(matrix(mu0[(cJK*cN+1):(cJK*cN+cJK)],nrow=cJ,ncol=cK,byrow=T),3))
    }
  }
  return(list(muhat=mu0,step=i))
}

#simulated data, small example
#we'll try to recover the base_betas
size<-rep(1,30)
propK<-matrix(runif(90,0,5),nrow=30,ncol=3)
propK<-propK/rowSums(propK)
base_betas<-matrix(c(2,3,5,2,6,2,4,3,3)/10,nrow=3,ncol=3)
propJ<-propK
for(i in 1:30){
  for(j in 1:3){
    noise<-runif(3,0,0.1)
    noise<-noise-mean(noise)
    propJ[i,j]<-crossprod(propK[i,],base_betas[,j]+noise) 
  }
}

#initialize and run
cJ<-dim(propJ)[2]
cK<-dim(propK)[2]
cN<-length(size)
cJK<-cJ*cK
rho<-10
mL<-fill_mL(rho,propJ,propK,size)
vE<-fill_vE(rho,mL,propJ,propK,size)
#uniform starting point
mu0<-c(rep(1/cJ,cJK*(cN+1)),rep(0.01,cK*cN),as.vector(t(propJ))-1/cJ,rep(0.01,cJK))
maxsteps<-100000
stepsize<-1/10000
eps<-10^-4
  
#runs 4 minutes without converging
data1<-date();data1
results<-proj_gradient_descent(cJ,cK,cN,mL,vE,mu0,stepsize,maxsteps,eps)
data2<-date();data1;data2
results$step
base_betas
round(matrix(results$muhat[(cJK*cN+1):(cJK*cN+cJK)],nrow=cJ,ncol=cK,byrow=T),3)
