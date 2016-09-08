library(Matrix)
library(microbenchmark)

inversion_m<-function(matA,matB){
  #matA*beta=matB
  beta<-solve.default(matA)%*%matB
  return(beta)
}

factorization_m<-function(matA,matB){
  #matA*beta=matB
  #matL'*matL*beta=matB
  #matL'*matC=matB
  #matL*beta=matC
  matL<-chol.default(matA)
  matC<-forwardsolve(t(matL),matB)
  beta<-backsolve(matL,matC)
  return(beta)
}

silly_data<-function(N,P){
  X<-matrix(rnorm(N*P,sd=10), nrow=N, ncol=P)
  truebeta<-rnorm(P,sd=10)
  W<-diag(sample(c(1,2),N,replace=T))
  error<-rnorm(N)
  Y = X%*%truebeta+error
  matA<-t(X)%*%W%*%X
  matB<-t(X)%*%W%*%Y
  return(list(matA=matA,matB=matB))
}

data<-silly_data(N=21,P=7)
matA<-data$matA;matB<-data$matB
p7<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),times=10);p7

data<-silly_data(N=150,P=50)
matA<-data$matA;matB<-data$matB
p50<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),times=10);p50

data<-silly_data(N=900,P=300)
matA<-data$matA;matB<-data$matB
p300<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),times=10);p300

data<-silly_data(N=3000,P=1000)
matA<-data$matA;matB<-data$matB
p1000<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),times=10);p1000

sparse_data<-function(N,P,prop=0.05){
  X<-matrix(rnorm(N*P,sd=10), nrow=N, ncol=P)
  mask = matrix(rbinom(N*P,1,prop), nrow=N)
  X = mask*X
  truebeta<-rnorm(P,sd=10)
  W<-diag(sample(c(1,2),N,replace=T))
  error<-rnorm(N)
  Y = X%*%truebeta+error
  matA<-Matrix(t(X)%*%W%*%X,sparse=T)
  matB<-Matrix(t(X)%*%W%*%Y,sparse=F)
  return(list(matA=matA,matB=matB))
}

sparse_m<-function(matA,matB){
  #matA*beta=matB
  beta<-solve(matA,matB)
  return(beta)
}

data<-sparse_data(N=900,P=300,prop=0.05)
matA<-data$matA;matB<-data$matB
p300_0.05<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),
                          sparse_m(matA,matB),times=10);p300_0.05

data<-sparse_data(N=3000,P=1000,prop=0.05)
matA<-data$matA;matB<-data$matB
p1000_0.05<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),
                           sparse_m(matA,matB),times=10);p1000_0.05

data<-sparse_data(N=3000,P=1000,prop=0.01)
matA<-data$matA;matB<-data$matB
p1000_0.01<-microbenchmark(inversion_m(matA,matB),factorization_m(matA,matB),
                           sparse_m(matA,matB),times=10);p1000_0.01
