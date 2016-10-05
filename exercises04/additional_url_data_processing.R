require(Matrix)
require(readr)

#Use original url_data_processing.R to obtain url_X.rds and url_y.rds

#then to scale X

rm(list=ls())
X<-readRDS("url_X.rds")
n <- nrow(X)
p<-ncol(X)
X_mean<-colMeans(X)
X@x<-(X@x)^2
X2_mean<-colMeans(X)
colVar <- X2_mean - (X_mean)^2
colSD<-sqrt(colVar * n/(n-1))
inv_colSD<-1/colSD
inv_colSD[inv_colSD==Inf]<-1
X<-readRDS("url_X.rds")
X<-X%*%Diagonal(p,inv_colSD)
saveRDS(X, file='url_X_scaled.rds')
rm(list=ls())

#then to split into training and test

rm(list=ls())
X<-readRDS("url_X_scaled.rds")
n<-nrow(X)
p<-ncol(X)
training_size<-round(n*0.8,0)
sampled_rows<-sample.int(n,size=training_size)
X_training<-X[sampled_rows,]
saveRDS(X_training, file='url_X_training.rds')
rm(X_training)
X_test<-X[-sampled_rows,]
saveRDS(X_test, file='url_X_test.rds')
rm(X_test)
rm(X)
y<-readRDS("url_y.rds")
y_training<-y[sampled_rows]
saveRDS(y_training, file='url_y_training.rds')
rm(y_training)
y_test<-y[-sampled_rows]
saveRDS(y_test, file='url_y_test.rds')
rm(list=ls())

#then to format by row

rm(list=ls())
X<-readRDS("url_X_training.rds")
X<-as(X,"RsparseMatrix")
X_dim<-dim(X)
saveRDS(X_dim, file='url_X_training_dim.rds')
rm(X_dim)
X_p<-X@p
saveRDS(X_p, file='url_X_training_p.rds')
rm(X_p)
X_j<-X@j
saveRDS(X_j, file='url_X_training_j.rds')
rm(X_j)
X_x<-X@x
saveRDS(X_x, file='url_X_training_x.rds')
rm(list=ls())
X<-readRDS("url_X_test.rds")
X<-as(X,"RsparseMatrix")
X_dim<-dim(X)
saveRDS(X_dim, file='url_X_test_dim.rds')
rm(X_dim)
X_p<-X@p
saveRDS(X_p, file='url_X_test_p.rds')
rm(X_p)
X_j<-X@j
saveRDS(X_j, file='url_X_test_j.rds')
rm(X_j)
X_x<-X@x
saveRDS(X_x, file='url_X_test_x.rds')
rm(list=ls())
y<-readRDS("url_y_training.rds")
y<-as.vector(y)
saveRDS(y, file='url_y_training_vec.rds')
rm(list=ls())
y<-readRDS("url_y_test.rds")
y<-as.vector(y)
saveRDS(y, file='url_y_test_vec.rds')
rm(list=ls())

