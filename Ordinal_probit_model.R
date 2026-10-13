

get_BPIC_ordinal_probit_model <- function(X, y) {

data<-data.frame(y,X)
train_data<- as.matrix( data)

x1<- as.matrix( train_data[,-1])
y1<- as.numeric(train_data[,1])


XTX<-crossprod(x1)


n<-nrow(x1)
p<-ncol(x1)


y1_num<-y1

a<-y1_num
b<-y1_num+1


category_position<-sapply(1:k,function(i) which(y1_num==i))   


# Prior specification for regression coefficients
mu_beta<-matrix(rep(0,p),ncol=1)
sigma_beta<- 10000 * diag(p)  
inv_sigma_beta<-solve_cal(sigma_beta)




#initial
beta_0<- - qnorm(  sum((y1==1))/n   )
beta<-matrix(c(beta_0,rep(0,(p-1)) ) )


# Empirical cumulative probabilities for threshold initialization
p1<-  sum((y1_num==1))/n   
p12<- sum((y1_num==1) | (y1_num==2))/n 
p123<- sum((y1_num==1) | (y1_num==2)| (y1_num==3))/n
p1234<- sum((y1_num==1) | (y1_num==2)| (y1_num==3)| (y1_num==4))/n


# Threshold (cutpoint) initialization
alpha<-array(0,k+1)
alpha[1]<- -Inf
alpha[2]<- 0
alpha[k+1]<- Inf

alpha[3]<- qnorm(p12)+beta_0
alpha[4]<- qnorm(p123)+beta_0



X_beta<- colSums(t(x1)*c(beta))




nsamples<-nmc/nthin
total_iteration<- nmc+nit

beta_final1<-matrix(NA,nsamples,p)
alpha_final<-matrix(NA,nsamples,k+1)



for( iter in 1:total_iteration)
{
  
  
  w<-qnorm(pnorm(alpha[a],X_beta , 1)+runif(length(y1))*(pnorm(alpha[b], X_beta , 1)-pnorm(alpha[a], X_beta,1)),X_beta,1)
  
  A<-inv_sigma_beta + XTX
  
  B<-inv_sigma_beta %*% mu_beta + colSums( (x1) * c(w) )
  
  var_beta<- solve_cal(A)
  mean_beta<- colSums(var_beta * c(B) )
  
  beta<-t(rmvnorm(1,mean_beta,var_beta))
  
  X_beta<- colSums(t(x1)*c(beta))
  
  
  
  for(i in 3:k)
  {
    
    m<-max(w[category_position[[i-1]]]) 
    n<-min(w[category_position[[i]]])
    
    alpha[i]<-runif(1,m,n)
    
  }
  
  
  
  
  if(iter>nit &(iter-nit)%%nthin==0)
    
  {
    beta_final1[(iter-nit)/nthin, ] <- beta
    alpha_final[(iter-nit)/nthin, ]<-alpha

  }
  
  
}




#########################
##DIC###
#########################

# Compute log-likelihood for each posterior draw
set.seed(12345)
log_likelihood_store <- numeric(nrow(beta_final1))
for (i in 1:nrow(beta_final1)) {
  beta_i  <- beta_final1[i, ]
  alpha_i <- alpha_final[i, ]
  X_beta_i <- colSums(t(x1) * c(beta_i))
  
  log_likelihood_store[i] <- sum(
    log(pnorm(alpha_i[y1_num + 1], X_beta_i, 1) -
          pnorm(alpha_i[y1_num],     X_beta_i, 1))
  )
}




alpha_bar <- apply(alpha_final, 2, mean)

beta_bar<-as.matrix( apply(beta_final1,2,mean) )

X_beta_bar<-x1 %*% beta_bar

set.seed(12345)
D_theta_bar<- -2*sum( log(pnorm(alpha_bar[y1_num+1],X_beta_bar,1)-pnorm(alpha_bar[y1_num],X_beta_bar,1)) )

D_bar<- -2*mean(log_likelihood_store)

p_D <- D_bar - D_theta_bar
BPIC <- D_bar + 2 * p_D

return(BPIC)
}

