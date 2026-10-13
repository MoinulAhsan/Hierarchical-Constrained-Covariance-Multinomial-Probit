
run_SSVS_Joint_S <- function(X, y, nit, nmc, nthin,
                    tau2,lamda2,proposal_step) {

data<-data.frame(y,X)
data<-na.omit(data)

x1<-data.frame(data[,-1])
y1 <- data$y
p<-ncol(x1)




train_data<-data.frame(y1,x1)

x1<- as.matrix( train_data[,-1])
y1<- as.numeric(train_data[,1])

n<-nrow(x1)
p<-ncol(x1)




y<-y1
x2<-x1


XTX<-crossprod(x2)

x_tilda<- diag((k-1)) %x% x2 


p<-ncol(x2)


# Get positions (indices) of observations in y for each category (1 to k)
category_position<-sapply(1:k,function(i)which(y==i))

# Get the number of observations in y for each category
category_frequency<-sapply(category_position,length)




mu_beta<- matrix(rep(0,(k-1)*p ),ncol=1)
sigma_beta<- 10 * diag( (k-1)*p )
inv_sigma_beta<-solve_cal(sigma_beta)

prior_beta_var<- diag(10,k-1)
inv_prior_beta_var<- solve(prior_beta_var)



#initial   #######################
beta<-matrix(rep(0,p*(k-1)),ncol=1)
beta_tilda<-beta



II<-rep(0,p) 



sigma<-diag( (k-1) )
sigma_inv<-solve_cal(sigma)






mean_w<- matrix( x_tilda %*% beta_tilda , ncol=k-1) 



#------------------------------------------
# Initialize w matrix for latent variables
#------------------------------------------

# w_mat: n x k matrix
# First column is 0 (baseline category)
w_mat<-cbind(matrix(0,nrow=n,ncol=1),matrix(NA,n,(k-1) ) )


mean_cat1<-c(mean_w[category_position[[1]],])
sd_cat1   <-sqrt( rep(diag(sigma),each=length(category_position[[1]]) ) )

# Sample truncated normal latent variables for category 1
w_mat[category_position[[1]],-1]<-qnorm(pnorm(-Inf,mean_cat1 , sd_cat1    )+runif(category_frequency[1]*(k-1))*(pnorm(0, mean_cat1 , sd_cat1)-pnorm(-Inf, mean_cat1,sd_cat1 )),mean_cat1,sd_cat1)


#------------------------------------------
# Fill in latent variables for remaining categories
#------------------------------------------
for(i in 2:k){
  
  # Truncation limits for the main diagonal
  lower_main <-0
  upper_main <-Inf
  
  # Linear predictor for observations in category i
  mean_main <- c(mean_w[category_position[[i]],(i-1)])
  sd_main <-sqrt(diag(sigma) )[i-1]
  
  # Sample truncated normal for the main diagonal of w_mat
  w_mat[category_position[[i]],i]<-qnorm(pnorm(lower_main,mean_main, sd_main)+runif(category_frequency[i]) *(pnorm(upper_main, mean_main , sd_main)-pnorm(lower_main,mean_main ,sd_main )), mean_main ,sd_main)
  
  
  # Truncation limits for off-diagonal latent variables
  lower_off <- -Inf
  upper_off <-(w_mat[category_position[[i]],i])
  
  # Conditional mean for off-diagonal elements
  cond_mean_off <- c(mean_w[category_position[[i]],-c(i-1)]) +  c(  t(  sigma[-(i-1),(i-1)] %*% solve_cal(sigma[(i-1),(i-1)])   %*% t(matrix( w_mat[category_position[[i]],(i)] - (mean_w[category_position[[i]],(i-1)])  ))   )  )
  # Conditional covariance for off-diagonal elements
  cond_cov_off <- sigma[-(i-1),-(i-1)] - sigma[-(i-1),(i-1)] %*% solve_cal(sigma[(i-1),(i-1)]) %*% (sigma[(i-1),-(i-1)])
  # Standard deviations for off-diagonal elements
  cond_sd_off <- sqrt( matrix(diag(cond_cov_off), nrow=category_frequency[i],ncol=ncol(cond_cov_off),byrow = T) )
  
  # Sample truncated normals for off-diagonal latent variables
  w_mat[category_position[[i]],c(-1,-i)]<-qnorm(pnorm(lower_off,cond_mean_off , cond_sd_off  )+runif(category_frequency[i]*(k-2))*(pnorm(upper_off, cond_mean_off , cond_sd_off  )-pnorm(lower_off, cond_mean_off , cond_sd_off  )),cond_mean_off ,cond_sd_off  )
  
}

w_tilda<-as.matrix( c( (w_mat[,c(2:k)]) ) )




#------------------------------------------
# Initialize parameters for latent variable covariance
#------------------------------------------

log_S <- rep(0, (k-1))     # Initial guess for log_S. log(S) is 0 means S is 1. so SD is 1
S <- rep(1, (k-1))    # standard deviation
R<- diag(1, (k-1))  # Start with identity matrix for R


# lamda2<-10  # Tuning parameter for proposal distribution 

accept_S<-0
accept_R<-0




compute_determinant <- function(R, i, j, r) {
  # Create a copy of R with r_ij replaced by r
  R_modified <- R
  R_modified[i, j] <- r
  R_modified[j, i] <- r
  
  # Calculate the determinant of the modified matrix
  det(R_modified)
}


find_interval <- function(R, i, j) {
  # Compute f(1), f(-1), and f(0) for the modified matrix
  f_1 <- compute_determinant(R, i, j, 1)
  f_neg1 <- compute_determinant(R, i, j, -1)
  f_0 <- compute_determinant(R, i, j, 0)
  
  # Calculate coefficients a, b, and c
  a <- (f_1 + f_neg1 - 2 * f_0) / 2
  b <- (f_1 - f_neg1) / 2
  c <- f_0
  
  # Calculate the discriminant
  discriminant <- b^2 - 4 * a * c
  
  # Find the roots if the discriminant is non-negative
  root1 <- (-b + sqrt(discriminant)) / (2 * a)
  root2 <- (-b - sqrt(discriminant)) / (2 * a)
  interval <- sort(c(root1, root2))
  return(interval)
  
}








nsamples<-nmc/nthin
total_iteration<- nmc+nit


#------------------------------------------
# Initialize storage arrays for MCMC / posterior samples
#------------------------------------------
beta_final_vector<- array(NA, c((k-1)*p, nsamples) )
Standardized_beta_final_vector<- array(NA, c((k-1)*p, nsamples) )


sigma2_final<-array(NA, c(k-1, k-1, nsamples))


S_final<-array(NA, c((k-1), nsamples) ) # Store standard deviations S for each iteration
R_final<-array(NA, c( (k-1), (k-1), nsamples)) # Store correlation matrices R for each iteration




X_BETA_tilda<-colSums(t(x_tilda)* c(beta))   

X_BETA<-matrix( X_BETA_tilda, ncol= (k-1) ) 




II_store<-matrix(NA,nrow=total_iteration,ncol=(k-1)*p)


#defining the matrix positions indicating the correlations 
correlation_positions<-expand.grid(1:((k-1)), (1):(k-1) )
correlation_positions<-correlation_positions[,c(2,1)]
#only choose the upper side of the diagonal
correlation_positions<-correlation_positions[which(correlation_positions[,2]>correlation_positions[,1]), ]




for( iter in 1 : total_iteration)
{
  
  
  ### w_mat update #####
  
  for(i in 2:k)
  {
    ## Conditional mean and variance
    mean_A_given_B<- X_BETA[,(i-1)] + t( sigma[(i-1),-(i-1)] %*% solve_cal(sigma[-(i-1),-(i-1)]) %*% t( w_mat[,-c(i,1)] - X_BETA[,-(i-1)]) )
    sigma_A_given_B<- sigma[(i-1),(i-1)] - sigma[(i-1),-(i-1)] %*% solve_cal(sigma[-(i-1),-(i-1)]) %*% as.matrix(sigma[(i-1),-(i-1)])
    
    
    ## ------------------------
    ## Update latent variable for *observed* category i
    ## ------------------------
    trunc_lower_obs <-apply(w_mat[category_position[[i]], -i], 1, max)
    trunc_upper_obs <-Inf
    
    mean_obs <- mean_A_given_B[category_position[[i]]]
    sd_obs   <- sqrt(sigma_A_given_B)
    
    w_mat[category_position[[i]],i]<-qnorm(pnorm(trunc_lower_obs,mean_obs, sd_obs)+runif(category_frequency[i])*(pnorm(trunc_upper_obs, mean_obs , sd_obs)-pnorm(trunc_lower_obs,mean_obs ,sd_obs)), mean_obs ,sd_obs)
    
    ## ------------------------
    ## Update latent variable for *non-category i* individuals
    ## ------------------------
    trunc_lower_non <- -Inf
    trunc_upper_non <-apply(w_mat[-category_position[[i]], ], 1, max )
    
    mean_non <- mean_A_given_B[-category_position[[i]]]
    sd_non   <- sqrt(sigma_A_given_B)
    
    w_mat[-category_position[[i]],i]<-qnorm(pnorm(trunc_lower_non,mean_non , sd_non )+runif(n-category_frequency[i])*(pnorm(trunc_upper_non, mean_non , sd_non )-pnorm(trunc_lower_non, mean_non , sd_non )),mean_non ,sd_non )
    
    
  }
  
  
  w_tilda<-as.matrix(c ( (w_mat[,c(2:k)]) ) ) 
  
  
  
  
  
  
  ###beta estimate ###
  
  
  ########################
  # SSVS combined all category for one covariate
  #######################
  for(bbb in 1:p ){
    aaa<-1:(k-1)
    gap<-(aaa-1)*p+bbb
    
    A<-sigma_inv * XTX[bbb,bbb] +inv_prior_beta_var 
    
    w_XB<-matrix(as.matrix(w_tilda)-(as.matrix(x_tilda[,-(gap)])%*%beta_tilda[-(gap)]))
    B<-sigma_inv %x% matrix(x2[,bbb],nrow = 1)    %*% w_XB
    
    
    kappa<-(1/sqrt(det(prior_beta_var %*% A)) ) * exp( .5*(  t(B)%*%solve_cal(A)%*%B ) )
    
    
    # II[bbb] = 1 => included in all categories
    # II[bbb] = 0 => excluded from all categories
    II[bbb]<-rbinom(n=1, size=1, prob= 1/(1+(1/kappa) )  )
    
    if(bbb==1)II[bbb]=1 # Force inclusion of first covariate (typically intercept)
    
    
    beta_tilda[gap]<-0
    
    # If covariate is included, sample coefficients jointly
    if(II[bbb]==1){
      
      var_beta<- solve_cal(A)
      mean_beta<- var_beta %*% B
      
      beta_tilda[gap]<-rmvnorm(1, mean_beta, (var_beta) )
    }
    
  }
  
  II_store[iter,]<-II
  
  
  
  
  
  beta<-beta_tilda
  
  
  
  
  #------------------------------------------
  # Compute updated linear predictor X * beta
  #------------------------------------------
  X_BETA_tilda<-colSums(t(x_tilda)* c(beta))   
  
  X_BETA<-matrix( X_BETA_tilda, ncol= (k-1) ) 
  
  
  
  # Residual latent variables (w - X*beta)
  W_Xbeta<-w_mat[,-c(1)] - X_BETA
  
  
  
  
  ##############################
  # Update S (standard deviations) simultaneously
  ##############################
  # Note:
  # - S follows a lognormal distribution
  # - Therefore, log(S) follows a normal distribution
  # - We use a constrained proposal to ensure sum(log(S)) = 0
  
  # Constraint matrix: 1^T (enforces sum-to-zero condition on log_S)
  ones_matrix <- matrix(1, ncol = (k-1)) # A=1^T
  
  # Proposal variance scaling parameter (tuning for acceptance rate)
  # tau2<-0.000625  # (adjustable; higher values → more variability)
  
  # Proposal distribution parameters
  mu_c<-log_S
  sigma_c<-tau2 *(diag(k-1)- (1/(k-1))*matrix(1,ncol=(k-1),nrow=(k-1)) )
  
  
  # Draw proposed log(S) satisfying constraint
  log_S_prop<-c(rmvn_sumzero_constraint(mu_vec =mu_c ,sigma_mat = sigma_c ))
  
  # Transform back to S (since S = exp(log_S))
  S_prop<-exp(log_S_prop)
  # Proposed covariance matrix for latent variables
  sigma_prop<-t(S_prop * R)*S_prop
  
  # Cholesky precision matrices
  mat_prop<-solve(t(chol(sigma_prop)))
  mat<-solve(t(chol(sigma)))
  
  # Proposed log-likelihood
  log_likelihood_star<-n*sum(log(diag(mat_prop))) + sum(dnorm( mat_prop %*% t(W_Xbeta) ,log=T))    #      log(det(mat_prop) ) =sum(log(diag(mat_prop)))
  # Current log-likelihood
  log_likelihood_current<-n*sum(log(diag(mat))) + sum(dnorm( mat %*% t(W_Xbeta) ,log=T))
  
  
  #------------------------------------------
  # Compute log-priors under sum-to-zero constraint
  #------------------------------------------
  
  mu_c_prior<-rep(0,(k-1))
  sigma_c_prior<-lamda2 *(diag(k-1)- (1/(k-1))*matrix(1,ncol=(k-1),nrow=(k-1)) )
  
  
  log_prior_star <-log_dmvn_sumzero_constraint(log(S_prop),mu_vec =mu_c_prior ,sigma_mat = sigma_c_prior  ) 
  log_prior_current <- log_dmvn_sumzero_constraint(log(S),mu_vec =mu_c_prior,sigma_mat = sigma_c_prior  )
  
  
  #------------------------------------------
  # Metropolis–Hastings acceptance step
  #------------------------------------------
  
  # Compute the acceptance ratio (posterior ratio)
  log_acceptance_ratio <- (log_likelihood_star + log_prior_star) - (log_likelihood_current + log_prior_current)
  acceptance_ratio <- exp(log_acceptance_ratio)
  
  
  # Accept or reject proposed S
  if (runif(1) < min(1, acceptance_ratio)) {
    S<-S_prop
    
    log_S<-log(S)
    
    sigma <-sigma_prop
    
    accept_S<-accept_S+1  # track acceptance count
    
  }
  
  
  
  
  ##############################
  # Update R (correlation matrix)
  ##############################
  
  for(ind in 1:nrow(correlation_positions)){
    
    i<-correlation_positions[ind,1]
    j<-correlation_positions[ind,2]
      
      
      #------------------------------------------
      # Step 1: Construct proposal grid for R[i,j]
      #------------------------------------------
      interval <- find_interval(R, i, j)  # valid interval for correlation
      
      # Create grid of candidate values near current R[i,j] 
      grid <- seq(max(R[i,j]-proposal_step,ceiling(interval[1]*100)/100+.01,-0.99), min(R[i,j]+proposal_step,floor(interval[2]*100)/100-.01,+0.99), by=.01)
      
      
      prop<-sample(grid, 1)   ## sample one from the grid
      
      R_prop<-R
      R_prop[i,j]<-prop
      R_prop[j,i]<-prop
      
      
      # Compute valid interval/grid for proposed value
      interval_prop <- find_interval(R_prop, i, j)
      grid_prop <- seq(max(R_prop[i,j]-proposal_step,ceiling(interval_prop[1]*100)/100+.01, -0.99), min(R_prop[i,j]+proposal_step,floor(interval_prop[2]*100)/100-.01,+0.99), by=.01)
      
      
      #------------------------------------------
      # Step 2: Compute likelihood for proposed vs current
      #------------------------------------------
      
      sigma_prop<-sigma
      sigma_prop[i,j]<-S[i]*R_prop[i,j]*S[j]
      sigma_prop[j,i]<-sigma_prop[i,j]
      
      
      mat_prop<-solve(t(chol(sigma_prop)))
      mat<-solve(t(chol(sigma)))
      
      log_likelihood_star<-n*sum(log(diag(mat_prop))) + sum(dnorm( mat_prop %*% t(W_Xbeta) ,log=T))    
      log_likelihood_current<-n*sum(log(diag(mat))) + sum(dnorm( mat %*% t(W_Xbeta) ,log=T))
      
      
      #------------------------------------------
      # Step 3: Proposal correction (due to grid sampling)
      #------------------------------------------
      
      log_prop_ratio<- -log(length(grid_prop))+log(length(grid))  
      
      #------------------------------------------
      # Step 4: Metropolis–Hastings acceptance step
      #------------------------------------------
      
      #Compute the acceptance ratio (posterior ratio)
      log_acceptance_ratio <- (log_likelihood_star ) - (log_likelihood_current)+log_prop_ratio
      
      acceptance_ratio <- exp(log_acceptance_ratio)
      
      accept_RR<-0
      
      if (runif(1) < min(1, acceptance_ratio)) {
        R<-R_prop
        
        sigma<-sigma_prop
        
        accept_RR<-1
        
      }
      
      accept_R<-accept_RR +accept_R
    }
    
  
  
  
  

  sigma_inv<-solve_cal(sigma)
  
  
  
  
  
  
  if(iter>nit &(iter-nit)%%nthin==0)
    
  {

    beta_final_vector[ ,  (iter-nit)/nthin]<-beta_tilda
    
    sigma2_final[, , (iter-nit)/nthin]<- sigma
    
    S_final[,(iter-nit)/nthin]<- S
    
    R_final[, , (iter-nit)/nthin]<- R
    
    
  }
  
  
}



# check MH acceptance rate
accept_prop_S<-accept_S/total_iteration
accept_prop_R<-accept_R/(total_iteration*sum(upper.tri(R))) 



# for combined category, rows of no_zero_beta will be all same 
# because if that covariate is significant then all the 3 category will be significant and not 0
no_zero_beta<-c()

for(i in (1:(p*(k-1))))
{
  no_zero_beta[i]<-(length(which(beta_final_vector[i,]==0)))
}




try_beta<-t(beta_final_vector)[,1:p]

try_beta[try_beta!=0]<-1

U_beta<-unique(try_beta)


AA<-dplyr::count(data.frame(try_beta), pick(everything()))

AA[order(AA[,(p+1)],decreasing = T),]


AA_sub <- colSums(AA[AA$n > nsamples * 0.01, ])
selected_X <- names(AA_sub)[AA_sub > 0 & names(AA_sub) != "n"]

idx <- as.numeric(sub("X", "", selected_X))
selected_names <- colnames(x2)[idx]

return(selected_names)
}




