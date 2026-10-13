


run_multinomial_bivariate_Conditional_S_mcmc<- function(X1,X2,y1,y2,
                                                         nit, nmc, nthin,
                                                         tau2,lamda2,
                                                         proposal_step) {
  
  
  
  # ---------------------------------------------------
  # Input:
  #   X              - design matrix (n x p)
  #   y              - response vector (n x 1, categorical 1:k)
  #   nit            - number of burn-in iterations
  #   nmc            - number of post-burn-in iterations to save
  #   nthin          - thinning interval
  #   tau2           - proposal variance for log(S) updates
  #   lamda2          - prior variance
  #   proposal_step  - proposal window size for correlation R
  #
  # Output:
  #   A list containing posterior draws and acceptance proportions
  # ---------------------------------------------------
  
  
  
  k1<-length(sort( unique(y1),decreasing = T )) # number of category for y1
  k2<-length(sort( unique(y2),decreasing = T )) # number of category for y2
  
  p1<-ncol(X1) 
  p2<-ncol(X2)
  
  x1<- X1
  x2<- X2
  
  
  data1<-data.frame(X1,y1)
  data2<-data.frame(X2,y2)
  
  
  n<-nrow(X1)
  
  
  
  
  y12<- cbind(y1,y2)
  
  category_position_y12 = as.list(rep(NA,k1*k2))
  dim(category_position_y12) = c(k1,k2)
  for (i in 1:k1)
  {
    for (j in 1:k2)
    {
      category_position_y12[[i,j]] = which(apply(y12, 1, function(a) { max(abs(a - c(i,j)))} )== 0) ## which observations have 1st y category = i AND 2nd y category = j
    }
  }
  
  
  
  
  
  
  big_x1<-diag((k1-1)) %x% x1 
  big_x2<-diag((k2-1)) %x% x2 
  
  x_tilda<-rbind(cbind(big_x1,matrix(0,ncol=ncol(big_x2),nrow=nrow(big_x1))),cbind(matrix(0,ncol=ncol(big_x1),nrow=nrow(big_x2)),big_x2))
  
  
  
  category_position_y1<-sapply(1:k1,function(i) which(y1==i)) 
  category_position_y2<-sapply(1:k2,function(i) which(y2==i)) 
  
  
  category_frequency_y1<-sapply(category_position_y1,length)
  category_frequency_y2<-sapply(category_position_y2,length)
  
  
  
  XTX_1<- crossprod(x1)
  XTX_2<- crossprod(x2)
  
  X1TX2<-crossprod(x1,x2)
  X2TX1<-t(X1TX2)
  
  X_tilda_T_X_tilda<-crossprod(x_tilda) 
  
  
  
  
  
  
  
  
  mu_beta<- matrix(rep(0,(((k1-1)*p1)+((k2-1)*p2)) ),ncol=1)
  sigma_beta<- 100 * diag((((k1-1)*p1)+((k2-1)*p2)))
  inv_sigma_beta<-solve_cal(sigma_beta)
  
  
  
  #initial   #######################
  
  beta<-matrix(rep(0,(((k1-1)*p1)+((k2-1)*p2))),ncol=1)
  
  beta_tilda<-beta
  
  
  
  
  sigma<-diag( (k1-1)+(k2-1) ) 
  sigma_inv<-solve_cal(sigma)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  mean_w<- matrix( x_tilda %*% beta_tilda , ncol=(k1-1)+(k2-1) )
  
  
  
  #------------------------------------------
  # Initialize w matrix for latent variables
  #------------------------------------------
  w_mat<-cbind(cbind(matrix(0,nrow=n,ncol=1),matrix(NA,n,(k1-1) ) ),cbind(matrix(0,nrow=n,ncol=1),matrix(NA,n,(k2-1) ) ))
  
  
  m1<-c(mean_w[category_position_y1[[1]],c(1:(k1-1))])
  v1<-sqrt( rep(diag(sigma[1:(k1-1),1:(k1-1)]),each=length(category_position_y1[[1]]) ) )
  
  
  w_mat[category_position_y1[[1]],-c(1,(k1+1):ncol(w_mat))]<-qnorm(pnorm(-Inf,m1 , v1 )+runif(category_frequency_y1[1]*(k1-1))*(pnorm(0, m1 , v1 )-pnorm(-Inf, m1,v1 )),m1,v1 )
  
  
  
  m2<-c(mean_w[category_position_y2[[1]],c(k1:(k1+k2-2))])
  v2<-sqrt( rep(diag(sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))]),each=length(category_position_y2[[1]]) ) )
  
  
  w_mat[category_position_y2[[1]],-c(1:(k1+1))]<-qnorm(pnorm(-Inf,m2 , v2 )+runif(category_frequency_y2[1]*(k2-1))*(pnorm(0, m2 , v2 )-pnorm(-Inf, m2,v2 )),m2,v2 )
  
  
  
  for(i in 2:k1){
    
    a<-0
    b<-Inf
    m<- c(mean_w[category_position_y1[[i]],(i-1)])
    v<-sqrt(diag(sigma[1:(k1-1),1:(k1-1)]) )[i-1]
    
    
    w_mat[category_position_y1[[i]],i]<-qnorm(pnorm(a,m, v)+runif(category_frequency_y1[i]) *(pnorm(b, m , v)-pnorm(a,m ,v )), m ,v)
    
    a<- -Inf
    b<-(w_mat[category_position_y1[[i]],i])
    
    
    m2<- c(mean_w[category_position_y1[[i]],-c((i-1),(k1:(k1+k2-2)))]) +  c(  t(  (sigma[1:(k1-1),1:(k1-1)])[-(i-1),(i-1)] %*% solve_cal((sigma[1:(k1-1),1:(k1-1)])[(i-1),(i-1)])   %*% t(matrix( w_mat[category_position_y1[[i]],(i)] - (mean_w[category_position_y1[[i]],(i-1)])  ))   )  )
    vv<- (sigma[1:(k1-1),1:(k1-1)])[-(i-1),-(i-1)] - (sigma[1:(k1-1),1:(k1-1)])[-(i-1),(i-1)] %*% solve_cal((sigma[1:(k1-1),1:(k1-1)])[(i-1),(i-1)]) %*% ((sigma[1:(k1-1),1:(k1-1)])[(i-1),-(i-1)])
    
    v2<- sqrt( matrix(diag(vv), nrow=category_frequency_y1[i],ncol=ncol(vv),byrow = T) )
    
    
    w_mat[category_position_y1[[i]],-c(1,i,(k1+1):(k1+k2) )]<-qnorm(pnorm(a,m2 , v2 )+runif(category_frequency_y1[i]*(k1-2))*(pnorm(b, m2 , v2 )-pnorm(a, m2 , v2 )),m2 ,v2 )
    
  }
  
  
  
  
  
  for(i in 2:k2){
    
    a<-0
    b<-Inf
    m<- c(mean_w[category_position_y2[[i]],(i+(k1-2))])
    v<-sqrt(diag(sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))]) )[i-1]
    
    
    w_mat[category_position_y2[[i]],(i+k1)]<-qnorm(pnorm(a,m, v)+runif(category_frequency_y2[i]) *(pnorm(b, m , v)-pnorm(a,m ,v )), m ,v)
    
    a<- -Inf
    b<-(w_mat[category_position_y2[[i]],(i+k1)])
    
    
    m2<- c(mean_w[category_position_y2[[i]],-c((i+(k1-2)),1:(k1-1))]) +  c(  t(  (sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[-(i-1),(i-1)] %*% solve_cal((sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),(i-1)])   %*% t(matrix( w_mat[category_position_y2[[i]],(i+k1)] - (mean_w[category_position_y2[[i]],(i+(k1-2))])  ))   )  )
    vv<- (sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[-(i-1),-(i-1)] - (sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[-(i-1),(i-1)] %*% solve_cal((sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),(i-1)]) %*% ((sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),-(i-1)])
    
    v2<- sqrt( matrix(diag(vv), nrow=category_frequency_y2[i],ncol=ncol(vv),byrow = T) )
    
    
    w_mat[category_position_y2[[i]],-c((i+k1),1:(k1+1))]<-qnorm(pnorm(a,m2 , v2 )+runif(category_frequency_y2[i]*(k2-2))*(pnorm(b, m2 , v2 )-pnorm(a, m2 , v2 )),m2 ,v2 )
    
  }
  
  w_tilda<-as.matrix( c( (w_mat[,-c(1,(k1+1))]) ) )
  
  
  
  
  
  
  
  k <- (k1-1)+(k2-1)  
  
  
  
  
  accept_S<-0
  accept_R<-0
  
  
  
  # Initialize S and R
  
  log_S <- rep(0, ((k1-1)+(k2-1)))     # Initial guess for log_S. log(S) is 0 tham means S is 1. so SD is 1
  S <- rep(1, ((k1-1)+(k2-1)))    # standard deviation
  R<- diag(1, ((k1-1)+(k2-1)))  # Start with identity matrix for R
  
  
  
  
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
  
  
  
  
  sigma2_final<-array(NA, c( (k1-1)+(k2-1), (k1-1)+(k2-1), nsamples))
  
  beta_final_vector<- array(NA, c( (((k1-1)*p1)+((k2-1)*p2)), nsamples) )
  
  Standardized_beta_final_vector<- array(NA, c( (((k1-1)*p1)+((k2-1)*p2)), nsamples) )
  
  
  S_final<-array(NA, c((k1-1)+(k2-1), nsamples) )
  R_final<-array(NA, c( (k1-1)+(k2-1), (k1-1)+(k2-1), nsamples))
  
  
  
  
  
  
  
  X_BETA_tilda<-colSums(t(x_tilda)* c(beta_tilda))   
  
  X_BETA<-matrix( X_BETA_tilda, ncol= (k1-1)+(k2-1) ) 
  
  
  
  #defining the matrix positions indicating the correlations 
  correlation_positions<-expand.grid(1:((k1+k2-2)), (1):(k1+k2-2) )
  correlation_positions<-correlation_positions[,c(2,1)]
  #only choose the upper side of the diagonal
  correlation_positions<-correlation_positions[which(correlation_positions[,2]>correlation_positions[,1]), ]
  
  
  
  for( iter in 1 : total_iteration)
  {
    
    
    
    ### ----- Update w_mat for first response (k1 categories) ----- ###
    
    for(i in 2:k1)
    {
      ## ---- Conditional distribution parameters for category i ----
      mean_A_given_B<- X_BETA[,(i-1)] + t( (sigma[1:(k1-1),1:(k1-1)])[(i-1),-(i-1)] %*% solve_cal((sigma[1:(k1-1),1:(k1-1)])[-(i-1),-(i-1)]) %*% t( w_mat[,-c(i,1,((k1+1):(k1+k2)))] - X_BETA[,-c((i-1),(k1:(k1+k2-2)))]) )
      sigma_A_given_B<- (sigma[1:(k1-1),1:(k1-1)])[(i-1),(i-1)] - (sigma[1:(k1-1),1:(k1-1)])[(i-1),-(i-1)] %*% solve_cal((sigma[1:(k1-1),1:(k1-1)])[-(i-1),-(i-1)]) %*% as.matrix((sigma[1:(k1-1),1:(k1-1)])[(i-1),-(i-1)])
      
      
      ## ---- Update latent variables for *observed* category i ----
      ## Lower truncation point: maximum of competing category latent variables
      trunc_lower_obs<-apply(w_mat[category_position_y1[[i]], -c(i,(k1+1):(k1+k2))], 1, max)
      trunc_upper_obs<-Inf ## No upper truncation (open interval)
      
      mean_obs<- mean_A_given_B[category_position_y1[[i]]]
      sd_obs<- sqrt(sigma_A_given_B)
      
      ## Sample from truncated normal: (a, ∞)
      w_mat[category_position_y1[[i]],i]<-qnorm(pnorm(trunc_lower_obs,mean_obs, sd_obs)+runif(category_frequency_y1[i])*(pnorm(trunc_upper_obs, mean_obs , sd_obs)-pnorm(trunc_lower_obs,mean_obs ,sd_obs)), mean_obs ,sd_obs)
      
      
      ## ---- Update latent variables for *non-category i* ----
      ## Upper truncation point: maximum of competing categories (must stay below)
      trunc_lower_non<- -Inf
      trunc_upper_non<-apply(w_mat[-category_position_y1[[i]], c((1:k1))], 1, max )
      
      mean_non<- mean_A_given_B[-category_position_y1[[i]]]
      sd_non<- sqrt(sigma_A_given_B)
      
      ## Sample from truncated normal: (-∞, b)
      w_mat[-category_position_y1[[i]],i]<-qnorm(pnorm(trunc_lower_non,mean_non , sd_non )+runif(n-category_frequency_y1[i])*(pnorm(trunc_upper_non, mean_non , sd_non )-pnorm(trunc_lower_non, mean_non , sd_non )),mean_non ,sd_non )
      
      
    }
    
    
    
    
    ### ----- Update w_mat for second response (k2 categories) ----- ###
    
    for(i in 2:k2)
    {
      
      ## ---- Conditional distribution parameters for category i ----
      mean_A_given_B<- X_BETA[,(i+(k1-2))] + t( (sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),-(i-1)] %*% solve_cal((sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[-(i-1),-(i-1)]) %*% t( w_mat[,-c((i+k1),(1:(k1+1)) )] - X_BETA[,-c((i+(k1-2)),(1:(k1-1)) )]) )
      sigma_A_given_B<- (sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),(i-1)] - (sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),-(i-1)] %*% solve_cal((sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[-(i-1),-(i-1)]) %*% as.matrix((sigma[(k1:(k1+k2-2)),(k1:(k1+k2-2))])[(i-1),-(i-1)])
      
      
      ## ---- Update latent variables for *observed* category i ----
      trunc_lower_obs<-apply(w_mat[category_position_y2[[i]], -c((i+k1),(1:k1))], 1, max)
      trunc_upper_obs<-Inf
      mean_obs<- mean_A_given_B[category_position_y2[[i]]]
      sd_obs<- sqrt(sigma_A_given_B)
      
      w_mat[category_position_y2[[i]],(i+k1)]<-qnorm(pnorm(trunc_lower_obs,mean_obs, sd_obs)+runif(category_frequency_y2[i])*(pnorm(trunc_upper_obs, mean_obs , sd_obs)-pnorm(trunc_lower_obs,mean_obs ,sd_obs)), mean_obs ,sd_obs)
      
      
      ## ---- Update latent variables for *non-category i* ----
      trunc_lower_non<- -Inf
      trunc_upper_non<-apply(w_mat[-category_position_y2[[i]], c((k1+1):(k1+k2))], 1, max )
      
      mean_non<- mean_A_given_B[-category_position_y2[[i]]]
      sd_non<- sqrt(sigma_A_given_B)
      
      w_mat[-category_position_y2[[i]],(i+k1)]<-qnorm(pnorm(trunc_lower_non,mean_non , sd_non )+runif(n-category_frequency_y2[i])*(pnorm(trunc_upper_non, mean_non , sd_non )-pnorm(trunc_lower_non, mean_non , sd_non )),mean_non ,sd_non )
      
      
    }
    
    
    
    
    w_tilda<- c(w_mat[,-c(1,(k1+1))])
    
    
    
    
    ### ----- Regression Coefficient Update (β) ----- ###
    
    
    
    ## -------------------------------------------------
    ## 1. Construct Posterior Precision Matrix for β
    ## -------------------------------------------------
    ## A_β = posterior precision = prior precision + data contribution
    A_beta<-rbind(cbind(sigma_inv[1:(k1-1),1:(k1-1)]%x%XTX_1,sigma_inv[1:(k1-1),(k1:(k1+k2-2))]%x%X1TX2),
                  cbind(sigma_inv[(k1:(k1+k2-2)),1:(k1-1)]%x%X2TX1,sigma_inv[(k1:(k1+k2-2)),(k1:(k1+k2-2))]%x%XTX_2))+ inv_sigma_beta
    
    ## -------------------------------------------------
    ## 2. Construct Posterior mean Vector (B_β)
    ## -------------------------------------------------
    ## This is the data contribution to the posterior mean of β.
    B_beta<-  colSums(  t(cbind(rbind((sigma_inv[1:(k1-1), 1:(k1-1)] %x% t(x1)),(sigma_inv[(k1:(k1+k2-2)), 1:(k1-1)] %x% t(x2))),
                                rbind((sigma_inv[1:(k1-1), (k1):(k1+k2-2)] %x% t(x1)),  (sigma_inv[(k1):(k1+k2-2), (k1):(k1+k2-2)] %x% t(x2))) ) ) *w_tilda   )
    
    
    
    ## -------------------------------------------------
    ## 3. Compute Posterior Mean and Covariance
    ## -------------------------------------------------
    var_beta <- solve_cal(A_beta)
    mu_beta <- colSums(var_beta  * c(B_beta))      
    
    
    
    ## -------------------------------------------------
    ## 4. Sample β from Multivariate Normal Posterior
    ## -------------------------------------------------
    beta_tilda<-  t(rmvnorm(1,mu_beta ,  var_beta ))
    
    beta<-beta_tilda
    
    
    
    ## -------------------------------------------------
    ## 5. Compute Xβ for Use in Latent Variable Update
    ## -------------------------------------------------
    X_BETA_tilda<-colSums( t(x_tilda)* c(beta_tilda) )        
    
    X_BETA<-matrix( X_BETA_tilda, ncol= (k1-1)+(k2-1) ) 
    
    
    
    
    
    ## Residual matrix: W - Xβ
    W_Xbeta<-w_mat[,-c(1,(k1+1))] - X_BETA
    
    
    
    
    
    ##############################
    # Update S (standard deviations) separately
    ##############################
    
    for(i in 1:(k1-1)){
      
      # Proposal: log-normal random walk for log(S[i])
      # Since S must be positive, we work on log(S)
      
      log_S_prop<-log_S
      
      delta_prop<-rnorm(1,0,sd=sqrt(tau2)) #proposed tau2 to adjust the acceptance proportion
      log_S_prop[i]<-log_S[i]+delta_prop
      log_S_prop[setdiff(1:(k1-1),i)]<-log_S[setdiff(1:(k1-1),i)] - (delta_prop/(k1-2))
      
      
      S_prop<-exp(log_S_prop)   # transform back to standard deviation scale
      
      
      # Update covariance matrix sigma with the proposed S[i]
      sigma_prop <-t(S_prop * R)*S_prop
      
      
      
      
      # Compute Cholesky transformation for likelihood calculation
      mat_prop<-solve(t(chol(sigma_prop)))
      mat<-solve(t(chol(sigma)))
      
      # Log-likelihood under proposed and current sigma
      log_likelihood_star<-n*sum(log(diag(mat_prop))) + sum(dnorm( mat_prop %*% t(W_Xbeta) ,log=T))    #      log(det(mat_prop) ) =sum(log(diag(mat_prop)))
      log_likelihood_current<-n*sum(log(diag(mat))) + sum(dnorm( mat %*% t(W_Xbeta) ,log=T))
      
      
      # Prior for log(S[i]): assume Normal(0, lamda2)
      log_prior_star <- sum(dnorm((log_S_prop[1:(k1-1)]),mean=0,sd=sqrt(lamda2),log=T))
      log_prior_current <- sum(dnorm((log_S[1:(k1-1)]),mean=0,sd=sqrt(lamda2),log=T))
      
      
      # Metropolis-Hastings acceptance ratio
      log_acceptance_ratio <- (log_likelihood_star + log_prior_star) - (log_likelihood_current + log_prior_current)
      acceptance_ratio <- exp(log_acceptance_ratio)
      
      # Accept or reject the proposed update
      if (runif(1) < min(1, acceptance_ratio)) {
        S<-S_prop
        
        log_S<-log_S_prop
        
        sigma <-sigma_prop
        
        accept_S<-accept_S+1 
        
      }
      
      
      
    }
    

    
    
    for(i in k1:ncol(W_Xbeta)){
      
      # Proposal: log-normal random walk for log(S[i])
      # Since S must be positive, we work on log(S)
      log_S_prop<-log_S
      
      delta_prop<-rnorm(1,0,sd=sqrt(tau2)) #proposed tau2 to adjust the acceptance proportion
      log_S_prop[i]<-log_S[i]+delta_prop
      log_S_prop[setdiff((k1:(k1+k2-2)),i)]<-log_S[setdiff((k1:(k1+k2-2)),i)] - (delta_prop/(k2-2))
      
      
      S_prop<-exp(log_S_prop)   # transform back to standard deviation scale
      
      
      # Update covariance matrix sigma with the proposed S[i]
      sigma_prop <-t(S_prop * R)*S_prop
      
      
      
      mat_prop<-solve(t(chol(sigma_prop)))
      mat<-solve(t(chol(sigma)))
      
      # Log-likelihood under proposed and current sigma
      log_likelihood_star<-n*sum(log(diag(mat_prop))) + sum(dnorm( mat_prop %*% t(W_Xbeta) ,log=T))    #      log(det(mat_prop) ) =sum(log(diag(mat_prop)))
      log_likelihood_current<-n*sum(log(diag(mat))) + sum(dnorm( mat %*% t(W_Xbeta) ,log=T))
      
      
      # Prior for log(S[i]): assume Normal(0, lamda2)
      log_prior_star <- sum(dnorm((log_S_prop[(k1:(k1+k2-2))]),mean=0,sd=sqrt(lamda2),log=T))
      log_prior_current <- sum(dnorm((log_S[(k1:(k1+k2-2))]),mean=0,sd=sqrt(lamda2),log=T))
      
      
      
      # Metropolis-Hastings acceptance ratio
      log_acceptance_ratio <- (log_likelihood_star + log_prior_star) - (log_likelihood_current + log_prior_current)
      acceptance_ratio <- exp(log_acceptance_ratio)
      
      # Accept or reject the proposed update
      if (runif(1) < min(1, acceptance_ratio)) {
        S<-S_prop
        
        log_S<-log_S_prop
        
        sigma <-sigma_prop
        
        accept_S<-accept_S+1  
        
      }
      
      
      
    }
    
    
    
    
    ################
    #simulate R
    ################
    
    
    for(ind in 1:nrow(correlation_positions)){
      
      i<-correlation_positions[ind,1]
      j<-correlation_positions[ind,2]
      
      
      #------------------------------------------
      # Step 1: Construct proposal grid for R[i,j]
      #------------------------------------------
      interval <- find_interval(R, i, j) # valid interval for correlation
      
      # Create grid of candidate values near current R[i,j]
      grid <- seq(max(R[i,j]-proposal_step,ceiling(interval[1]*100)/100,-0.99), min(R[i,j]+proposal_step,floor(interval[2]*100)/100,+0.99), by=.01)
      
      
      prop<-sample(grid, 1)   ## sample one from the grid
      
      R_prop<-R
      R_prop[i,j]<-prop
      R_prop[j,i]<-prop
      
      # Compute valid interval/grid for proposed value
      interval_prop <- find_interval(R_prop, i, j)
      grid_prop <- seq(max(R_prop[i,j]-proposal_step,ceiling(interval_prop[1]*100)/100, -0.99), min(R_prop[i,j]+proposal_step,floor(interval_prop[2]*100)/100,+0.99), by=.01)
      
      
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
      
      
      
      beta_final_vector[ ,  (iter-nit)/nthin]<-beta
      
      # loglikelihood_store_new[(iter-nit)/nthin]<-unlist(log_likelihood_mnp_2response(category_position_y12,k1,k2,mu1=X_BETA[,1:(k1-1)],mu2=X_BETA[,((k1):ncol(X_BETA))],sigma,r=250,Halton_choice=FALSE)[2])
      
      
      sigma2_final[, , (iter-nit)/nthin]<- sigma
      
      S_final[,(iter-nit)/nthin]<- S
      
      R_final[, , (iter-nit)/nthin]<- R
      
      
      
      
    }
    
    
    
  }
  
  
  
  
  # check MH acceptance rate
  accept_prop_S<-accept_S/(total_iteration*ncol(W_Xbeta)) 
  accept_prop_R<-accept_R/(total_iteration*sum(upper.tri(R))) 
  
  
  
  result=list(beta_final_vector=beta_final_vector,sigma2_final=sigma2_final,S_final=S_final,R_final=R_final,accept_prop_S=accept_prop_S,accept_prop_R=accept_prop_R)
  
  
  return(result)
  
}













