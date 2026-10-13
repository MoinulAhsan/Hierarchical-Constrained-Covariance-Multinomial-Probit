

run_multinomial_univariate_mcmc <- function(X, y,
                                            nit, nmc, nthin,
                                            tau2,lamda2,
                                            proposal_step = NULL,
                                            kappa = NULL,
                                            S_update_method = "joint",
                                            R_update_method = "grid") {
  
  
  
  # ---------------------------------------------------
  # Input:
  #   X              - design matrix (n x p)
  #   y              - response vector (n x 1, categorical 1:k)
  #   nit            - number of burn-in iterations
  #   nmc            - number of post-burn-in iterations to save
  #   nthin          - thinning interval
  #   tau2           - proposal variance for log(S) updates. refer to as tau2_phi in the manuscript
  #   lamda2          - prior variance
  #   proposal_step  - grid proposal window (required for R_update_method = "grid")
  #   kappa          - shifted-beta concentration (required for R_update_method = "shifted_beta")
  #   S_update_method - "joint" or "conditional"
  #   R_update_method - "grid" or "shifted_beta"
  #
  # Output:
  #   A list containing posterior draws and acceptance proportions
  # ---------------------------------------------------
  
  
  
  
  
  S_update_method <- match.arg(S_update_method, c("joint", "conditional"))
  R_update_method <- match.arg(R_update_method, c("grid", "shifted_beta"))
  
  # Require only the tuning parameter used by the selected R update.
  if (R_update_method == "grid") {
    if (is.null(proposal_step) || !is.numeric(proposal_step) ||
        length(proposal_step) != 1L || !is.finite(proposal_step) ||
        proposal_step <= 0) {
      stop('For R_update_method = "grid", provide a positive finite proposal_step.')
    }
  } else {
    if (is.null(kappa) || !is.numeric(kappa) ||
        length(kappa) != 1L || !is.finite(kappa) || kappa <= 2) {
      stop('For R_update_method = "shifted_beta", provide kappa > 2.')
    }
  }
  
  k<-length(unique(y)) # number of category
  
  x2<-X
  
  data<-data.frame(x2,y)
  
  
  
  
  x2<- as.matrix( data[,-(ncol(x2)+1)])
  y<- as.numeric(data[,ncol(data)])
  n<-length(y)
  p<-ncol(x2)
  
  
  
  XTX<-crossprod(x2)
  
  x_tilda<- diag((k-1)) %x% x2 # Construct block-diagonal matrix
  
  
  # Get positions (indices) of observations in y for each category (1 to k)
  category_position<-sapply(1:k,function(i)which(y==i))
  
  # Get the number of observations in y for each category
  category_frequency<-sapply(category_position,length)
  
  
  
  
  
  
  
  #------------------------------------------
  # Initialize prior mean and covariance for beta
  #------------------------------------------
  mu_beta<- matrix(rep(0,(k-1)*p ),ncol=1)
  sigma_beta<- 100 * diag( (k-1)*p )
  inv_sigma_beta<-solve_cal(sigma_beta)
  
  
  #------------------------------------------
  # Initialize parameters for MCMC / iterative algorithm
  #------------------------------------------
  
  # Initial beta vector (all zeros)
  beta<-matrix(rep(0,p*(k-1)),ncol=1)
  
  beta_tilda<-beta
  
  # Initial covariance matrix for categories (identity matrix)
  sigma<-diag( (k-1) )
  sigma_inv<-solve_cal(sigma) ## Inverse of sigma using Cholesky-based function
  
  
  
  
  
  
  # mean_w: n x (k-1) matrix of linear predictors
  mean_w<- matrix( x_tilda %*% beta_tilda , ncol=(k-1) )
  
  
  
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
  beta_final<- array(NA, c(p,(k-1), nsamples) )
  beta_final_vector<- array(NA, c((k-1)*p, nsamples) )
  Standardized_beta_final_vector<- array(NA, c((k-1)*p, nsamples) )
  
  
  sigma2_final<-array(NA, c(k-1, k-1, nsamples))
  
  
  S_final<-array(NA, c((k-1), nsamples) ) # Store standard deviations S for each iteration
  R_final<-array(NA, c( (k-1), (k-1), nsamples)) # Store correlation matrices R for each iteration
  
  
  
  
  # initial linear predictor
  X_BETA_tilda<-colSums(t(x_tilda)* c(beta))   
  
  X_BETA<-matrix( X_BETA_tilda, ncol= (k-1) ) 
  
  
  #defining the matrix positions indicating the correlations 
  correlation_positions<-expand.grid(1:((k-1)), (1):(k-1) )
  correlation_positions<-correlation_positions[,c(2,1)]
  #only choose the upper side of the diagonal
  correlation_positions<-correlation_positions[which(correlation_positions[,2]>correlation_positions[,1]), ]
  
  
  
  for( iter in 1 : total_iteration){
    
    
    
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
    
    A_beta<-sigma_inv %x% XTX + inv_sigma_beta # Compute posterior precision matrix for beta
    
    B_beta<- c(  matrix( colSums(t( sigma_inv %x% t(x2) ) * c(w_tilda)) ) + inv_sigma_beta %*%mu_beta  )     # Compute posterior mean vector for beta
    
    var1<- solve_cal(A_beta)          # Posterior covariance (inverse of precision)
    mu1<- colSums(var1 * B_beta)      # Posterior mean
    
    
    
    
    
    # Sample new beta from multivariate normal posterior
    beta_tilda<-  t(rmvnorm(1,mu1,  var1))
    
    beta<-beta_tilda
    
    
    
    #------------------------------------------
    # Compute updated linear predictor X * beta
    #------------------------------------------
    X_BETA_tilda<-colSums(t(x_tilda)* c(beta))   
    
    X_BETA<-matrix( X_BETA_tilda, ncol= (k-1) ) 
    
    
    
    # Residual latent variables (w - X*beta)
    W_Xbeta<-w_mat[,-c(1)] - X_BETA
    
    
    # Update S using the selected algorithm
    if (S_update_method == "conditional") {
      S_update <- update_conditional_S(
        W_Xbeta = W_Xbeta,
        S = S,
        log_S = log_S,
        sigma = sigma,
        R = R,
        tau2 = tau2,
        lamda2 = lamda2,
        accept_S = accept_S
      )
    } else {
      S_update <- update_S_jointly(
        log_S = log_S,
        S = S,
        sigma = sigma,
        R = R,
        W_Xbeta = W_Xbeta,
        tau2 = tau2,
        lamda2 = lamda2,
        rmvn_sumzero_constraint = rmvn_sumzero_constraint,
        log_dmvn_sumzero_constraint = log_dmvn_sumzero_constraint,
        accept_S = accept_S
      )
    }
    
    log_S <- S_update$log_S
    S <- S_update$S
    sigma <- S_update$sigma
    accept_S <- S_update$accept_S
    
    # Update R using the selected algorithm
    if (R_update_method == "grid") {
      R_update <- update_R_approach1(
        R = R,
        S = S,
        sigma = sigma,
        W_Xbeta = W_Xbeta,
        correlation_positions = correlation_positions,
        proposal_step = proposal_step,
        find_interval = find_interval,
        accept_R = accept_R
      )
    } else {
      R_update <- update_R_approach2(
        R = R,
        S = S,
        sigma = sigma,
        W_Xbeta = W_Xbeta,
        correlation_positions = correlation_positions,
        kappa = kappa,
        find_interval = find_interval,
        accept_R = accept_R
      )
    }
    
    R <- R_update$R
    sigma <- R_update$sigma
    accept_R <- R_update$accept_R
    
    sigma_inv<-solve_cal(sigma)  
    
    
    
    
    
    if(iter>nit &(iter-nit)%%nthin==0)
      
    {
      
      
      beta_final[ , , (iter-nit)/nthin] <- matrix( beta ,nrow = p )
      
      beta_final_vector[ ,  (iter-nit)/nthin]<-beta
      
      sigma2_final[, , (iter-nit)/nthin]<- sigma
      
      S_final[,(iter-nit)/nthin]<- S
      
      R_final[, , (iter-nit)/nthin]<- R
      
      
      
    }
    
    
  }
  
  
  
  
  # check MH acceptance rate
  accept_prop_S<-accept_S/total_iteration 
  accept_prop_R<-accept_R/(total_iteration*sum(upper.tri(R))) 
  
  
  result=list(beta_final_vector=beta_final_vector,sigma2_final=sigma2_final,S_final=S_final,R_final=R_final,accept_prop_S=accept_prop_S,accept_prop_R=accept_prop_R)
  
  
  
  return(result)
  
}




