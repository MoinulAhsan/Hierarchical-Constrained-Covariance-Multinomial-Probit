###############
# Probability Prediction

new_probability_mnp_1response = function(category,mu,sigma,r=250,Halton_choice=FALSE)
  
{
  k = ncol(sigma) + 1 ## how many categories 
  
  ## now run the loop for category 1 to k observations
  
  above=rep(0,k-1)
  
  j=category
  
  
  muj=mu
  Aj= -diag(1,k-1) ## construct I_(k-1) matrix and make it negative
  if (j > 1) Aj[,j-1]=rep(1,k-1) ## make the (j-1)-th column all entries equal to 1 for category 2 and higher
  trunpt=as.vector(-Aj%*%t(muj)) ## set the lower bounds
  Lj=t(chol(Aj%*%sigma%*%t(Aj))) ## find the lower triangular cholesky factor
  logl = ghkvec(Lj,trunpt,above,r,HALTON=Halton_choice,0)+1.0e-50
  #  note:  ghkvec does an entire vector of n probs each with different truncation points but the
  #            same cov matrix.  
  
  
  return(logl)
}

################################ 


# ##################

## multinomial log likelihood single response
## assume k categories, p covariates (including intercept)
## input: y is a n by 1 response vector
## input:  mu = X*beta is n by (k-1) matrix as we excluded the first category 
## input: sigma is a (k-1) by (k-1) matrix
## Before running MCMC, create a list for positions of observations belonging to different y categories as:
## category_position = sapply(1:k,function(a){which(y==a)})
## It will save significant computations in log-likelihood computations at every iteration within MCMC since
## we do not have to search which y values are from a certain category

log_likelihood_mnp_1response = function(category_position,mu,sigma,r=250,Halton_choice=FALSE)
  
{
  k = ncol(sigma) + 1 ## how many categories 
  
  logl = array(0,k)
  
  ## now run the loop for category 1 to k observations
  
  above=rep(0,k-1)
  
  for (j in 1:k)
  {
    muj=mu[category_position[[j]],]
    Aj= -diag(1,k-1) ## construct I_(k-1) matrix and make it negative
    if (j > 1) Aj[,j-1]=rep(1,k-1) ## make the (j-1)-th column all entries equal to 1 for category 2 and higher
    trunpt=as.vector(-Aj%*%t(muj)) ## set the lower bounds
    Lj=t(chol(Aj%*%sigma%*%t(Aj))) ## find the lower triangular cholesky factor
    logl[j] = sum(log(ghkvec(Lj,trunpt,above,r,HALTON=Halton_choice,0)+1.0e-50))
    #  note:  ghkvec does an entire vector of n probs each with different truncation points but the
    #            same cov matrix.  
    
  }
  return(list(category_specific_log_likelihood = logl,total_log_likelihood = sum(logl)))
}

################################ 



log_likelihood_mnp_1response_WAIC <- function(category_position, mu, sigma,
                                              r = 250, Halton_choice = FALSE) {
  k <- ncol(sigma) + 1
  logl <- vector("list", k)
  above <- rep(0, k - 1)
  
  for (j in 1:k) {
    muj <- mu[category_position[[j]], ]
    Aj <- -diag(1, k - 1)
    if (j > 1) Aj[, j - 1] <- 1
    
    trunpt <- as.vector(-Aj %*% t(muj))
    Lj <- t(chol(Aj %*% sigma %*% t(Aj)))
    
    logl[[j]] <- log(ghkvec(Lj, trunpt, above, r,
                            HALTON = Halton_choice, 0) + 1e-50)
  }
  
  list(
    observation_log_likelihood = unlist(logl),
    total_log_likelihood = sum(unlist(logl))
  )
}






################################################

## This code is necessary when
## working with MVN distributions
## derived using sum = 0 constraint.
## Under such linear constraints, dispersion matrix
## of MVN becomes singular. So, the usual commands
## such as dmvnorm or rmvnorm no longer works.

##function to generate a vector from MVN with sum = 0 constraints
## need to use eigenvalues and eigenvectors of dispersion matrix

rmvn_sumzero_constraint = function(mu_vec, sigma_mat) 
{
  eig = eigen(sigma_mat,only.values = F)
  eig_val = eig$values ## eigenvalues
  eig_vec = eig$vectors ## eigenvectors
  
  out = matrix(mu_vec,ncol=1) + eig_vec%*%matrix(rnorm(nrow(sigma_mat),0,sqrt(eig_val)),ncol=1)
  return(out)
}


## function to calculate the log(density) of MVN random_vec with sum = 0 constraint
## since dispersion matrix is singular, we cannot take inverse, we have to use g-inverse
## only the part in the exponential is calculated
## the part with normalizing constant not calculated
## since it will cancel in the MH ratio

log_dmvn_sumzero_constraint = function(random_vec, mu_vec, sigma_mat){-0.5*sum(c(random_vec-mu_vec)*(ginv(sigma_mat)%*%matrix(random_vec-mu_vec,ncol=1)))}

################################################


## Computes the inverse of a positive-definite matrix x using Cholesky decomposition
solve_cal<-function(x) {chol2inv(chol(x))} 





















#################################
# Bivariate
# ###############################



## multinomial log likelihood two responses
## assume k1 and k2 categories, p1 and p2 covariates (including intercept)
## input: y is a n by 2 response matrix
## input:  mu1 = X1*beta1 is n by (k1-1) matrix as we excluded the first category
## input:  mu2 = X2*beta2 is n by (k2-1) matrix as we excluded the first category
## input: sigma is a (k1+k2-2) by (k1+k2-2) matrix
## Before running MCMC, create a 2D list for positions of observations belonging to different y categories for both responses as:

# library(bayesm)

log_likelihood_mnp_2response = function(pos_y,k1,k2,mu1,mu2,sigma,r=250,Halton_choice=FALSE)
  
{
  mu_combined = cbind(mu1,mu2)
  
  logl = matrix(0,k1,k2)
  
  ## Now, run the loop for k1*k2 category combinations
  
  above=rep(0,k1+k2-2)
  
  for (i in 1:k1)
  {
    for (j in 1:k2)
    {
      if (length(pos_y[[i,j]])==0) next  ## skip any category combination that never occurred in the data
      muij=mu_combined[pos_y[[i,j]],]
      Aij= -diag(1,k1+k2-2) ## construct I_(k1-1+k2-1) matrix and make it negative
      if (i > 1) Aij[1:(k1-1),i-1]=rep(1,k1-1) ## make the left block (i-1)-th column top part equal to 1 if first response category is 2 or higher
      if (j > 1) Aij[k1:(k1+k2-2),(k1-1+j-1)]=rep(1,k2-1) ## make the right block (j-1)-th column bottom part equal to 1 if second response category is 2 or higher
      trunpt=as.vector(-Aij%*%t(muij)) ## set the lower bounds
      Lij=t(chol(Aij%*%sigma%*%t(Aij))) ## find the lower triangular cholesky factor
      logl[i,j] = sum(log(ghkvec(Lij,trunpt,above,r,HALTON=Halton_choice,0)+1.0e-50))
      #  note:  ghkvec does an entire vector of n probs each with different truncation points but the
      #            same cov matrix.
    }
  }
  #print(logl)
  return(list(category_combination_specific_log_likelihood = logl,total_log_likelihood = sum(logl)))
}

################################ 






##############################
log_likelihood_mnp_2response_WAIC <- function(
    pos_y,
    k1,
    k2,
    mu1,
    mu2,
    sigma,
    r = 250,
    Halton_choice = FALSE
) {
  
  mu_combined <- cbind(mu1, mu2)
  
  ## total number of observations
  n <- nrow(mu_combined)
  
  ## observation-specific log-likelihood
  obs_logl <- rep(NA_real_, n)
  
  above <- rep(0, k1 + k2 - 2)
  
  for (i in 1:k1) {
    
    for (j in 1:k2) {
      
      if (length(pos_y[[i, j]]) == 0)
        next
      
      muij <- mu_combined[pos_y[[i, j]], , drop = FALSE]
      
      Aij <- -diag(1, k1 + k2 - 2)
      
      ## First response
      if (i > 1)
        Aij[1:(k1 - 1), i - 1] <- 1
      
      ## Second response
      if (j > 1)
        Aij[k1:(k1 + k2 - 2), (k1 - 1 + j - 1)] <- 1
      
      trunpt <- as.vector(-Aij %*% t(muij))
      
      Lij <- t(chol(Aij %*% sigma %*% t(Aij)))
      
      ## one log-likelihood per observation
      obs_logl[pos_y[[i, j]]] <-
        log(
          ghkvec(
            Lij,
            trunpt,
            above,
            r,
            HALTON = Halton_choice,
            0
          ) + 1e-50
        )
    }
  }
  
  return(
    list(
      observation_log_likelihood = obs_logl,
      total_log_likelihood = sum(obs_logl)
    )
  )
}



###############################





###############
# Probability Predection


new_probability_mnp_2response = function(category,k1,k2,mu1,mu2,sigma,r=250,Halton_choice=FALSE)
  
{
  mu_combined = cbind(mu1,mu2)
  
  logl = matrix(0,k1,k2)
  
  ## Now, run the loop for k1*k2 category combinations
  
  above=rep(0,k1+k2-2)
  i=ceiling(category/k2)
  j=category-((i-1)*k2)
  
  # if (length(pos_y[[i,j]])==0) next  ## skip any category combination that never occurred in the data
  muij=mu_combined
  Aij= -diag(1,k1+k2-2) ## construct I_(k1-1+k2-1) matrix and make it negative
  if (i > 1) Aij[1:(k1-1),i-1]=rep(1,k1-1) ## make the left block (i-1)-th column top part equal to 1 if first response category is 2 or higher
  if (j > 1) Aij[k1:(k1+k2-2),(k1-1+j-1)]=rep(1,k2-1) ## make the right block (j-1)-th column bottom part equal to 1 if second response category is 2 or higher
  trunpt=as.vector(-Aij%*%t(muij)) ## set the lower bounds
  Lij=t(chol(Aij%*%sigma%*%t(Aij))) ## find the lower triangular cholesky factor
  logl = ghkvec(Lij,trunpt,above,r,HALTON=Halton_choice,0)+1.0e-50
  #  note:  ghkvec does an entire vector of n probs each with different truncation points but the
  #            same cov matrix.
  
  
  return(logl)
}

################################ 


















# MCMC helper functions
#########################################

update_conditional_S<- function(
    W_Xbeta,
    S,
    log_S,
    sigma,
    R,
    tau2,
    lamda2,
    accept_S
) {
  
  n <- nrow(W_Xbeta)
  k_minus_1 <- length(S)
  
  for (i in 1:k_minus_1) {
    
    #------------------------------------------
    # Proposal
    #------------------------------------------
    log_S_prop <- log_S
    
    delta_prop <- rnorm(
      1,
      mean = 0,
      sd = sqrt(tau2)
    )
    
    log_S_prop[i] <- log_S[i] + delta_prop
    log_S_prop[-i] <- log_S[-i] -
      delta_prop / (k_minus_1 - 1)
    
    S_prop <- exp(log_S_prop)
    
    #------------------------------------------
    # Proposed covariance matrix
    #------------------------------------------
    sigma_prop <- t(S_prop * R) * S_prop
    
    #------------------------------------------
    # Likelihood
    #------------------------------------------
    mat_prop <- solve(t(chol(sigma_prop)))
    mat      <- solve(t(chol(sigma)))
    
    log_likelihood_star <-
      n * sum(log(diag(mat_prop))) +
      sum(dnorm(mat_prop %*% t(W_Xbeta), log = TRUE))
    
    log_likelihood_current <-
      n * sum(log(diag(mat))) +
      sum(dnorm(mat %*% t(W_Xbeta), log = TRUE))
    
    #------------------------------------------
    # Prior
    #------------------------------------------
    log_prior_star <-
      sum(
        dnorm(
          log_S_prop,
          mean = 0,
          sd = sqrt(lamda2),
          log = TRUE
        )
      )
    
    log_prior_current <-
      sum(
        dnorm(
          log_S,
          mean = 0,
          sd = sqrt(lamda2),
          log = TRUE
        )
      )
    
    #------------------------------------------
    # MH ratio
    #------------------------------------------
    log_acceptance_ratio <-
      (log_likelihood_star + log_prior_star) -
      (log_likelihood_current + log_prior_current)
    
    if (log(runif(1)) < min(0, log_acceptance_ratio)) {
      
      S <- S_prop
      log_S <- log_S_prop
      sigma <- sigma_prop
      
      accept_S <- accept_S + 1
    }
  }
  
  return(list(
    S = S,
    log_S = log_S,
    sigma = sigma,
    accept_S = accept_S
  ))
}




update_S_jointly <- function(
    log_S,
    S,
    sigma,
    R,
    W_Xbeta,
    tau2,
    lamda2,
    rmvn_sumzero_constraint,
    log_dmvn_sumzero_constraint,
    accept_S
) {
  
  k_minus_1 <- length(log_S)
  n <- nrow(W_Xbeta)
  
  ##############################
  # Proposal distribution
  ##############################
  
  mu_c <- log_S
  
  sigma_c <- tau2 * (
    diag(k_minus_1) -
      (1 / k_minus_1) *
      matrix(1,
             nrow = k_minus_1,
             ncol = k_minus_1)
  )
  
  log_S_prop <- c(
    rmvn_sumzero_constraint(
      mu_vec = mu_c,
      sigma_mat = sigma_c
    )
  )
  
  S_prop <- exp(log_S_prop)
  
  sigma_prop <- t(S_prop * R) * S_prop
  
  ##############################
  # Log-likelihood
  ##############################
  
  mat_prop <- solve(t(chol(sigma_prop)))
  mat      <- solve(t(chol(sigma)))
  
  log_likelihood_star <-
    n * sum(log(diag(mat_prop))) +
    sum(dnorm(mat_prop %*% t(W_Xbeta), log = TRUE))
  
  log_likelihood_current <-
    n * sum(log(diag(mat))) +
    sum(dnorm(mat %*% t(W_Xbeta), log = TRUE))
  
  ##############################
  # Prior
  ##############################
  
  mu_c_prior <- rep(0, k_minus_1)
  
  sigma_c_prior <- lamda2 * (
    diag(k_minus_1) -
      (1 / k_minus_1) *
      matrix(1,
             nrow = k_minus_1,
             ncol = k_minus_1)
  )
  
  log_prior_star <-
    log_dmvn_sumzero_constraint(
      log_S_prop,
      mu_vec = mu_c_prior,
      sigma_mat = sigma_c_prior
    )
  
  log_prior_current <-
    log_dmvn_sumzero_constraint(
      log_S,
      mu_vec = mu_c_prior,
      sigma_mat = sigma_c_prior
    )
  
  ##############################
  # MH ratio
  ##############################
  
  log_acceptance_ratio <-
    (log_likelihood_star + log_prior_star) -
    (log_likelihood_current + log_prior_current)
  
  if (log(runif(1)) < min(0, log_acceptance_ratio)) {
    
    log_S <- log_S_prop
    S<-exp(log_S)
    sigma <- sigma_prop
    accept_S <- accept_S + 1
    
  }
  
  return(list(
    log_S = log_S,
    S = S,
    sigma = sigma,
    accept_S = accept_S
  ))
}





update_R_approach1 <- function(
    R,
    S,
    sigma,
    W_Xbeta,
    correlation_positions,
    proposal_step,
    find_interval,
    accept_R
) {
  
  n <- nrow(W_Xbeta)
  
  for (ind in 1:nrow(correlation_positions)) {
    
    i <- correlation_positions[ind, 1]
    j <- correlation_positions[ind, 2]
    
    #------------------------------------------
    # Proposal distribution
    #------------------------------------------
    interval <- find_interval(R, i, j)
    
    grid <- seq(
      max(R[i, j] - proposal_step,
          ceiling(interval[1] * 100) / 100,
          -0.99),
      min(R[i, j] + proposal_step,
          floor(interval[2] * 100) / 100,
          0.99),
      by = 0.01
    )
    
    prop <- sample(grid, 1)
    
    R_prop <- R
    R_prop[i, j] <- prop
    R_prop[j, i] <- prop
    
    interval_prop <- find_interval(R_prop, i, j)
    
    grid_prop <- seq(
      max(R_prop[i, j] - proposal_step,
          ceiling(interval_prop[1] * 100) / 100,
          -0.99),
      min(R_prop[i, j] + proposal_step,
          floor(interval_prop[2] * 100) / 100,
          0.99),
      by = 0.01
    )
    
    #------------------------------------------
    # Proposed covariance matrix
    #------------------------------------------
    sigma_prop <- sigma
    sigma_prop[i, j] <- S[i] * R_prop[i, j] * S[j]
    sigma_prop[j, i] <- sigma_prop[i, j]
    
    #------------------------------------------
    # Log-likelihood
    #------------------------------------------
    mat_prop <- solve(t(chol(sigma_prop)))
    mat      <- solve(t(chol(sigma)))
    
    log_likelihood_star <-
      n * sum(log(diag(mat_prop))) +
      sum(dnorm(mat_prop %*% t(W_Xbeta), log = TRUE))
    
    log_likelihood_current <-
      n * sum(log(diag(mat))) +
      sum(dnorm(mat %*% t(W_Xbeta), log = TRUE))
    
    #------------------------------------------
    # Proposal correction
    #------------------------------------------
    log_prop_ratio <- log(length(grid)) - log(length(grid_prop))
    
    #------------------------------------------
    # MH step
    #------------------------------------------
    log_acceptance_ratio <-
      log_likelihood_star -
      log_likelihood_current +
      log_prop_ratio
    
    if (log(runif(1)) < min(0, log_acceptance_ratio)) {
      
      R <- R_prop
      sigma <- sigma_prop
      accept_R <- accept_R + 1
      
    }
  }
  
  return(list(
    R = R,
    sigma = sigma,
    accept_R = accept_R
  ))
}








get_beta_parameters <- function(kappa, interval, r_current) {
  
  L <- interval[1]
  U <- interval[2]
  
  alpha <- ((kappa - 1) * L +(2 - kappa) * r_current -U) /(L - U)
  
  beta <- kappa - alpha
  
  return(c(alpha, beta))
}


update_R_approach2 <- function(
    R,
    S,
    sigma,
    W_Xbeta,
    correlation_positions,
    kappa,
    find_interval,
    accept_R
) {
  
  n <- nrow(W_Xbeta)
  
  
  for (ind in 1:nrow(correlation_positions)) {
    
    i <- correlation_positions[ind,1]
    j <- correlation_positions[ind,2]
    
    interval <- find_interval(R, i, j)
    
    beta_dist_par <- get_beta_parameters(
      kappa,
      interval,
      R[i,j]
    )
    
    prop <- interval[1] +
      (interval[2] - interval[1]) *
      rbeta(
        1,
        shape1 = beta_dist_par[1],
        shape2 = beta_dist_par[2]
      )
    
    R_prop <- R
    R_prop[i,j] <- prop
    R_prop[j,i] <- prop
    
    interval_prop <- find_interval(
      R_prop,
      i,
      j
    )
    
    beta_dist_par_prop <- get_beta_parameters(
      kappa,
      interval_prop,
      R_prop[i,j]
    )
    
    sigma_prop <- t(S * R_prop) * S
    
    mat_prop <- solve(t(chol(sigma_prop)))
    mat      <- solve(t(chol(sigma)))
    
    log_likelihood_star <-
      n * sum(log(diag(mat_prop))) +
      sum(dnorm(mat_prop %*% t(W_Xbeta), log = TRUE))
    
    log_likelihood_current <-
      n * sum(log(diag(mat))) +
      sum(dnorm(mat %*% t(W_Xbeta), log = TRUE))
    
    log_q_current_to_prop <-
      dbeta(
        (prop - interval[1]) /
          (interval[2] - interval[1]),
        shape1 = beta_dist_par[1],
        shape2 = beta_dist_par[2],
        log = TRUE
      ) -
      log(interval[2] - interval[1])
    
    log_q_prop_to_current <-
      dbeta(
        (R[i,j] - interval_prop[1]) /
          (interval_prop[2] - interval_prop[1]),
        shape1 = beta_dist_par_prop[1],
        shape2 = beta_dist_par_prop[2],
        log = TRUE
      ) -
      log(interval_prop[2] - interval_prop[1])
    
    log_prop_ratio <-
      log_q_prop_to_current -
      log_q_current_to_prop
    
    log_acceptance_ratio <-
      log_likelihood_star -
      log_likelihood_current +
      log_prop_ratio
    
    if (log(runif(1)) < min(0, log_acceptance_ratio)) {
      
      R <- R_prop
      sigma <- sigma_prop
      accept_R <- accept_R + 1
      
    }
  }
  
  return(list(
    R = R,
    sigma = sigma,
    accept_R = accept_R
  ))
}









compare_mnp_mcmc <- function(fit, X, y, true_beta, true_sigma,
                             r = 250, Halton_choice = FALSE,
                             seed = 12345) {
  
  X <- as.matrix(X)
  k <- length(unique(y))
  p <- ncol(X)
  n <- nrow(X)
  
  category_position <- lapply(
    1:k, function(j) which(y == j)
  )
  
  # Automatically extract MCMC output
  beta <- fit$beta_final_vector
  sigma <- fit$sigma2_final
  
  # Observation-level log-likelihood
  get_loglik <- function(beta_vector, sigma_matrix) {
    
    X_BETA <- X %*% matrix(
      beta_vector, nrow = p, ncol = k - 1
    )
    
    log_likelihood_mnp_1response_WAIC(
      category_position,
      mu = X_BETA,
      sigma = sigma_matrix,
      r = r,
      Halton_choice = Halton_choice
    )$observation_log_likelihood
  }
  
  set.seed(seed)
  
  # True log-likelihood
  True_loglik <- sum(get_loglik(true_beta, true_sigma))
  
  # Posterior pointwise log-likelihood
  loglikelihood_WAIC <- t(
    vapply(seq_len(ncol(beta)), function(i) {
      
      get_loglik(beta[, i], sigma[, , i])
      
    }, numeric(n))
  )
  
  # Total log-likelihood
  loglikelihood <- rowSums(loglikelihood_WAIC)
  
  # Posterior mean parameters
  beta_bar <- rowMeans(beta)
  sigma_bar <- apply(sigma, c(1, 2), mean)
  
  # DIC and BPIC
  D_bar <- -2 * mean(loglikelihood)
  
  D_theta_bar <- -2 * sum(
    get_loglik(beta_bar, sigma_bar)
  )
  
  p_D <- D_bar - D_theta_bar
  
  DIC <- D_bar + p_D
  BPIC <- D_bar + 2 * p_D
  
  # 95% HPD interval
  chain <- coda::mcmc(loglikelihood)
  
  HPDI <- coda::HPDinterval(
    chain, prob = 0.95
  )[1, ]
  
  # WAIC
  WAIC <- loo::waic(
    loglikelihood_WAIC
  )$estimates["waic", "Estimate"]
  
  # Geweke Z-score
  Geweke <- coda::geweke.diag(
    chain, frac1 = 0.2, frac2 = 0.4
  )$z
  
  # Effective sample size
  ESS <- coda::effectiveSize(chain)
  
  # Final summary
  summary_table <- data.frame(
    True = round(True_loglik, 2),
    HPDI = sprintf("(%.2f, %.2f)", HPDI[1], HPDI[2]),
    BPIC = round(BPIC, 2),
    WAIC = round(as.numeric(WAIC), 2),
    Geweke_Z = round(as.numeric(Geweke), 3),
    ESS = round(as.numeric(ESS), 2)
  )
  
  return(summary_table)
}










# ==========================================================
# Function to calculate TVD, KL, and JSD
# Using direct philentropy functions
# ==========================================================

calculate_distances <- function(estimated_PP_list, true_PP_list) {
  
  stopifnot(length(estimated_PP_list) == length(true_PP_list))
  
  n_cv <- length(estimated_PP_list)
  
  points_TVD <- numeric(n_cv)
  points_KL  <- numeric(n_cv)
  points_JSD <- numeric(n_cv)
  
  for (j in seq_len(n_cv)) {
    
    estimated <- estimated_PP_list[[j]]
    P <- as.matrix(true_PP_list[[j]])
    
    # Dimensions:
    # observations x categories x posterior samples
    
    stopifnot(length(dim(estimated)) == 3)
    
    n_obs <- dim(estimated)[1]
    n_categories <- dim(estimated)[2]
    n_samples <- dim(estimated)[3]
    
    stopifnot(
      nrow(P) == n_obs,
      ncol(P) == n_categories
    )
    
    # Store distances for each observation and posterior sample
    out_TVD <- matrix(NA_real_, n_obs, n_samples)
    out_KL  <- matrix(NA_real_, n_obs, n_samples)
    out_JSD <- matrix(NA_real_, n_obs, n_samples)
    
    for (s in seq_len(n_samples)) {
      
      Q <- matrix(
        estimated[, , s],
        nrow = n_obs,
        ncol = n_categories
      )
      
      for (i in seq_len(n_obs)) {
        
        # P = true probability distribution
        # Q = estimated probability distribution
        
        p <- as.numeric(P[i, ])
        q <- as.numeric(Q[i, ])
        
        # --------------------------------------------------
        # Total Variation Distance
        # --------------------------------------------------
        
        out_TVD[i, s] <- philentropy::manhattan(
          P = p,
          Q = q,
          testNA = FALSE
        ) / 2
        
        # --------------------------------------------------
        # KL Divergence: KL(P || Q)
        # Natural logarithm
        # --------------------------------------------------
        
        out_KL[i, s] <- philentropy::kullback_leibler_distance(
          P = p,
          Q = q,
          unit = "log",
          epsilon = 1e-5,
          testNA = FALSE
        )
        
        # --------------------------------------------------
        # Jensen-Shannon Divergence
        # Natural logarithm
        # --------------------------------------------------
        
        out_JSD[i, s] <- philentropy::jensen_shannon(
          P = p,
          Q = q,
          unit = "log",
          testNA = FALSE
        )
        
      }
    }
    
    # Average over observations and posterior samples
    
    points_TVD[j] <- mean(out_TVD)
    points_KL[j]  <- mean(out_KL)
    points_JSD[j] <- mean(out_JSD)
    
    cat("Completed CV:", j, "of", n_cv, "\n")
  }
  
  return(list(
    TVD = points_TVD,
    KL  = points_KL,
    JSD = points_JSD
  ))
}





















compare_bivariate_mnp_mcmc <- function(fit, X1, X2, y1, y2, true_beta, true_sigma,
                                       k1, k2, r = 250, Halton_choice = FALSE,
                                       seed = 12345) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  n <- nrow(X1)
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  stopifnot(nrow(X2) == n, length(y1) == n, length(y2) == n)
  
  beta <- fit$beta_final_vector
  sigma <- fit$sigma2_final
  
  big_x1 <- diag(k1-1) %x% X1
  big_x2 <- diag(k2-1) %x% X2
  
  x_tilda <- rbind(
    cbind(big_x1, matrix(0,nrow=nrow(big_x1),ncol=ncol(big_x2))),
    cbind(matrix(0,nrow=nrow(big_x2),ncol=ncol(big_x1)),big_x2)
  )
  
  y12 <- cbind(y1,y2)
  category_position_y12 <- as.list(rep(NA,k1*k2))
  dim(category_position_y12) <- c(k1,k2)
  
  for(i in 1:k1){
    for(j in 1:k2){
      category_position_y12[[i,j]] <- which(y1 == i & y2 == j)
    }
  }
  
  get_X_BETA <- function(beta_vector){
    X_BETA_tilda <- colSums(t(x_tilda)*c(beta_vector))
    matrix(X_BETA_tilda,ncol=(k1-1)+(k2-1))
  }
  
  get_loglik <- function(beta_vector,sigma_matrix){
    
    X_BETA <- get_X_BETA(beta_vector)
    
    log_likelihood_mnp_2response_WAIC(
      pos_y = category_position_y12,
      k1 = k1,
      k2 = k2,
      mu1 = X_BETA[,1:(k1-1),drop=FALSE],
      mu2 = X_BETA[,k1:(k1+k2-2),drop=FALSE],
      sigma = sigma_matrix,
      r = r,
      Halton_choice = Halton_choice
    )$observation_log_likelihood
  }
  
  set.seed(seed)
  
  # True log-likelihood
  True_loglik <- sum(get_loglik(true_beta,true_sigma))
  
  # Posterior pointwise log-likelihood
  n_draws <- ncol(beta)
  
  loglikelihood_WAIC <- t(
    vapply(seq_len(n_draws),function(i){
      get_loglik(beta[,i],sigma[,,i])
    },numeric(n))
  )
  
  # Total log-likelihood
  loglikelihood <- rowSums(loglikelihood_WAIC)
  
  # Posterior mean parameters
  beta_bar <- rowMeans(beta)
  sigma_bar <- apply(sigma,c(1,2),mean)
  
  # DIC and BPIC
  D_bar <- -2*mean(loglikelihood)
  D_theta_bar <- -2*sum(get_loglik(beta_bar,sigma_bar))
  
  p_D <- D_bar-D_theta_bar
  
  DIC <- D_bar+p_D
  BPIC <- D_bar+2*p_D
  
  # 95% HPD interval
  chain <- coda::mcmc(loglikelihood)
  HPDI <- coda::HPDinterval(chain,prob=0.95)[1,]
  
  # WAIC
  WAIC <- loo::waic(loglikelihood_WAIC)$estimates["waic","Estimate"]
  
  # Geweke Z-score
  Geweke <- coda::geweke.diag(
    chain,frac1=0.2,frac2=0.4
  )$z
  
  # Effective sample size
  ESS <- coda::effectiveSize(chain)
  
  # Final summary
  summary_table <- data.frame(
    True = round(True_loglik,2),
    HPDI = sprintf("(%.2f, %.2f)",HPDI[1],HPDI[2]),
    DIC = round(DIC,2),
    BPIC = round(BPIC,2),
    WAIC = round(as.numeric(WAIC),2),
    Geweke_Z = round(as.numeric(Geweke),3),
    ESS = round(as.numeric(ESS),2)
  )
  
  return(summary_table)
}
