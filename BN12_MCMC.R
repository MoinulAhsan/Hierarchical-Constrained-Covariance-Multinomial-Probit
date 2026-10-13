# BN12: Burgette and Nordheim (2012) trace-constrained MNP MCMC
# Required packages: mvtnorm (rmvnorm), LaplacesDemon (rinvwishart).
# X must include an intercept if the model requires one; y must be 1,...,k.

run_multinomial_univariate_mcmc_BN12 <- function(X, y, nit, nmc, nthin,
                                                 BN12_sigma_beta2 = 10000,
                                                 BN12_df = NULL) {
  X <- as.matrix(X)
  y <- as.integer(y)
  BN12_n <- nrow(X)
  BN12_p <- ncol(X)
  BN12_k <- length(unique(y))
  BN12_q <- BN12_k - 1L
  
  if (!is.numeric(X) || anyNA(X) || anyNA(y) || length(y) != BN12_n)
    stop("BN12 requires a numeric X without NAs and a response y of matching length.")
  if (BN12_k < 3L || !identical(sort(unique(y)), seq_len(BN12_k)))
    stop("BN12 requires at least three response categories coded 1,...,k.")
  if (any(!is.finite(c(nit, nmc, nthin))) || nit < 0 || nmc < 1 ||
      nthin < 1 || any(c(nit, nmc, nthin) != as.integer(c(nit, nmc, nthin))) ||
      nmc %% nthin != 0)
    stop("BN12 requires nit >= 0, nmc >= 1 and nthin >= 1 as integers, with nmc divisible by nthin.")
  if (is.null(BN12_df)) BN12_df <- BN12_k + 1
  if (!is.finite(BN12_sigma_beta2) || BN12_sigma_beta2 <= 0 ||
      !is.finite(BN12_df) || BN12_df <= BN12_q - 1)
    stop("Invalid BN12 prior hyperparameters.")
  if (!requireNamespace("mvtnorm", quietly = TRUE) ||
      !requireNamespace("LaplacesDemon", quietly = TRUE))
    stop("Install mvtnorm and LaplacesDemon before running BN12.")
  
  BN12_inv <- function(M) chol2inv(chol(M))
  BN12_tr <- function(M) sum(diag(M))
  BN12_XTX <- crossprod(X)
  BN12_x_tilda <- diag(BN12_q) %x% X
  BN12_positions <- lapply(seq_len(BN12_k), function(j) which(y == j))
  BN12_frequency <- lengths(BN12_positions)
  
  BN12_mu_beta <- matrix(0, BN12_p * BN12_q, 1)
  BN12_inv_sigma_beta <- diag(1 / BN12_sigma_beta2, BN12_p * BN12_q)
  
  BN12_beta <- matrix(0, BN12_p * BN12_q, 1)
  BN12_sigma <- diag(BN12_q)
  BN12_sigma_inv <- BN12_inv(BN12_sigma)
  BN12_mean_w <- matrix(BN12_x_tilda %*% BN12_beta, ncol = BN12_q)
  BN12_w_mat <- cbind(0, matrix(NA_real_, BN12_n, BN12_q))
  
  # Initialize Gaussian latent utilities.
  BN12_m1 <- c(BN12_mean_w[BN12_positions[[1]], , drop = FALSE])
  BN12_v1 <- sqrt(rep(diag(BN12_sigma), each = BN12_frequency[1]))
  BN12_w_mat[BN12_positions[[1]], -1] <- qnorm(
    pnorm(-Inf, BN12_m1, BN12_v1) +
      runif(BN12_frequency[1] * BN12_q) *
      (pnorm(0, BN12_m1, BN12_v1) - pnorm(-Inf, BN12_m1, BN12_v1)),
    BN12_m1, BN12_v1)
  
  for (i in 2:BN12_k) {
    BN12_idx <- BN12_positions[[i]]
    BN12_m <- BN12_mean_w[BN12_idx, i - 1]
    BN12_v <- sqrt(BN12_sigma[i - 1, i - 1])
    BN12_w_mat[BN12_idx, i] <- qnorm(
      pnorm(0, BN12_m, BN12_v) + runif(length(BN12_idx)) *
        (1 - pnorm(0, BN12_m, BN12_v)), BN12_m, BN12_v)
    
    BN12_cond_m <- c(BN12_mean_w[BN12_idx, -(i - 1), drop = FALSE]) +
      c(t(BN12_sigma[-(i - 1), i - 1, drop = FALSE] %*%
            BN12_inv(BN12_sigma[i - 1, i - 1, drop = FALSE]) %*%
            t(matrix(BN12_w_mat[BN12_idx, i] - BN12_m))))
    BN12_cond_cov <- BN12_sigma[-(i - 1), -(i - 1), drop = FALSE] -
      BN12_sigma[-(i - 1), i - 1, drop = FALSE] %*%
      BN12_inv(BN12_sigma[i - 1, i - 1, drop = FALSE]) %*%
      BN12_sigma[i - 1, -(i - 1), drop = FALSE]
    BN12_cond_sd <- sqrt(matrix(diag(BN12_cond_cov), nrow = length(BN12_idx),
                                ncol = BN12_q - 1L, byrow = TRUE))
    BN12_bound <- BN12_w_mat[BN12_idx, i]
    BN12_w_mat[BN12_idx, -c(1, i)] <- qnorm(
      runif(length(BN12_idx) * (BN12_q - 1L)) *
        pnorm(BN12_bound, BN12_cond_m, BN12_cond_sd),
      BN12_cond_m, BN12_cond_sd)
  }
  
  BN12_nsamples <- nmc / nthin
  BN12_total_iteration <- nit + nmc
  BN12_beta_final <- array(NA_real_, c(BN12_q, BN12_p, BN12_nsamples))
  BN12_Trace_beta <- matrix(NA_real_, BN12_p * BN12_q, BN12_nsamples)
  BN12_Trace_sigma <- array(NA_real_, c(BN12_q, BN12_q, BN12_nsamples))
  BN12_Trace_beta_our_scale <- matrix(NA_real_, BN12_p * BN12_q, BN12_nsamples)
  BN12_Trace_sigma_our_scale <- array(NA_real_, c(BN12_q, BN12_q, BN12_nsamples))
  
  BN12_start <- proc.time()[["elapsed"]]
  for (iter in seq_len(BN12_total_iteration)) {
    BN12_X_BETA <- X %*% matrix(BN12_beta, nrow = BN12_p, ncol = BN12_q)
    
    # 1. Update latent utilities, conditional on the other utilities.
    for (i in 2:BN12_k) {
      BN12_j <- i - 1
      BN12_other <- setdiff(seq_len(BN12_q), BN12_j)
      BN12_cond_mean <- c(BN12_X_BETA[, BN12_j] +
                            t(BN12_sigma[BN12_j, BN12_other, drop = FALSE] %*%
                                BN12_inv(BN12_sigma[BN12_other, BN12_other, drop = FALSE]) %*%
                                t(BN12_w_mat[, BN12_other + 1, drop = FALSE] -
                                    BN12_X_BETA[, BN12_other, drop = FALSE])))
      BN12_cond_var <- c(BN12_sigma[BN12_j, BN12_j] -
                           BN12_sigma[BN12_j, BN12_other, drop = FALSE] %*%
                           BN12_inv(BN12_sigma[BN12_other, BN12_other, drop = FALSE]) %*%
                           BN12_sigma[BN12_other, BN12_j, drop = FALSE])
      BN12_sd <- sqrt(BN12_cond_var)
      BN12_idx <- BN12_positions[[i]]
      BN12_not_idx <- which(y != i)
      
      BN12_lower <- apply(BN12_w_mat[BN12_idx, -i, drop = FALSE], 1, max)
      BN12_w_mat[BN12_idx, i] <- qnorm(
        pnorm(BN12_lower, BN12_cond_mean[BN12_idx], BN12_sd) +
          runif(length(BN12_idx)) *
          (1 - pnorm(BN12_lower, BN12_cond_mean[BN12_idx], BN12_sd)),
        BN12_cond_mean[BN12_idx], BN12_sd)
      
      BN12_upper <- apply(BN12_w_mat[BN12_not_idx, , drop = FALSE], 1, max)
      BN12_w_mat[BN12_not_idx, i] <- qnorm(
        runif(length(BN12_not_idx)) *
          pnorm(BN12_upper, BN12_cond_mean[BN12_not_idx], BN12_sd),
        BN12_cond_mean[BN12_not_idx], BN12_sd)
    }
    
    # 2. BN12 working scale alpha and Gaussian beta update.
    BN12_alpha <- sqrt(BN12_tr(BN12_sigma_inv) /
                         rchisq(1, BN12_df * BN12_q))
    BN12_w_mat1 <- BN12_w_mat * BN12_alpha
    BN12_w_tilda <- matrix(c(BN12_w_mat[, -1, drop = FALSE]) * BN12_alpha,
                           ncol = 1)
    
    BN12_A_beta <- BN12_sigma_inv %x% BN12_XTX + BN12_inv_sigma_beta
    BN12_B_beta <- (BN12_sigma_inv %x% t(X)) %*% BN12_w_tilda +
      BN12_inv_sigma_beta %*% BN12_mu_beta
    BN12_var1 <- BN12_inv(BN12_A_beta)
    BN12_mu1 <- BN12_var1 %*% BN12_B_beta
    
    BN12_resid <- matrix(c(BN12_w_tilda - BN12_x_tilda %*% BN12_mu1),
                         ncol = BN12_q)
    BN12_aa <- sum((BN12_resid %*% BN12_sigma_inv) * BN12_resid)
    BN12_alpha <- sqrt(
      (BN12_aa + c(t(BN12_mu1) %*% BN12_inv_sigma_beta %*% BN12_mu1) +
         BN12_tr(BN12_sigma_inv)) /
        rchisq(1, (BN12_df + BN12_n) * BN12_q))
    BN12_beta_tilda <- matrix(mvtnorm::rmvnorm(1, c(BN12_mu1),
                                               BN12_alpha^2 * BN12_var1), ncol = 1)
    
    # 3. BN12 inverse-Wishart covariance update under the trace constraint.
    BN12_r <- matrix(c(BN12_w_tilda - BN12_x_tilda %*% BN12_beta_tilda),
                     ncol = BN12_q)
    BN12_A <- crossprod(BN12_r)
    BN12_sigma_tilda <- LaplacesDemon::rinvwishart(
      BN12_df + BN12_n, BN12_A + diag(BN12_q))
    BN12_alpha <- sqrt(BN12_tr(BN12_sigma_tilda) / BN12_q)
    BN12_sigma <- BN12_sigma_tilda / BN12_alpha^2
    BN12_sigma_inv <- BN12_inv(BN12_sigma)
    BN12_beta <- BN12_beta_tilda / BN12_alpha
    BN12_w_mat <- BN12_w_mat1 / BN12_alpha
    
    # 4. Retain posterior draws.
    if (iter > nit && (iter - nit) %% nthin == 0L) {
      BN12_t <- (iter - nit) / nthin
      BN12_beta_final[, , BN12_t] <- matrix(BN12_beta, nrow = BN12_q, byrow = TRUE)
      BN12_Trace_beta[, BN12_t] <- c(BN12_beta)
      BN12_Trace_sigma[, , BN12_t] <- BN12_sigma
      
      # Rescale to zero-sum log-S constraint for parameter comparisons.
      BN12_c_m <- prod(diag(BN12_sigma))^(1 / BN12_q)
      BN12_Trace_beta_our_scale[, BN12_t] <- c(BN12_beta) / sqrt(BN12_c_m)
      BN12_Trace_sigma_our_scale[, , BN12_t] <- BN12_sigma / BN12_c_m
    }
  }
  
  BN12_result <- list(
    BN12_beta_final = BN12_beta_final,
    BN12_Trace_beta = BN12_Trace_beta,
    BN12_Trace_sigma = BN12_Trace_sigma,
    BN12_Trace_beta_our_scale = BN12_Trace_beta_our_scale,
    BN12_Trace_sigma_our_scale = BN12_Trace_sigma_our_scale,
    beta_final_vector = BN12_Trace_beta,
    sigma2_final = BN12_Trace_sigma,
    BN12_runtime_seconds = proc.time()[["elapsed"]] - BN12_start
  )
  return(BN12_result)
}
