



data {
  int<lower=0> J;  // Number of time intervals.
  vector[J] hPriorSh;  // Shape parameters for the gamma prior distribution of the baseline hazard.
  real c0;  // Rate parameter for the gamma prior distribution of the baseline hazard.
  real eta0;  // A constant used for calculating the prior distribution of the baseline hazard.
  real kappa0;  // A constant used for calculating the prior distribution of the baseline hazard.
  
  int<lower=0> n;  // Number of observations in the dataset.
  int<lower=0> p;  // Dimensionality of the covariates.
  matrix[n, p] X;  // Matrix of covariates for the observed dataset.
  
  
  real lambda;  // gaussian prior scale
  
  matrix[n, J] R_tilde_minus_D_tilde;  // Matrix indicating risk set minus event set for each observation across intervals.
  matrix[n, J] D_tilde;  // Matrix indicating which intervals an observation has an event.
}


parameters {
  vector[p] beta;  // regression coefficients
  //real tau_downweight;  // Downweighting factor for the synthetic dataset in the estimation process.
  vector<lower=0>[J] h_seq;  // parameter with a gamma prior
}

model {
  /**
   * This model block represents a Bayesian Cox Proportional Hazards model with the inclusion of a customized 
   * hazard function and a log catalytic prior for the coefficients based on synthetic data.
   *
   * The likelihood component of the model is constructed using observed data, while the prior is constructed using 
   * both observed and synthetic data.
   *
   * The total target log posterior (unnormalized) is updated with contributions from both the likelihood and the prior.
   */

  matrix[n, J] exp_xbeta_mat = rep_matrix(exp(X * beta), J);  // Matrix where each column is the exponential of X multiplied by beta.
  vector[J] first_sum;  // Vector to store the summation terms for the risk set minus event set.
  matrix[n, J] h_mat = rep_matrix(h_seq', n);  // Replicating the hazard sequence across `n` rows.
  matrix[n, J] h_exp_xbeta_mat = -h_mat .* exp_xbeta_mat;  // Matrix storing product of hazard sequence and the exponential transformation of X and beta.
  vector[J] second_sum;  // Vector to store the summation terms for the event set.
  real logprior =  lambda*dot_self(beta)/2;
  
  
  for (j in 1:J) {
    first_sum[j] = sum(exp_xbeta_mat[, j] .* R_tilde_minus_D_tilde[, j]);  // Summing over the risk set minus event set for the `j-th` interval.
    second_sum[j] = sum(log1m_exp(h_exp_xbeta_mat[, j]) .* D_tilde[, j]);  // Summing over the event set for the `j-th` interval using the log1m_exp transformation.
  }
  
  target += sum(-h_seq .* first_sum + second_sum)-logprior;  // Update the target log posterior with the likelihood component.
  //target = target-1/2*lambda*dot_self(beta);  // Update the target log posterior with the log prior component.
  
  
  h_seq ~ gamma(hPriorSh, c0);  // Prior distribution for the hazard sequence based on a Gamma distribution.
}


