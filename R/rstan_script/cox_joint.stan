functions {
  /**
   * Calculates the log prior for the coefficients (`coefs`) in a Cox model with a customized hazard constant.
   *
   * This function computes the log prior using the synthetic data and a specified hazard constant.
   * The contribution from each synthetic observation is calculated, and then scaled by a downweighting factor.
   *
   * @param syn_size Number of synthetic samples.
   * @param coefs Vector of model coefficients.
   * @param syn_x Matrix of synthetic covariates.
   * @param syn_time Vector of synthetic survival times.
   * @param syn_event Vector of synthetic event indicators.
   * @param syn_time_hazard_constant Customized hazard constant values evaluated at each synthetic survival time.
   * @param tau Positive number that downweighting the synthetic dataset in the estimation process.
   *
   * @return A scalar representing the log prior for the given `coefs` using the synthetic data.
   */
  real log_of_cat_prior(vector coefs,
                        real tau,
                        int syn_size,
                        matrix syn_x,
                        vector syn_time,
                        vector syn_event,
                        real hazard_constant
                        ) {

    vector[syn_size] syn_lp = syn_x * coefs;
    vector[syn_size] exp_syn_lp = exp(syn_lp);

    real syn_prior = 0;
    for (i in 1:syn_size) {
      syn_prior = syn_prior + (syn_event[i] * syn_lp[i] - exp_syn_lp[i] * syn_time[i] * hazard_constant);
    }
    return syn_prior * tau / syn_size;
  }
}


data {
  int<lower=0> J;  // Number of time intervals.
  vector[J] alpha_of_h_j_list;  // Shape parameters for the gamma prior distribution of the baseline hazard.
  real c_0;  // Rate parameter for the gamma prior distribution of the baseline hazard.

  int<lower=0> obs_size;  // Number of observations in the dataset.
  int<lower=0> col_num;  // Dimensionality of the covariates.
  matrix[obs_size, col_num] obs_x;  // Matrix of covariates for the observed dataset.

  int<lower=0> syn_size;  // Number of synthetic observations.
  matrix[syn_size, col_num] syn_x;  // Matrix of covariates for the synthetic dataset.
  vector[syn_size] syn_time;  // Vector of synthetic survival times.
  vector[syn_size] syn_event;  // Vector of censoring indicator.
  real<lower=0> hazard_constant;  // Customized hazard constant value.

  matrix[obs_size, J] adj_risk_set;  // Matrix indicating risk set minus event set for each observation across intervals.
  matrix[obs_size, J] failure_set;  // Matrix indicating which intervals an observation has an event.

  real alpha;
  real gamma;
  real kappa;

}


parameters {
  vector[col_num] coefs;  // regression coefficients
  vector<lower=0>[J] h_j;  // parameter with a gamma prior
  real<lower=0> tau; // downweighting factor for synthetic data
}

model {
  /**
   * This model block represents a Bayesian Cox Proportional Hazards model with the inclusion of a customized
   * hazard constant and a log catalytic prior for the coefficients based on synthetic data.
   *
   * The likelihood component of the model is constructed using observed data, while the prior is constructed using
   * both observed and synthetic data.
   *
   * The total target log posterior (unnormalized) is updated with contributions from both the likelihood and the prior.
   */
  vector[J] adj_risk_set_sum;  // Vector to store the summation terms for the risk set minus event set.
  vector[J] event_set_sum;  // Vector to store the summation terms for the event set.
  matrix[obs_size, J] exp_obs_lp = rep_matrix(exp(obs_x * coefs), J);  // Matrix where each column is the exponential of obs_x multiplied by coefs.
  matrix[obs_size, J] h_seq = rep_matrix(h_j', obs_size);  // Replicating the hazard sequence across `obs_size` rows.
  matrix[obs_size, J] h_seq_exp_obs_lp = -h_seq .* exp_obs_lp;  // Matrix storing product of hazard sequence and the exponential transformation of obs_x and coefs.

  for (j in 1:J) {
    // Summing over the adjusted risk set for the `j-th` interval.
    adj_risk_set_sum[j] = sum(exp_obs_lp[, j] .* adj_risk_set[, j]);
    // Summing over the event set for the `j-th` interval using the log1m_exp transformation.
    event_set_sum[j] = sum(log1m_exp(h_seq_exp_obs_lp[, j]) .* failure_set[, j]);
  }

  // Update the target log posterior with the likelihood component.
  target += sum(-h_j .* adj_risk_set_sum + event_set_sum);

  // Update the target log posterior with the log of catalytic prior component.
  target += log_of_cat_prior(coefs,
                             tau,
                             syn_size,
                             syn_x,
                             syn_time,
                             syn_event,
                             hazard_constant);

  // Update the target log posterior with the log of tau likelihood
  target += (col_num + alpha - 1) * log(tau) - (kappa + 1/gamma) * tau;

  // Prior distribution for the hazard sequence based on a Gamma distribution.
  h_j ~ gamma(alpha_of_h_j_list, c_0);
}


