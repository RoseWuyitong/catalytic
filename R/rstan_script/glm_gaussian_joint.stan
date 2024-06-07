data {
  int<lower=1> obs_n;
  int<lower=1> syn_n;
  int<lower=1> col_num;
  matrix[obs_n, col_num] obs_x;
  matrix[syn_n, col_num] syn_x;
  vector[obs_n] obs_y;
  vector[syn_n] syn_y;
  real kappa;
  real alpha;
  real gamma;
}

parameters {
  real<lower=0> tau;
  vector[col_num+1] beta;
  real<lower=0> sigma;
}

model {
  vector[obs_n] obs_mu = beta[1] + obs_x * beta[2:(col_num+1)];
  vector[syn_n] syn_mu = beta[1] + syn_x * beta[2:(col_num+1)];
  real syn_weight = tau/syn_n;

  // Likelihood
  for (i in 1:obs_n)
    target += normal_lpdf(obs_y[i] | obs_mu[i], sigma);

  // Catalytic prior
  for (i in 1:syn_n)
      target += syn_weight * normal_lpdf(syn_y[i] | syn_mu[i], sigma);

  // Gamma function for tau
  target += ((col_num+alpha)/2-1)*log(tau)-tau*(kappa+1/gamma);
}
