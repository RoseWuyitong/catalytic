data {
  int<lower=1> obs_n;
  int<lower=1> syn_n;
  int<lower=1> x_col_n;
  matrix[obs_n, x_col_n] obs_x;
  matrix[syn_n, x_col_n] syn_x;
  vector[obs_n] obs_y;
  vector[syn_n] syn_y;
  real<lower=0> tau;
}

parameters {
  vector[x_col_n+1] beta;
}

model {
  vector[obs_n] obs_mu = beta[1] + obs_x * beta[2:(x_col_n+1)];
  vector[syn_n] syn_mu = beta[1] + syn_x * beta[2:(x_col_n+1)];
  real obs_sigma = pow(sd(obs_mu), 2);
  real syn_sigma = pow(sd(syn_mu), 2);

  // Likelihood
  for (i in 1:obs_n)
    target += normal_lpdf(obs_y[i] | obs_mu[i], obs_sigma);

  // Catalytic prior
  for (i in 1:syn_n)
      target += tau/syn_n * normal_lpdf(syn_y[i] | syn_mu[i], syn_sigma);
}
