data {
  int<lower=0> obs_n;
  int<lower=0> syn_n;
  int<lower=0> x_col_n;
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
  vector[N] obs_mu = beta[1] + obs_x * beta[2:x_col_n+1];
  vector[M] syn_mu = beta[1] + syn_x * beta[2:x_col_n+1];

  // Likelihood
  for (i in 1:N)
    target += normal_lpdf(obs_y[i] | obs_mu[i]);

  // Catalytic prior
  for (i in 1:M)
      target += tau/syn_n * normal_lpdf(syn_y[i] | syn_mu[i]);
}
