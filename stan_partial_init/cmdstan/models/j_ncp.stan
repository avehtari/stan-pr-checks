data {
  int<lower=0> N;
  int<lower=0> n_grp;
  vector[N] y;
  array[N] int<lower=1, upper=n_grp> grp_id;
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> sigma_residual;
  vector[n_grp] raw;
}
model {
  raw ~ std_normal();
  y ~ normal(mu + sigma * raw[grp_id], sigma_residual);
}
