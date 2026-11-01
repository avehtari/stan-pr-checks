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
  vector[n_grp] re = mu + sigma * raw;
  re ~ normal(mu, sigma);
  target += log(sigma) * n_grp;
  y ~ normal(re[grp_id], sigma_residual);
}
