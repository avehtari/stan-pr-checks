parameters {
  real mu;
  real<lower=0> sigma;
  real<offset=mu, multiplier=sigma> x;
}
model {
  mu ~ std_normal();
  sigma ~ std_normal();
  x ~ normal(mu, sigma);
}
