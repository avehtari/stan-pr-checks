parameters {
  real<lower=0> sigma;
  real<multiplier=sigma> x;
}
model {
  sigma ~ std_normal();
  x ~ normal(0, sigma);
}
