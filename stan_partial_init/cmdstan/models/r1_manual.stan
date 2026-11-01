parameters {
  real<lower=0> sigma;
  real x_raw;
}
transformed parameters {
  real x = x_raw * sigma;
}
model {
  sigma ~ std_normal();
  x ~ normal(0, sigma);
  target += log(sigma);
}
