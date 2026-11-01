data {
  int n;
  vector[n] y;
}
parameters {
  real<lower=0> sigma;
  vector<multiplier=sigma>[n] x;
}
model {
  sigma ~ std_normal();
  x ~ normal(0, sigma);
  y ~ normal(x, 0.01);
}
