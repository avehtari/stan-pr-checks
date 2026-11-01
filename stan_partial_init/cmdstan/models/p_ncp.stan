data {
  int n;
  vector[n] y;
}
parameters {
  real<lower=0> sigma;
  vector[n] raw;
}
model {
  sigma ~ std_normal();
  raw ~ std_normal();
  y ~ normal(sigma * raw, 0.01);
}
