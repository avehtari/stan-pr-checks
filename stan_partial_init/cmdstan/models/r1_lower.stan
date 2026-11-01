parameters {
  real a;
  real<lower=a> z;
}
model {
  a ~ normal(0, 5);
  z ~ normal(a + 1, 1);
}
