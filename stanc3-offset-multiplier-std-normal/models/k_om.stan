transformed data {
  int n = 200;
}
parameters {
  real m;
  real<lower=0> s;
  vector<offset=m, multiplier=s>[n] pr;
}
model {
  m ~ normal(-0.55, 1.0);
  s ~ normal(0, 0.5);
  pr ~ normal(m, s);
}
