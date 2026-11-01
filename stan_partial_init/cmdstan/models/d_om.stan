data {
  int<lower=0> I;
  array[I] int N;
  array[I] int y;
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector<offset=mu, multiplier=sigma>[I] logit_p;
}
model {
  mu ~ student_t(3, 0, 2.5);
  sigma ~ normal(0, 1);
  logit_p ~ normal(mu, sigma);
  y ~ binomial_logit(N, logit_p);
}
