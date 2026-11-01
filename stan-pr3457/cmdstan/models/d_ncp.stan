data {
  int<lower=0> I;
  array[I] int N;
  array[I] int y;
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector[I] raw;
}
model {
  mu ~ student_t(3, 0, 2.5);
  sigma ~ normal(0, 1);
  raw ~ std_normal();
  y ~ binomial_logit(N, mu + sigma * raw);
}
