// Cases for the stanc3 offset/multiplier std_normal prototype.
data {
  int<lower=1> N;
  real<lower=0> s_data;
  vector<lower=0>[N] s_vec_data;
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector<lower=0>[N] sigma_vec;
  vector<offset=mu, multiplier=sigma>[N] a;         // tilde: rewrite
  vector<multiplier=sigma>[N] b;                    // multiplier only: rewrite
  vector<offset=mu>[N] c;                           // offset only: rewrite
  real<offset=mu, multiplier=sigma> d;              // scalar, normal_lpdf: rewrite
  vector<offset=mu, multiplier=sigma_vec>[N] e;     // vector multiplier: rewrite
  vector<offset=mu, multiplier=s_data>[N] f;        // data multiplier, lupdf: rewrite
  vector<offset=mu, multiplier=s_vec_data>[N] g;    // data vector multiplier, lpdf
  array[N] real<offset=mu, multiplier=sigma> h;     // array: rewrite
  vector<offset=mu, multiplier=sigma>[N] k;         // other scale: no rewrite
  matrix<offset=mu, multiplier=sigma>[2, 2] m;      // matrix: no rewrite
}
model {
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
  sigma_vec ~ normal(0, 1);
  a ~ normal(mu, sigma);
  b ~ normal(0, sigma);
  c ~ normal(mu, 1);
  target += normal_lpdf(d | mu, sigma);
  e ~ normal(mu, sigma_vec);
  f ~ normal(mu, s_data);
  target += normal_lpdf(g | mu, s_vec_data);
  h ~ normal(mu, sigma);
  k ~ normal(mu, 2 * sigma);
  to_vector(m) ~ normal(mu, sigma);
}
