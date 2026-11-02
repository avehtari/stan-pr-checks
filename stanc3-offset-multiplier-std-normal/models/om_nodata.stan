// Cases for the stanc3 offset/multiplier std_normal prototype, without a
// data block (for the C++ probe om_probe_main.hpp).
functions {
  real f(real s) {
    return s;
  }
}
transformed data {
  int N = 3;
  real s_data = 0.7;
  vector[3] s_vec_data = [0.5, 1, 2]';
  int K = 2;
}
parameters {
  real mu;
  real<lower=0> sigma;
  vector<lower=0>[N] sigma_vec;
  vector<offset=mu, multiplier=sigma>[N] a;        // tilde: rewrite
  vector<multiplier=sigma>[N] b;                   // multiplier only
  vector<offset=mu>[N] c;                          // offset only
  real<offset=mu, multiplier=sigma> d;             // normal_lpdf, scalar
  vector<offset=mu, multiplier=sigma_vec>[N] e;    // vector multiplier
  vector<offset=mu, multiplier=s_data>[N] g;       // data multiplier, tilde
  real<multiplier=s_data> g2;                      // data multiplier, lupdf
  vector<offset=mu, multiplier=s_vec_data>[N] g3;  // data vector multiplier
  array[N] real<multiplier=K> h;                   // int multiplier
  real<offset=mu, multiplier=f(sigma)> p;          // user function, profile
  vector<multiplier=sigma>[N] q;                   // in a loop and an if
  real<offset=mu, multiplier=sigma> t;             // truncation
  real<multiplier=sigma> z;                        // two statements
  tuple(real<multiplier=sigma>, real) tp;          // tuple: no rewrite
  vector<offset=mu, multiplier=sigma>[N] k;        // other scale: no rewrite
  matrix<offset=mu, multiplier=sigma>[2, 2] m;     // matrix: no rewrite
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
  g ~ normal(mu, s_data);
  target += normal_lupdf(g2 | 0, s_data);
  target += normal_lpdf(g3 | mu, s_vec_data);
  h ~ normal(0, K);
  profile("p") {
    p ~ normal(mu, f(sigma));
  }
  for (i in 1:N) {
    if (i > 1) {
      q ~ normal(0, sigma);
    }
    q[i] ~ normal(0, sigma);
  }
  t ~ normal(mu, sigma) T[-1000, ];
  z ~ normal(0, sigma);
  z ~ normal(0, sigma);
  tp.1 ~ normal(0, sigma);
  tp.2 ~ normal(0, 1);
  k ~ normal(mu, 2 * sigma);
  to_vector(m) ~ normal(mu, sigma);
}
