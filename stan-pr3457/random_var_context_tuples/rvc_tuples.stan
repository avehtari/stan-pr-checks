parameters {
  real a;
  array[2] tuple(real<lower=10>, real<upper=0>) at;
}
model {
  a ~ normal(0, 1);
  for (i in 1:2) {
    at[i].1 ~ normal(11, 1);
    at[i].2 ~ normal(-1, 1);
  }
}
