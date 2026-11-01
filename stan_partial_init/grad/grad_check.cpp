// Gradient of log p(mu, log_sigma, raw) for the hierarchical normal model
//   y_i ~ normal(x_i, 1),  x_i = mu + sigma * raw_i,  i = 1..N
// in two forms that are equal in exact arithmetic:
//   cov:   x declared <offset=mu, multiplier=sigma>, x ~ normal(mu, sigma)
//          (offset_multiplier_constrain adds N log(sigma) to lp)
//   ncp:   raw ~ std_normal(), x = mu + sigma * raw (no change of variables)
// Both use propto = true, as the sampler does.  The reference is the
// analytic gradient in long double.  Prints the gradient error per
// component in the unconstrained (sampler) parameterization.
#include <stan/math.hpp>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using stan::math::var;
using Vec = Eigen::Matrix<var, Eigen::Dynamic, 1>;
using VecD = Eigen::VectorXd;

// theta = (mu, log_sigma, raw_1..raw_N)
double grad_cov(const VecD& theta, const VecD& y, VecD& g) {
  const int N = y.size();
  Vec th = theta.cast<var>();
  var lp = 0;
  var mu = th(0);
  var sigma = stan::math::lb_constrain(th(1), 0.0, lp);
  Vec raw = th.tail(N);
  Vec x = stan::math::offset_multiplier_constrain(raw, mu, sigma, lp);
  var target = lp;
  target += stan::math::normal_lpdf<true>(x, mu, sigma);
  target += stan::math::normal_lpdf<true>(y, x, 1.0);
  target.grad();
  g.resize(theta.size());
  for (int i = 0; i < theta.size(); ++i) g(i) = th(i).adj();
  double v = target.val();
  stan::math::recover_memory();
  return v;
}

double grad_ncp(const VecD& theta, const VecD& y, VecD& g) {
  const int N = y.size();
  Vec th = theta.cast<var>();
  var lp = 0;
  var mu = th(0);
  var sigma = stan::math::lb_constrain(th(1), 0.0, lp);
  Vec raw = th.tail(N);
  var target = lp;
  target += stan::math::std_normal_lpdf<true>(raw);
  Vec x = stan::math::add(mu, stan::math::multiply(sigma, raw));
  target += stan::math::normal_lpdf<true>(y, x, 1.0);
  target.grad();
  g.resize(theta.size());
  for (int i = 0; i < theta.size(); ++i) g(i) = th(i).adj();
  double v = target.val();
  stan::math::recover_memory();
  return v;
}

// analytic gradient, long double
void grad_ref(const VecD& theta, const VecD& y, std::vector<long double>& g) {
  const int N = y.size();
  long double mu = theta(0);
  long double s = theta(1);
  long double sigma = std::exp(s);
  g.assign(theta.size(), 0.0L);
  g[1] = 1.0L;  // log Jacobian of sigma = exp(s)
  for (int i = 0; i < N; ++i) {
    long double r = theta(2 + i);
    long double x = mu + sigma * r;
    long double dLdx = (long double)y(i) - x;
    g[0] += dLdx;
    g[1] += sigma * r * dLdx;
    g[2 + i] = -r + sigma * dLdx;
  }
}

// usage: grad_check [mu]   (default: mu from uniform(-2, 2))
int main(int argc, char** argv) {
  const int N = 1000;
  std::mt19937 rng(20261031);
  std::normal_distribution<double> nd(0.0, 1.0);
  std::uniform_real_distribution<double> ud(-2.0, 2.0);
  VecD y(N);
  for (int i = 0; i < N; ++i) y(i) = 0.5 * nd(rng) + nd(rng);
  VecD theta(N + 2);
  theta(0) = ud(rng);
  for (int i = 0; i < N; ++i) theta(2 + i) = ud(rng);
  if (argc > 1) theta(0) = std::atof(argv[1]);
  std::printf("mu = %g\n", theta(0));

  std::printf("%8s | %-29s | %-29s | %-29s\n", "sigma", "d/dmu: ref  cov  ncp",
              "d/dlog_sigma: ref  cov  ncp", "raw: max abs err cov  ncp");
  for (int e = -1; e >= -30; --e) {
    theta(1) = e * std::log(10.0);
    VecD gc, gn;
    std::vector<long double> gr;
    grad_cov(theta, y, gc);
    grad_ncp(theta, y, gn);
    grad_ref(theta, y, gr);
    double raw_c = 0, raw_n = 0;
    for (int i = 0; i < N; ++i) {
      raw_c = std::max(raw_c, (double)std::fabs(gc(2 + i) - gr[2 + i]));
      raw_n = std::max(raw_n, (double)std::fabs(gn(2 + i) - gr[2 + i]));
    }
    std::printf("%8.0e | %9.3g %9.3g %9.3g | %9.4g %9.4g %9.4g | %9.2g %9.2g\n",
                std::exp(theta(1)), (double)gr[0], gc(0), gn(0), (double)gr[1],
                gc(1), gn(1), raw_c, raw_n);
  }
  return 0;
}
