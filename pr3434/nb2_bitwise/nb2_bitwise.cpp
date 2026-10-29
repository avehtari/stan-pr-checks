// nb2_bitwise.cpp --- neg_binomial_2_lpmf value and every adjoint in hex,
// for a bitwise comparison of two trees (claim: "evaluates log(mu + phi)
// once ... The results are bitwise identical"; commit 1c8afb741a against
// its parent 0377416157).
//
// Data: 1000 instances, the seed of bench_nb2.cpp (mu from e^-2 to e^4, phi
// from e^-2 to e^5, Poisson counts), and 1000 more with phi from 1e-8 to
// 1e15 and counts 0 to 1000.  Calls: vector mu and phi with mu var, phi var,
// both var; scalar phi var with vector mu var; and every instance alone.
// Output: one line per printed number, `<call> <index> <hex>`.
//
// Usage: nb2_bitwise > out.txt; then compare two outputs with cmp or diff.

#include <stan/math/rev.hpp>

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using stan::math::var;

namespace {

struct Data {
  std::vector<int> n;
  std::vector<double> mu, phi;
};

Data make_data() {
  Data d;
  std::mt19937_64 rng(20261020);
  std::uniform_real_distribution<double> lmu(-2, 4), lphi(-2, 5);
  for (int i = 0; i < 1000; ++i) {
    d.mu.push_back(std::exp(lmu(rng)));
    d.phi.push_back(std::exp(lphi(rng)));
    std::poisson_distribution<int> p(d.mu.back());
    d.n.push_back(p(rng));
  }
  std::uniform_real_distribution<double> lphi2(std::log(1e-8),
                                               std::log(1e15));
  std::uniform_int_distribution<int> cnt(0, 1000);
  for (int i = 0; i < 1000; ++i) {
    d.mu.push_back(std::exp(lmu(rng)));
    d.phi.push_back(std::exp(lphi2(rng)));
    d.n.push_back(cnt(rng));
  }
  return d;
}

void print_vec(const char* call, const std::vector<var>& v) {
  for (size_t i = 0; i < v.size(); ++i) {
    std::printf("%s %zu %a\n", call, i, v[i].adj());
  }
}

}  // namespace

int main() {
  const Data d = make_data();
  const size_t N = d.n.size();
  {
    std::vector<var> mu(d.mu.begin(), d.mu.end());
    var lp = stan::math::neg_binomial_2_lpmf(d.n, mu, d.phi);
    lp.grad();
    std::printf("mu_var value %a\n", lp.val());
    print_vec("mu_var d_mu", mu);
    stan::math::recover_memory();
  }
  {
    std::vector<var> phi(d.phi.begin(), d.phi.end());
    var lp = stan::math::neg_binomial_2_lpmf(d.n, d.mu, phi);
    lp.grad();
    std::printf("phi_var value %a\n", lp.val());
    print_vec("phi_var d_phi", phi);
    stan::math::recover_memory();
  }
  {
    std::vector<var> mu(d.mu.begin(), d.mu.end());
    std::vector<var> phi(d.phi.begin(), d.phi.end());
    var lp = stan::math::neg_binomial_2_lpmf(d.n, mu, phi);
    lp.grad();
    std::printf("both_var value %a\n", lp.val());
    print_vec("both_var d_mu", mu);
    print_vec("both_var d_phi", phi);
    stan::math::recover_memory();
  }
  {
    std::vector<var> mu(d.mu.begin(), d.mu.end());
    var phi = 3.7;
    var lp = stan::math::neg_binomial_2_lpmf(d.n, mu, phi);
    lp.grad();
    std::printf("scalar_phi value %a\n", lp.val());
    print_vec("scalar_phi d_mu", mu);
    std::printf("scalar_phi d_phi %a\n", phi.adj());
    stan::math::recover_memory();
  }
  std::printf("double value %a\n",
              stan::math::neg_binomial_2_lpmf(d.n, d.mu, d.phi));
  for (size_t i = 0; i < N; ++i) {
    var mu = d.mu[i];
    var phi = d.phi[i];
    var lp = stan::math::neg_binomial_2_lpmf(d.n[i], mu, phi);
    lp.grad();
    std::printf("single %zu %a %a %a\n", i, lp.val(), mu.adj(), phi.adj());
    stan::math::recover_memory();
  }
  return 0;
}
