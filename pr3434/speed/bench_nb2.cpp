// Time neg_binomial_2_lpmf with std::vector<int> n and vector mu, phi of
// size N = 1000: double; mu var; phi var; both var (value and gradient).
// Prints ns per element and, for the bitwise check, the value and the
// first adjoints in hex.  Usage: bench_nb2 <label> <rounds>
#include <stan/math/rev.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using stan::math::var;
constexpr int N = 1000;
constexpr int REPS = 2000;

struct data {
  std::vector<int> n;
  std::vector<double> mu, phi;
};

data make_data() {
  std::mt19937_64 rng(20261020);
  std::uniform_real_distribution<double> lmu(-2, 4), lphi(-2, 5);
  data d;
  for (int i = 0; i < N; ++i) {
    d.mu.push_back(std::exp(lmu(rng)));
    d.phi.push_back(std::exp(lphi(rng)));
    std::poisson_distribution<int> p(d.mu.back());
    d.n.push_back(p(rng));
  }
  return d;
}

// Kind 0: double, 1: mu var, 2: phi var, 3: both var.
template <int Kind>
double once(const data& d, double bump, double* adj) {
  if constexpr (Kind == 0) {
    std::vector<double> mu(d.mu);
    mu[0] += bump;
    return stan::math::neg_binomial_2_lpmf(d.n, mu, d.phi);
  } else {
    std::vector<var> mu_v(d.mu.begin(), d.mu.end());
    std::vector<var> phi_v(d.phi.begin(), d.phi.end());
    std::vector<double> mu_d(d.mu);
    std::vector<double> phi_d(d.phi);
    mu_d[0] += bump;
    var lp;
    if constexpr (Kind == 1) {
      lp = stan::math::neg_binomial_2_lpmf(d.n, mu_v, phi_d);
    } else if constexpr (Kind == 2) {
      lp = stan::math::neg_binomial_2_lpmf(d.n, mu_d, phi_v);
    } else {
      lp = stan::math::neg_binomial_2_lpmf(d.n, mu_v, phi_v);
    }
    lp.grad();
    const double res = lp.val();
    adj[0] = mu_v[0].adj() + mu_v[N - 1].adj();
    adj[1] = phi_v[0].adj() + phi_v[N - 1].adj();
    stan::math::recover_memory();
    return res;
  }
}

template <int Kind>
void run(const char* label, int rounds, const data& d) {
  static const char* names[] = {"double", "mu_var", "phi_var", "both_var"};
  double adj[2] = {0, 0};
  const double v0 = once<Kind>(d, 0.0, adj);
  std::printf("VALUE %s %s %a %a %a\n", label, names[Kind], v0, adj[0],
              adj[1]);
  for (int r = 0; r < rounds; ++r) {
    double acc = 0;
    const int reps = Kind == 0 ? REPS : REPS / 4;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < reps; ++i) {
      acc += once<Kind>(d, acc * 1e-300, adj);
    }
    const auto t1 = std::chrono::steady_clock::now();
    std::printf("RESULT %s %s %d %.3f sink %g\n", label, names[Kind], r + 1,
                std::chrono::duration<double, std::nano>(t1 - t0).count()
                    / (static_cast<double>(reps) * N),
                acc);
  }
}

int main(int argc, char** argv) {
  const char* label = argc > 1 ? argv[1] : "?";
  const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
  const data d = make_data();
  run<0>(label, rounds, d);
  run<1>(label, rounds, d);
  run<2>(label, rounds, d);
  run<3>(label, rounds, d);
  return 0;
}
