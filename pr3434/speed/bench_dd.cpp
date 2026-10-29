// bench_dd.cpp --- cost of digamma_diff at small shapes, alone and in the
// gradients that use it (stan-dev/math#3153, decision 12). Built four times
// by run_bench_dd.sh: develop (overlay of the 23 develop headers), branch,
// and the branch with candidate digamma_diff C1 or C2 (dd_candidates.hpp
// describes them). The develop binary still has the branch digamma_diff
// (its kernels dd_* show the noise between binaries); its distributions use
// the plain digamma differences.
//
// Method: STAN_MATH_DISTRIBUTIONS.md section 4. The loop carries
// bump = acc * 1e-300, which changes no result bit and stops the compiler
// from hoisting the call.
//
// Output: TIME <kernel> <ns per call>

#include <stan/math/rev.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using stan::math::var;

namespace {

constexpr int BLOCK = 256;

template <typename F>
void time_it(const char* name, int reps, F&& f) {
  double acc = 0.0;
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = acc * 1e-300;
    for (int i = 0; i < BLOCK; ++i) {
      acc += f(i, bump);
    }
  }
  auto t1 = std::chrono::steady_clock::now();
  const double ns
      = 1e9 * std::chrono::duration<double>(t1 - t0).count() / (reps * BLOCK);
  std::printf("TIME %s %.3f\n", name, ns);
  std::printf("# %s acc=%.17g\n", name, acc);
}

std::vector<double> grid(double lo, double hi, bool log_scale) {
  std::vector<double> v(BLOCK);
  for (int i = 0; i < BLOCK; ++i) {
    const double t = (i % 97) / 97.0;
    v[i] = log_scale ? lo * std::pow(hi / lo, t) : lo + (hi - lo) * t;
  }
  return v;
}

std::vector<int> counts(int lo, int hi) {
  std::vector<int> v(BLOCK);
  for (int i = 0; i < BLOCK; ++i) {
    v[i] = lo + (i * 7) % (hi - lo + 1);
  }
  return v;
}

template <typename F>
double grad1(F&& f, double a) {
  var va = a;
  var lp = f(va);
  lp.grad();
  const double out = lp.val() + va.adj();
  stan::math::recover_memory();
  return out;
}

template <typename F>
double grad2(F&& f, double a, double b) {
  var va = a;
  var vb = b;
  var lp = f(va, vb);
  lp.grad();
  const double out = lp.val() + va.adj() + vb.adj();
  stan::math::recover_memory();
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  const int reps = argc > 1 ? std::atoi(argv[1]) : 200;
  const int n = reps * 20;
  const auto x_small = grid(0.5, 10.0, false);
  const auto x_large = grid(10.0, 1e4, true);
  const auto y_large = grid(10.0, 1e3, true);
  const auto n_small = counts(0, 8);
  const auto n_mid = counts(9, 100);
  const auto mu = grid(0.5, 20.0, false);

#ifdef HAS_DIGAMMA_DIFF
  time_it("dd_small_int", n, [&](int i, double b) {
    return stan::math::digamma_diff(x_small[i] + b, n_small[i]);
  });
#endif
  time_it("pair_small_int", n, [&](int i, double b) {
    const double x = x_small[i] + b;
    return stan::math::digamma(x + n_small[i]) - stan::math::digamma(x);
  });
#ifdef HAS_DIGAMMA_DIFF
  time_it("dd_mid_int", n, [&](int i, double b) {
    return stan::math::digamma_diff(x_small[i] + b, n_mid[i]);
  });
#endif
  time_it("pair_mid_int", n, [&](int i, double b) {
    const double x = x_small[i] + b;
    return stan::math::digamma(x + n_mid[i]) - stan::math::digamma(x);
  });
#ifdef HAS_DIGAMMA_DIFF
  time_it("dd_half", n, [&](int i, double b) {
    return stan::math::digamma_diff(x_small[i] + b, 0.5);
  });
  time_it("dd_large_x_small_int", n, [&](int i, double b) {
    return stan::math::digamma_diff(x_large[i] + b, n_small[i]);
  });
#endif

  time_it("nb2_grad_phi_small_n_small", n, [&](int i, double b) {
    return grad1(
        [&](const var& phi) {
          return stan::math::neg_binomial_2_lpmf(n_small[i], mu[i], phi);
        },
        x_small[i] + b);
  });
  time_it("nb2_grad_phi_small_n_mid", n, [&](int i, double b) {
    return grad1(
        [&](const var& phi) {
          return stan::math::neg_binomial_2_lpmf(n_mid[i], mu[i], phi);
        },
        x_small[i] + b);
  });
  time_it("nb2_grad_phi_large_n_small", n, [&](int i, double b) {
    return grad1(
        [&](const var& phi) {
          return stan::math::neg_binomial_2_lpmf(n_small[i], mu[i], phi);
        },
        x_large[i] + b);
  });
  time_it("nb_grad_alpha_small", n, [&](int i, double b) {
    return grad2(
        [&](const var& alpha, const var& beta) {
          return stan::math::neg_binomial_lpmf(n_small[i], alpha, beta);
        },
        x_small[i] + b, 0.7);
  });
  time_it("nb2log_grad_phi_small", n, [&](int i, double b) {
    return grad1(
        [&](const var& phi) {
          return stan::math::neg_binomial_2_log_lpmf(n_small[i], 0.4, phi);
        },
        x_small[i] + b);
  });
  time_it("bb_grad_small", n, [&](int i, double b) {
    return grad2(
        [&](const var& a, const var& bb) {
          return stan::math::beta_binomial_lpmf(n_small[i] + 6, 20, a, bb);
        },
        x_small[i] + b, x_small[(i + 31) % BLOCK]);
  });
  time_it("lbeta_grad_small", n, [&](int i, double b) {
    return grad2(
        [&](const var& x, const var& y) { return stan::math::lbeta(x, y); },
        x_small[i] + b, x_small[(i + 31) % BLOCK]);
  });
  time_it("lbeta_grad_mixed", n, [&](int i, double b) {
    return grad2(
        [&](const var& x, const var& y) { return stan::math::lbeta(x, y); },
        x_small[i] + b, y_large[i]);
  });

  // one phi, 1000 counts from 0 to 20: the usual model
  std::vector<int> ns(1000);
  for (int j = 0; j < 1000; ++j) {
    ns[j] = (j * 13) % 21;
  }
  time_it("nb2_vec1000_grad_phi_small", reps, [&](int i, double b) {
    return grad1(
        [&](const var& phi) {
          return stan::math::neg_binomial_2_lpmf(ns, 4.0, phi);
        },
        x_small[i] + b);
  });
  return 0;
}
