// probe_cmp_beta_speed.cpp --- speed and Level 1/2 spot checks for the
// incomplete beta gradient roots, develop against stable-inc-beta.
//
// Compiled twice from this one source, once against each tree. Both
// binaries run in the same exclusive job, so they share the machine.
//
// Part 1, function level: grad_reg_inc_beta in double, var, fvar<double>
//         and fvar<fvar<double>>, and inc_beta_dda alone in double.
// Part 2, distribution level, what a model pays: N = 1000 observations,
//         shared var shapes, value plus reverse sweep, ns per observation.
// Part 3, Level 1 and 2 spot checks: the gradients of inc_beta, beta_cdf
//         and beta_lcdf at the defect points of the sweep, printed in hex
//         for scoring against the mpmath reference.
//
// Timing rules from STAN_MATH_DISTRIBUTIONS.md section 4: an accumulator
// dependent `bump` is added to every input so the loop is not
// loop-invariant; a baseline is subtracted; every result is also given in
// cycles at the nominal clock so a removed call is visible.
//
// The input grid stays where develop is finite and terminates quickly:
// a, b in [0.5, 60], z in [0.05, 0.95].

#include <stan/math/mix.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

using stan::math::fvar;
using stan::math::var;

constexpr int BLOCK = 1024;

struct Inputs {
  std::vector<double> a, b, z;
};

Inputs make_inputs() {
  Inputs in;
  for (int i = 0; i < BLOCK; ++i) {
    const double u = static_cast<double>(i) / (BLOCK - 1);
    const double v = static_cast<double>((i * 7919) % BLOCK) / (BLOCK - 1);
    const double w = static_cast<double>((i * 104729) % BLOCK) / (BLOCK - 1);
    in.a.push_back(std::exp(std::log(0.5) + u * std::log(120.0)));
    in.b.push_back(std::exp(std::log(0.5) + v * std::log(120.0)));
    in.z.push_back(0.05 + 0.90 * w);
  }
  return in;
}

double nominal_ghz() {
  if (FILE* f = std::fopen("/proc/cpuinfo", "r")) {
    char line[512];
    while (std::fgets(line, sizeof(line), f)) {
      const char* at = std::strstr(line, "@ ");
      if (std::strncmp(line, "model name", 10) == 0 && at != nullptr) {
        std::fclose(f);
        return std::atof(at + 2);
      }
    }
    std::fclose(f);
  }
  return 0.0;
}

template <typename F>
double time_kernel(const char* name, const Inputs& in, F kernel, int reps,
                   double baseline_ns, double ghz) {
  double acc = 0.0;
  for (int i = 0; i < BLOCK; ++i) {
    double v;
    try {
      v = kernel(in.a[i], in.b[i], in.z[i]);
    } catch (const std::exception&) {
      v = 0.0;
    }
    if (std::isfinite(v)) {
      acc += v;
    }
  }
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = std::isfinite(acc) ? acc * 1e-300 : 0.0;
    for (int i = 0; i < BLOCK; ++i) {
      double v;
      try {
        v = kernel(in.a[i] + bump, in.b[i] + bump, in.z[i] + bump);
      } catch (const std::exception&) {
        v = 0.0;
      }
      if (std::isfinite(v)) {
        acc += v;
      }
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count()
                    / (static_cast<double>(reps) * BLOCK);
  const double net = ns - baseline_ns;
  std::printf("%-44s %10.1f %10.1f %10.0f   (acc %.3g)\n", name, ns, net,
              net * ghz, acc);
  return ns;
}

// ── function-level kernels ─────────────────────────────────────────────
template <typename T>
T pair_sum(const T& a, const T& b, const T& z) {
  using stan::math::digamma;
  T g1 = 0, g2 = 0;
  stan::math::grad_reg_inc_beta(g1, g2, a, b, z, digamma(a), digamma(b),
                                digamma(a + b), stan::math::beta(a, b));
  return g1 + g2;
}
double pair_d(double a, double b, double z) { return pair_sum(a, b, z); }
double pair_v(double a, double b, double z) {
  stan::math::nested_rev_autodiff nested;
  var av = a, bv = b, zv = z;
  var r = pair_sum(av, bv, zv);
  r.grad();
  return r.val() + av.adj() * 0.0;
}
double pair_fd(double a, double b, double z) {
  fvar<double> av = a, bv = b, zv = z;
  av.d_ = 1.0;
  return pair_sum(av, bv, zv).val_;
}
double pair_ffd(double a, double b, double z) {
  fvar<fvar<double>> av = a, bv = b, zv = z;
  av.d_ = 1.0;
  return pair_sum(av, bv, zv).val_.val_;
}
double dda_d(double a, double b, double z) {
  return stan::math::inc_beta_dda(a, b, z, stan::math::digamma(a),
                                  stan::math::digamma(a + b));
}

// ── distribution level ─────────────────────────────────────────────────
template <typename Y, typename F>
void time_dist(const char* name, double alpha, double beta,
               const std::vector<Y>& y, F fn, int reps, double ghz) {
  const int n_obs = static_cast<int>(y.size());
  double acc = 0.0;
  {
    stan::math::nested_rev_autodiff nested;
    var a = alpha, b = beta;
    var lp = fn(y, a, b);
    lp.grad();
    if (std::isfinite(a.adj())) {
      acc += a.adj();
    }
  }
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    stan::math::nested_rev_autodiff nested;
    // a non-finite gradient (lccdf at an underflowed 1 - I) must not
    // poison the shape; the guard on acc keeps the bump finite
    const double bump = std::isfinite(acc) ? acc * 1e-300 : 0.0;
    var a = alpha + bump, b = beta + bump;
    var lp = fn(y, a, b);
    lp.grad();
    if (std::isfinite(a.adj())) {
      acc += a.adj();
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double ns_per_obs
      = std::chrono::duration<double, std::nano>(t1 - t0).count()
        / (static_cast<double>(reps) * n_obs);
  std::printf("%-26s a=%-4g b=%-4g %9.1f ns/obs %9.0f cycles/obs   (grad %.6g)\n",
              name, alpha, beta, ns_per_obs, ns_per_obs * ghz,
              acc / (reps + 1));
}

// ── Level 1 and 2 spot checks ──────────────────────────────────────────
template <typename F>
void spot(const char* name, double a, double b, double y, F fn) {
  double val = std::nan(""), ga = std::nan(""), gb = std::nan("");
  try {
    stan::math::nested_rev_autodiff nested;
    var av = a, bv = b;
    var r = fn(y, av, bv);
    r.grad();
    val = r.val();
    ga = av.adj();
    gb = bv.adj();
  } catch (const std::exception& e) {
    std::printf("SPOT %-12s %a %a %a threw: %s\n", name, a, b, y, e.what());
    return;
  }
  std::printf("SPOT %-12s %a %a %a %a %a %a   # (%g, %g, %g) val %.6g da %.6g db %.6g\n",
              name, a, b, y, val, ga, gb, a, b, y, val, ga, gb);
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const double ghz = nominal_ghz();
  std::printf("# probe_cmp_beta_speed; nominal clock %.2f GHz; BLOCK=%d\n", ghz,
              BLOCK);

  const Inputs in = make_inputs();

  std::printf("\n## function level, ns per call\n");
  std::printf("%-44s %10s %10s %10s\n", "kernel", "ns/call", "net ns",
              "net cycles");
  const double base = time_kernel(
      "baseline (a + b + z)", in,
      [](double a, double b, double z) { return a + b + z; }, 300, 0.0, ghz);
  time_kernel("3 digamma + beta only", in,
              [](double a, double b, double z) {
                return stan::math::digamma(a) + stan::math::digamma(b)
                       + stan::math::digamma(a + b) + stan::math::beta(a, b)
                       + 0.0 * z;
              },
              300, base, ghz);
  time_kernel("grad_reg_inc_beta double (+3 digamma, beta)", in, pair_d, 100,
              base, ghz);
  time_kernel("grad_reg_inc_beta var, value+grad", in, pair_v, 20, base, ghz);
  time_kernel("grad_reg_inc_beta fvar<double>", in, pair_fd, 40, base, ghz);
  time_kernel("grad_reg_inc_beta fvar<fvar<double>>", in, pair_ffd, 20, base,
              ghz);
  time_kernel("inc_beta_dda double (+2 digamma)", in, dda_d, 100, base, ghz);

  std::printf("\n## distribution level, N = 1000 observations, shared var "
              "shapes, value + reverse sweep\n");
  const int N = 1000, R = 20;
  for (double s : {5.0, 50.0}) {
    // keep 1 - I representable for lccdf: y within +-4 sd of the mean
    const double half = std::min(0.48, 4.0 * std::sqrt(0.25 / (2.0 * s + 1.0)));
    std::vector<double> yb(N);
    for (int i = 0; i < N; ++i) {
      yb[i] = 0.5 - half + 2.0 * half * static_cast<double>(i) / (N - 1);
    }
    time_dist("beta_cdf", s, s, yb,
              [](const std::vector<double>& y, const var& a, const var& b) {
                return stan::math::beta_cdf(y, a, b);
              },
              R, ghz);
    time_dist("beta_lcdf", s, s, yb,
              [](const std::vector<double>& y, const var& a, const var& b) {
                return stan::math::beta_lcdf(y, a, b);
              },
              R, ghz);
    time_dist("beta_lccdf", s, s, yb,
              [](const std::vector<double>& y, const var& a, const var& b) {
                return stan::math::beta_lccdf(y, a, b);
              },
              R, ghz);
  }
  std::vector<int> n_obs(N);
  for (int i = 0; i < N; ++i) {
    n_obs[i] = i % 60;
  }
  time_dist("neg_binomial_cdf", 5.0, 0.3, n_obs,
            [](const std::vector<int>& n, const var& a, const var& b) {
              return stan::math::neg_binomial_cdf(n, a, b);
            },
            R, ghz);
  time_dist("neg_binomial_lcdf", 5.0, 0.3, n_obs,
            [](const std::vector<int>& n, const var& a, const var& b) {
              return stan::math::neg_binomial_lcdf(n, a, b);
            },
            R, ghz);
  std::vector<double> yt(N);
  for (int i = 0; i < N; ++i) {
    yt[i] = -5.0 + 10.0 * static_cast<double>(i) / (N - 1);
  }
  time_dist("student_t_lcdf (nu = a)", 3.0, 1.0, yt,
            [](const std::vector<double>& y, const var& a, const var& b) {
              return stan::math::student_t_lcdf(y, a, 0.0, b);
            },
            R, ghz);

  std::printf("\n## Level 1 and 2 spot checks: name a b y value d/da d/db\n");
  auto f_inc = [](double y, const var& a, const var& b) {
    return stan::math::inc_beta(a, b, y);
  };
  auto f_cdf = [](double y, const var& a, const var& b) {
    return stan::math::beta_cdf(y, a, b);
  };
  auto f_lcdf = [](double y, const var& a, const var& b) {
    return stan::math::beta_lcdf(y, a, b);
  };
  auto f_lccdf = [](double y, const var& a, const var& b) {
    return stan::math::beta_lccdf(y, a, b);
  };
  const double P[][3] = {{600, 600, 0.5},     {352, 590, 0.758},
                         {1.145, 6786, 1.8e-4}, {446, 5, 0.76},
                         {10, 1, 0.5},          {2, 3, 0.25},
                         {2500, 2, 0.999},      {1, 100, 0.01}};
  for (auto& p : P) {
    spot("inc_beta", p[0], p[1], p[2], f_inc);
    spot("beta_cdf", p[0], p[1], p[2], f_cdf);
    spot("beta_lcdf", p[0], p[1], p[2], f_lcdf);
    spot("beta_lccdf", p[0], p[1], p[2], f_lccdf);
  }
  std::printf("\n# done\n");
  return 0;
}
