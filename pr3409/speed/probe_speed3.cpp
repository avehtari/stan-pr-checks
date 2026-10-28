// probe_speed3.cpp --- speed of the incomplete beta gradient roots and
// their consumers on three trees: develop, the stable-inc-beta branch, and
// the branch with the reverse-mode / one-sided refactor (stan-dev/math#3409
// commit d2eeacd676).
//
// Compiled once per tree from this one source. Build with -DHAVE_GRAD_RIB
// on the tree that has grad_reg_inc_beta<GRAD_RIB>; the other trees use the
// equivalent older calls.
//
// Part 1, function level, ns per call.
// Part 2, distribution level: N = 1000 observations, shared shapes, value
//         plus reverse sweep, ns per observation.
// Part 3, distribution level with fvar<var> shapes: value, tangent, and the
//         reverse sweep of the tangent (one Hessian-vector product).
// Part 4, spot checks at the defect points (SPOT lines, hex), for
//         score_spot_beta.py.
//
// Timing rules from STAN_MATH_DISTRIBUTIONS.md section 4: an accumulator
// dependent bump on every input, a baseline subtracted, cycles printed.
// Grid: a, b in [0.5, 60], z in [0.05, 0.95], where develop terminates.

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
  std::printf("K %-46s %10.1f %10.1f %10.0f   (acc %.3g)\n", name, ns, net,
              net * ghz, acc);
  return ns;
}

// ── function-level kernels ──
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
  return r.val() + 0.0 * av.adj();
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
// d/da alone, in double, including the digamma calls it needs
double one_sided_d(double a, double b, double z) {
  using stan::math::digamma;
#ifdef HAVE_GRAD_RIB
  return stan::math::grad_reg_inc_beta<stan::math::GRAD_RIB::G1>(
      a, b, z, digamma(a), 0.0, digamma(a + b), 0.0);
#else
  return stan::math::inc_beta_dda(a, b, z, digamma(a), digamma(a + b));
#endif
}
double inc_beta_vvd(double a, double b, double z) {
  stan::math::nested_rev_autodiff nested;
  var av = a, bv = b;
  var r = stan::math::inc_beta(av, bv, z);
  r.grad();
  return r.val() + av.adj() + bv.adj();
}
double inc_beta_vdd(double a, double b, double z) {
  stan::math::nested_rev_autodiff nested;
  var av = a;
  var r = stan::math::inc_beta(av, b, z);
  r.grad();
  return r.val() + av.adj();
}
double inc_beta_fv(double a, double b, double z) {
  // one Hessian-vector product of inc_beta in (a, b)
  stan::math::nested_rev_autodiff nested;
  var av = a, bv = b;
  fvar<var> af(av, 1.0), bf(bv, 1.0);
  fvar<var> zf = z;
  fvar<var> r = stan::math::inc_beta(af, bf, zf);
  r.d_.grad();
  return r.val_.val() + av.adj() + bv.adj();
}

// ── distribution level, var ──
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
  std::printf("D %-30s a=%-5g b=%-5g %9.1f ns/obs %9.0f cyc/obs (grad %.6g)\n",
              name, alpha, beta, ns_per_obs, ns_per_obs * ghz,
              acc / (reps + 1));
}

// ── distribution level, fvar<var>: value, tangent, reverse of the tangent ──
template <typename Y, typename F>
void time_dist_fv(const char* name, double alpha, double beta,
                  const std::vector<Y>& y, F fn, int reps, double ghz) {
  const int n_obs = static_cast<int>(y.size());
  double acc = 0.0;
  auto once = [&](double bump) {
    stan::math::nested_rev_autodiff nested;
    var av = alpha + bump, bv = beta + bump;
    fvar<var> a(av, 1.0), b(bv, 1.0);
    fvar<var> lp = fn(y, a, b);
    lp.d_.grad();
    if (std::isfinite(av.adj())) {
      acc += av.adj();
    }
  };
  once(0.0);
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    once(std::isfinite(acc) ? acc * 1e-300 : 0.0);
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double ns_per_obs
      = std::chrono::duration<double, std::nano>(t1 - t0).count()
        / (static_cast<double>(reps) * n_obs);
  std::printf("F %-30s a=%-5g b=%-5g %9.1f ns/obs %9.0f cyc/obs (hvp %.6g)\n",
              name, alpha, beta, ns_per_obs, ns_per_obs * ghz,
              acc / (reps + 1));
}

// ── spot checks ──
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
  std::printf(
      "SPOT %-12s %a %a %a %a %a %a   # (%g, %g, %g) val %.6g da %.6g db "
      "%.6g\n",
      name, a, b, y, val, ga, gb, a, b, y, val, ga, gb);
}

std::vector<double> beta_y(double s, int n) {
  // keep 1 - I representable for lccdf: y within +-4 sd of the mean
  const double half = std::min(0.48, 4.0 * std::sqrt(0.25 / (2.0 * s + 1.0)));
  std::vector<double> y(n);
  for (int i = 0; i < n; ++i) {
    y[i] = 0.5 - half + 2.0 * half * static_cast<double>(i) / (n - 1);
  }
  return y;
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const double ghz = nominal_ghz();
  std::printf("# probe_speed3; nominal clock %.2f GHz; BLOCK=%d\n", ghz, BLOCK);
  const Inputs in = make_inputs();

  std::printf("\n## function level: name ns/call net-ns net-cycles\n");
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
  time_kernel("grad_reg_inc_beta both, double", in, pair_d, 100, base, ghz);
  time_kernel("d/da only, double (+2 digamma)", in, one_sided_d, 100, base,
              ghz);
  time_kernel("grad_reg_inc_beta both, var, val+grad", in, pair_v, 20, base,
              ghz);
  time_kernel("grad_reg_inc_beta both, fvar<double>", in, pair_fd, 40, base,
              ghz);
  time_kernel("grad_reg_inc_beta both, fvar<fvar<double>>", in, pair_ffd, 20,
              base, ghz);
  time_kernel("inc_beta(var, var, double), val+grad", in, inc_beta_vvd, 40,
              base, ghz);
  time_kernel("inc_beta(var, double, double), val+grad", in, inc_beta_vdd, 40,
              base, ghz);
  time_kernel("inc_beta fvar<var> (a, b), one HVP", in, inc_beta_fv, 20, base,
              ghz);

  std::printf("\n## distribution level, var, N = 1000, shared shapes\n");
  const int N = 1000, R = 20;
  using VD = std::vector<double>;
  for (double s : {5.0, 50.0}) {
    const VD yb = beta_y(s, N);
    time_dist("beta_cdf", s, s, yb,
              [](const VD& y, const var& a, const var& b) {
                return stan::math::beta_cdf(y, a, b);
              },
              R, ghz);
    time_dist("beta_lcdf", s, s, yb,
              [](const VD& y, const var& a, const var& b) {
                return stan::math::beta_lcdf(y, a, b);
              },
              R, ghz);
    time_dist("beta_lccdf", s, s, yb,
              [](const VD& y, const var& a, const var& b) {
                return stan::math::beta_lccdf(y, a, b);
              },
              R, ghz);
    time_dist("beta_lcdf, alpha only", s, s, yb,
              [s](const VD& y, const var& a, const var&) {
                return stan::math::beta_lcdf(y, a, s);
              },
              R, ghz);
    time_dist("beta_proportion_lcdf", 0.5, 2 * s, yb,
              [](const VD& y, const var& a, const var& b) {
                return stan::math::beta_proportion_lcdf(y, a, b);
              },
              R, ghz);
  }
  std::vector<int> n_obs(N);
  for (int i = 0; i < N; ++i) {
    n_obs[i] = i % 60;
  }
  using VI = std::vector<int>;
  time_dist("neg_binomial_cdf", 5.0, 0.3, n_obs,
            [](const VI& n, const var& a, const var& b) {
              return stan::math::neg_binomial_cdf(n, a, b);
            },
            R, ghz);
  time_dist("neg_binomial_lcdf", 5.0, 0.3, n_obs,
            [](const VI& n, const var& a, const var& b) {
              return stan::math::neg_binomial_lcdf(n, a, b);
            },
            R, ghz);
  time_dist("neg_binomial_lccdf", 5.0, 0.3, n_obs,
            [](const VI& n, const var& a, const var& b) {
              return stan::math::neg_binomial_lccdf(n, a, b);
            },
            R, ghz);
  time_dist("neg_binomial_2_cdf", 16.0, 5.0, n_obs,
            [](const VI& n, const var& a, const var& b) {
              return stan::math::neg_binomial_2_cdf(n, a, b);
            },
            R, ghz);
  VD yt(N);
  for (int i = 0; i < N; ++i) {
    yt[i] = -5.0 + 10.0 * static_cast<double>(i) / (N - 1);
  }
  time_dist("student_t_lcdf (nu = a)", 3.0, 1.0, yt,
            [](const VD& y, const var& a, const var& b) {
              return stan::math::student_t_lcdf(y, a, 0.0, b);
            },
            R, ghz);
  time_dist("student_t_cdf (nu = a)", 3.0, 1.0, yt,
            [](const VD& y, const var& a, const var& b) {
              return stan::math::student_t_cdf(y, a, 0.0, b);
            },
            R, ghz);

  std::printf("\n## distribution level, fvar<var>, N = 1000, one HVP\n");
  using FV = fvar<var>;
  for (double s : {5.0, 50.0}) {
    const VD yb = beta_y(s, N);
    time_dist_fv("beta_lcdf", s, s, yb,
                 [](const VD& y, const FV& a, const FV& b) {
                   return stan::math::beta_lcdf(y, a, b);
                 },
                 5, ghz);
    time_dist_fv("beta_cdf", s, s, yb,
                 [](const VD& y, const FV& a, const FV& b) {
                   return stan::math::beta_cdf(y, a, b);
                 },
                 5, ghz);
  }
  time_dist_fv("neg_binomial_lcdf", 5.0, 0.3, n_obs,
               [](const VI& n, const FV& a, const FV& b) {
                 return stan::math::neg_binomial_lcdf(n, a, b);
               },
               5, ghz);

  std::printf("\n## spot checks: name a b y value d/da d/db\n");
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
  const double P[][3] = {{600, 600, 0.5},       {352, 590, 0.758},
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
