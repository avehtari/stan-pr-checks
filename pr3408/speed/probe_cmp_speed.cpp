// Bundle: root speed, develop against the PR (Triton job 20356756).
// Supports the speed table of #3408 Part 1.  Not run again here.
// probe_cmp_speed.cpp --- axis (c): speed, develop against branch.
//
// Compiled twice from this one source, once against each tree. Both
// binaries run in the same exclusive job, so they share the machine.
//
// Part 1, function level: the two roots in double, var, fvar<double> and
//         fvar<fvar<double>>.
// Part 2, distribution level, what a model pays: N = 1000 observations,
//         one shared var shape, value plus reverse sweep timed together,
//         reported as ns per observation.
//
// Timing rules from STAN_MATH_DISTRIBUTIONS.md section 4: an accumulator
// dependent `bump` is added to every input so the loop is not
// loop-invariant; a baseline is subtracted; every result is also given in
// cycles at the nominal clock so a removed call is visible.
//
// The input grid stays inside the region where develop terminates and is
// finite: a in [0.5, 60], z in [0.1, 40]. Develop hangs for z above about
// 700 in the fvar configurations and is NaN above a = 171.6.

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
  std::vector<double> a, z;
};

Inputs make_inputs() {
  Inputs in;
  for (int i = 0; i < BLOCK; ++i) {
    const double u = static_cast<double>(i) / (BLOCK - 1);
    const double v = static_cast<double>((i * 7919) % BLOCK) / (BLOCK - 1);
    in.a.push_back(std::exp(std::log(0.5) + u * (std::log(60.0) - std::log(0.5))));
    in.z.push_back(std::exp(std::log(0.1) + v * (std::log(40.0) - std::log(0.1))));
  }
  return in;
}

double nominal_ghz() {
  // "model name : Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz"
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
    const double v = kernel(in.a[i], in.z[i]);
    if (std::isfinite(v)) {
      acc += v;
    }
  }
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    const double bump = std::isfinite(acc) ? acc * 1e-300 : 0.0;
    for (int i = 0; i < BLOCK; ++i) {
      const double v = kernel(in.a[i] + bump, in.z[i] + bump);
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
double upper_d(double a, double z) {
  return stan::math::grad_reg_inc_gamma(a, z, stan::math::tgamma(a),
                                        stan::math::digamma(a));
}
double upper_v(double a, double z) {
  stan::math::nested_rev_autodiff nested;
  var av = a, zv = z;
  var r = stan::math::grad_reg_inc_gamma(av, zv, stan::math::tgamma(av),
                                         stan::math::digamma(av));
  r.grad();
  return r.val() + av.adj() * 0.0;
}
double upper_fd(double a, double z) {
  fvar<double> av = a, zv = z;
  av.d_ = 1.0;
  return stan::math::grad_reg_inc_gamma(av, zv, stan::math::tgamma(av),
                                        stan::math::digamma(av))
      .val_;
}
double upper_ffd(double a, double z) {
  fvar<fvar<double>> av = a, zv = z;
  av.d_ = 1.0;
  return stan::math::grad_reg_inc_gamma(av, zv, stan::math::tgamma(av),
                                        stan::math::digamma(av))
      .val_.val_;
}
double lower_d(double a, double z) {
  return stan::math::grad_reg_lower_inc_gamma(a, z);
}
double lower_v(double a, double z) {
  stan::math::nested_rev_autodiff nested;
  var av = a, zv = z;
  var r = stan::math::grad_reg_lower_inc_gamma(av, zv);
  r.grad();
  return r.val() + av.adj() * 0.0;
}
double lower_fd(double a, double z) {
  fvar<double> av = a, zv = z;
  av.d_ = 1.0;
  return stan::math::grad_reg_lower_inc_gamma(av, zv).val_;
}
double lower_ffd(double a, double z) {
  fvar<fvar<double>> av = a, zv = z;
  av.d_ = 1.0;
  return stan::math::grad_reg_lower_inc_gamma(av, zv).val_.val_;
}

// ── distribution level ─────────────────────────────────────────────────
template <typename F>
void time_dist(const char* name, double alpha, F fn, int n_obs, int reps,
               double ghz) {
  std::vector<double> y(n_obs);
  for (int i = 0; i < n_obs; ++i) {
    const double u = static_cast<double>(i) / (n_obs - 1);
    y[i] = alpha * std::exp(std::log(0.05) + u * (std::log(3.0) - std::log(0.05)));
  }
  double acc = 0.0;
  // warm up once
  {
    stan::math::nested_rev_autodiff nested;
    var a = alpha;
    var lp = fn(y, a);
    lp.grad();
    acc += a.adj();
  }
  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < reps; ++r) {
    stan::math::nested_rev_autodiff nested;
    var a = alpha + acc * 1e-300;
    var lp = fn(y, a);
    lp.grad();
    if (std::isfinite(a.adj())) {
      acc += a.adj();
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double ns_per_obs
      = std::chrono::duration<double, std::nano>(t1 - t0).count()
        / (static_cast<double>(reps) * n_obs);
  std::printf("%-30s alpha=%-5g %10.1f ns/obs %10.0f cycles/obs   (grad %.6g)\n",
              name, alpha, ns_per_obs, ns_per_obs * ghz, acc / (reps + 1));
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const double ghz = nominal_ghz();
  std::printf("# probe_cmp_speed; nominal clock %.2f GHz; BLOCK=%d\n", ghz, BLOCK);

  const Inputs in = make_inputs();

  std::printf("\n## function level, ns per call\n");
  std::printf("%-44s %10s %10s %10s\n", "kernel", "ns/call", "net ns",
              "net cycles");
  const double base = time_kernel(
      "baseline (a + z)", in, [](double a, double z) { return a + z; }, 300,
      0.0, ghz);
  time_kernel("tgamma + digamma only", in,
              [](double a, double z) {
                return stan::math::tgamma(a) + stan::math::digamma(a) + 0.0 * z;
              },
              300, base, ghz);
  time_kernel("grad_reg_inc_gamma double (+tgamma,digamma)", in, upper_d, 300,
              base, ghz);
  time_kernel("grad_reg_inc_gamma var, value+grad", in, upper_v, 50, base,
              ghz);
  time_kernel("grad_reg_inc_gamma fvar<double>", in, upper_fd, 100, base, ghz);
  time_kernel("grad_reg_inc_gamma fvar<fvar<double>>", in, upper_ffd, 50, base,
              ghz);
  time_kernel("grad_reg_lower_inc_gamma double", in, lower_d, 300, base, ghz);
  time_kernel("grad_reg_lower_inc_gamma var, value+grad", in, lower_v, 50,
              base, ghz);
  time_kernel("grad_reg_lower_inc_gamma fvar<double>", in, lower_fd, 100, base,
              ghz);
  time_kernel("grad_reg_lower_inc_gamma fvar<fvar<double>>", in, lower_ffd, 50,
              base, ghz);

  std::printf("\n## distribution level, N = 1000 observations, one shared var "
              "shape, value + reverse sweep\n");
  const int N = 1000, R = 40;
  for (double alpha : {5.0, 50.0}) {
    time_dist("gamma_cdf", alpha,
              [](const std::vector<double>& y, const var& a) {
                return stan::math::gamma_cdf(y, a, 1.0);
              },
              N, R, ghz);
    time_dist("gamma_lcdf", alpha,
              [](const std::vector<double>& y, const var& a) {
                return stan::math::gamma_lcdf(y, a, 1.0);
              },
              N, R, ghz);
    time_dist("gamma_lccdf", alpha,
              [](const std::vector<double>& y, const var& a) {
                return stan::math::gamma_lccdf(y, a, 1.0);
              },
              N, R, ghz);
    time_dist("chi_square_lccdf (nu = 2 alpha)", alpha,
              [](const std::vector<double>& y, const var& a) {
                return stan::math::chi_square_lccdf(y, 2.0 * a);
              },
              N, R, ghz);
    time_dist("inv_gamma_lcdf", alpha,
              [](const std::vector<double>& y, const var& a) {
                return stan::math::inv_gamma_lcdf(y, a, 1.0);
              },
              N, R, ghz);
  }
  std::printf("\n# done\n");
  return 0;
}
