// Branch normal_lpdf against the develop copy (normal_lpdf_dev) for the
// type cases of cases.hpp and the input cases below: the same value or the
// same exception text, gradients within TOL relative.  Reverse mode for all
// type cases, fvar<double> tangents for the cases without varmat.
// Also counts heap allocations of at least 8 N bytes during one gradient
// evaluation (N = 1000, regular input) with a malloc interposer.
// Output: CHECK, HEAP and CHECK_SUMMARY lines.

#include "cases.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

// ── malloc counting ──
extern "C" void* __libc_malloc(size_t);
static size_t g_threshold = static_cast<size_t>(-1);
static int g_big = 0;
extern "C" void* malloc(size_t n) {
  if (n >= g_threshold) {
    ++g_big;
  }
  return __libc_malloc(n);
}

using Eigen::VectorXd;
using stan::math::var;
using fd = stan::math::fvar<double>;

constexpr double TOL = 1e-12;
static int n_fail = 0, n_ok = 0;
volatile double sink = 0;

struct inputs {
  std::string name;
  VectorXd y, mu, sigma;
};

struct outcome {
  std::string what;
  double val = 0;
  std::vector<double> grad, cond;
};

template <bool Dev, bool P, typename Y, typename M, typename S>
auto call(const Y& y, const M& mu, const S& s) {
  if constexpr (Dev) {
    return stan::math::normal_lpdf_dev<P>(y, mu, s);
  } else {
    return stan::math::normal_lpdf<P>(y, mu, s);
  }
}

// Scale of the rounding error of each gradient entry: the sum of the
// absolute values of its terms (scalar arguments) or 0 (vector arguments).
static std::vector<double> conds(const inputs& in, int ky, int km, int ks) {
  const Eigen::Index n = std::max({in.y.size(), in.mu.size(), in.sigma.size()});
  auto at = [](const VectorXd& v, Eigen::Index i) {
    return v.size() == 1 ? v[0] : v[i];
  };
  double c_loc = 0, c_scale = 0;
  auto fits = [n](const VectorXd& v) { return v.size() == 1 || v.size() == n; };
  if (in.y.size() > 0 && in.mu.size() > 0 && in.sigma.size() > 0
      && fits(in.y) && fits(in.mu) && fits(in.sigma)) {
    for (Eigen::Index i = 0; i < n; ++i) {
      const double s = at(in.sigma, i);
      const double z = (at(in.y, i) - at(in.mu, i)) / s;
      c_loc += std::fabs(z / s);
      c_scale += (z * z + 1) / s;
    }
  }
  std::vector<double> c;
  auto add = [&](int k, const VectorXd& v, double cs) {
    if (!nc::is_ad_kind(k)) {
      return;
    }
    if (nc::is_scalar_kind(k)) {
      c.push_back(cs);
    } else {
      for (Eigen::Index i = 0; i < v.size(); ++i) {
        c.push_back(0);
      }
    }
  };
  add(ky, in.y, c_loc);
  add(km, in.mu, c_loc);
  add(ks, in.sigma, c_scale);
  return c;
}

template <bool Dev, bool P, int KY, int KM, int KS>
outcome run_rev(const inputs& in) {
  outcome o;
  try {
    auto y = nc::make<KY, var>(in.y);
    auto mu = nc::make<KM, var>(in.mu);
    auto s = nc::make<KS, var>(in.sigma);
    auto lp = call<Dev, P>(y, mu, s);
    if constexpr (std::is_same_v<decltype(lp), var>) {
      lp.grad();
      o.val = lp.val();
    } else {
      o.val = lp;
    }
    for (const VectorXd& a :
         {nc::adjoints(y), nc::adjoints(mu), nc::adjoints(s)}) {
      o.grad.insert(o.grad.end(), a.data(), a.data() + a.size());
    }
    o.what = "value";
  } catch (const std::exception& e) {
    o.what = e.what();
  }
  stan::math::recover_memory();
  return o;
}

// fvar<double>: value and the tangent with d_ = 1 on argument dir.
template <bool Dev, bool P, int KY, int KM, int KS>
outcome run_fwd(const inputs& in) {
  outcome o;
  o.what = "value";
  for (int dir = 0; dir < 3; ++dir) {
    const int k = dir == 0 ? KY : (dir == 1 ? KM : KS);
    if (!nc::is_ad_kind(k)) {
      continue;
    }
    try {
      auto y = nc::make<KY, fd>(in.y, dir == 0);
      auto mu = nc::make<KM, fd>(in.mu, dir == 1);
      auto s = nc::make<KS, fd>(in.sigma, dir == 2);
      auto lp = call<Dev, P>(y, mu, s);
      if constexpr (std::is_same_v<decltype(lp), fd>) {
        o.val = lp.val_;
        o.grad.push_back(lp.d_);
      } else {
        o.val = lp;
      }
    } catch (const std::exception& e) {
      o.what = e.what();
      return o;
    }
  }
  return o;
}

static bool close(double a, double b, double scale, double& rel) {
  if (std::isnan(a) || std::isnan(b)) {
    return std::isnan(a) && std::isnan(b);
  }
  if (std::isinf(a) || std::isinf(b)) {
    return a == b;
  }
  const double r = std::fabs(a - b) / scale;
  rel = std::max(rel, r);
  return r <= TOL;
}

static void compare(const char* tc, const inputs& in, const char* mode,
                    bool propto, const outcome& ref, const outcome& br,
                    const std::vector<double>& cond) {
  bool ok = ref.what == br.what;
  double rel = 0, rel_plain = 0;
  if (ok && ref.what == "value") {
    ok = close(br.val, ref.val, std::max(1.0, std::fabs(ref.val)), rel)
         && br.grad.size() == ref.grad.size();
    for (size_t i = 0; ok && i < ref.grad.size(); ++i) {
      const double c = i < cond.size() ? cond[i] : 0;
      close(br.grad[i], ref.grad[i], std::max(1.0, std::fabs(ref.grad[i])),
            rel_plain);
      ok = close(br.grad[i], ref.grad[i],
                 std::max({1.0, std::fabs(ref.grad[i]), c}), rel)
           && ok;
    }
  }
  n_fail += !ok;
  n_ok += ok;
  std::printf("CHECK %s %-22s %s propto=%d %s rel %.2g rel_plain %.2g | %s\n",
              tc, in.name.c_str(), mode, propto, ok ? "OK" : "FAIL", rel,
              rel_plain,
              ref.what == "value" ? "value" : ref.what.substr(0, 60).c_str());
  if (!ok) {
    std::printf("  dev: %s val %.17g\n  br:  %s val %.17g\n", ref.what.c_str(),
                ref.val, br.what.c_str(), br.val);
  }
}

static std::vector<inputs> make_inputs(int ky, int km, int ks, int N) {
  std::mt19937 rng(42);
  std::normal_distribution<double> z(0, 1);
  const bool ys = nc::is_scalar_kind(ky), ms = nc::is_scalar_kind(km),
             ss = nc::is_scalar_kind(ks);
  VectorXd mu(N), y(N), sigma(N);
  for (int i = 0; i < N; ++i) {
    mu[i] = 0.3 + 2 * z(rng);
    y[i] = mu[i] + 0.7 * z(rng);
    sigma[i] = 0.7 * std::exp(0.3 * z(rng));
  }
  inputs base{"regular", ys ? y.head(1) : y, ms ? mu.head(1) : mu,
              ss ? VectorXd::Constant(1, 0.7) : sigma};
  const int py = ys ? 0 : 3, pm = ms ? 0 : 5, ps = ss ? 0 : 4;
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<inputs> v;
  auto add = [&](const char* name, auto f) {
    inputs c = base;
    c.name = name;
    f(c);
    v.push_back(c);
  };
  add("regular", [](inputs&) {});
  add("small_sigma", [](inputs& c) { c.sigma *= 1e-3; });
  add("sigma_inf", [&](inputs& c) { c.sigma[ps] = inf; });
  add("y_nan", [&](inputs& c) { c.y[py] = nan; });
  add("y_inf", [&](inputs& c) { c.y[py] = inf; });
  add("y_neg_inf", [&](inputs& c) { c.y[py] = -inf; });
  add("mu_inf", [&](inputs& c) { c.mu[pm] = -inf; });
  add("mu_nan", [&](inputs& c) { c.mu[pm] = nan; });
  add("sigma_zero", [&](inputs& c) { c.sigma[ps] = 0.0; });
  add("sigma_neg", [&](inputs& c) { c.sigma[ps] = -1.0; });
  add("sigma_nan", [&](inputs& c) { c.sigma[ps] = nan; });
  add("y_nan_sigma_neg", [&](inputs& c) {
    c.y[py] = nan;
    c.sigma[ps] = -1.0;
  });
  add("mu_inf_sigma_neg", [&](inputs& c) {
    c.mu[pm] = inf;
    c.sigma[ps] = -1.0;
  });
  add("y_inf_mu_inf", [&](inputs& c) {
    c.y[py] = inf;
    c.mu[pm] = inf;
  });
  add("y_inf_sigma_inf", [&](inputs& c) {
    c.y[py] = inf;
    c.sigma[ps] = inf;
  });
  add("overflow", [&](inputs& c) {
    c.y[py] = 1e300;
    c.mu[pm] = -1e300;
  });
  const int n_vec = !ys + !ms + !ss;
  if (n_vec > 0) {
    add("size_zero", [&](inputs& c) {
      if (!ys) c.y.resize(0);
      if (!ms) c.mu.resize(0);
      if (!ss) c.sigma.resize(0);
    });
    if (n_vec < 3) {
      add("size_zero_bad_scalar", [&](inputs& c) {
        if (!ys) c.y.resize(0);
        if (!ms) c.mu.resize(0);
        if (!ss) c.sigma.resize(0);
        if (ms) {
          c.mu[0] = nan;
        } else if (ss) {
          c.sigma[0] = -1.0;
        } else {
          c.y[0] = nan;
        }
      });
    }
  }
  if (n_vec >= 2) {
    add("size_mismatch", [&](inputs& c) {
      if (!ss) {
        c.sigma.conservativeResize(c.sigma.size() - 1);
      } else {
        c.mu.conservativeResize(c.mu.size() - 1);
      }
    });
  }
  return v;
}

template <bool P, int KY, int KM, int KS>
int heap_count_one(bool dev, const inputs& in) {
  int count = 0;
  for (int rep = 0; rep < 3; ++rep) {
    auto y = nc::make<KY, var>(in.y);
    auto mu = nc::make<KM, var>(in.mu);
    auto s = nc::make<KS, var>(in.sigma);
    g_big = 0;
    g_threshold = 8 * static_cast<size_t>(in.y.size() > 1
                                              ? in.y.size()
                                              : std::max(in.mu.size(),
                                                         in.sigma.size()));
    double v = 0;
    if (dev) {
      auto lp = call<true, P>(y, mu, s);
      if constexpr (std::is_same_v<decltype(lp), var>) {
        lp.grad();
        v = lp.val();
      } else {
        v = lp;
      }
    } else {
      auto lp = call<false, P>(y, mu, s);
      if constexpr (std::is_same_v<decltype(lp), var>) {
        lp.grad();
        v = lp.val();
      } else {
        v = lp;
      }
    }
    g_threshold = static_cast<size_t>(-1);
    count = g_big;
    sink = sink + v;
    stan::math::recover_memory();
  }
  return count;  // count of the last (warm) repetition
}

template <size_t I>
void run_case() {
  constexpr nc::type_case tc = nc::type_cases[I];
  constexpr int KY = tc.y, KM = tc.mu, KS = tc.sigma;
  const auto ins = make_inputs(KY, KM, KS, 1000);
  for (const inputs& in : ins) {
    const auto cond = conds(in, KY, KM, KS);
    compare(tc.name, in, "rev", true, run_rev<true, true, KY, KM, KS>(in),
            run_rev<false, true, KY, KM, KS>(in), cond);
    compare(tc.name, in, "rev", false, run_rev<true, false, KY, KM, KS>(in),
            run_rev<false, false, KY, KM, KS>(in), cond);
    if constexpr (!nc::is_varmat_kind(KY) && !nc::is_varmat_kind(KM)
                  && !nc::is_varmat_kind(KS)) {
      // tangents: one entry per autodiff argument; scale from its terms
      std::vector<double> fc;
      double c_loc = 0, c_scale = 0;
      {
        const auto full = conds(in, nc::VS, nc::VS, nc::VS);
        if (full.size() == 3) {
          c_loc = full[0];
          c_scale = full[2];
        }
      }
      if (nc::is_ad_kind(KY)) fc.push_back(c_loc);
      if (nc::is_ad_kind(KM)) fc.push_back(c_loc);
      if (nc::is_ad_kind(KS)) fc.push_back(c_scale);
      compare(tc.name, in, "fwd", true, run_fwd<true, true, KY, KM, KS>(in),
              run_fwd<false, true, KY, KM, KS>(in), fc);
      compare(tc.name, in, "fwd", false, run_fwd<true, false, KY, KM, KS>(in),
              run_fwd<false, false, KY, KM, KS>(in), fc);
    }
  }
  const inputs& reg = ins[0];
  std::printf("HEAP %s propto=1 dev %d branch %d\n", tc.name,
              heap_count_one<true, KY, KM, KS>(true, reg),
              heap_count_one<true, KY, KM, KS>(false, reg));
  std::printf("HEAP %s propto=0 dev %d branch %d\n", tc.name,
              heap_count_one<false, KY, KM, KS>(true, reg),
              heap_count_one<false, KY, KM, KS>(false, reg));
}

template <size_t... I>
void run_all(std::index_sequence<I...>) {
  (run_case<I>(), ...);
}

int main() {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  // self-test of the malloc counter
  g_big = 0;
  g_threshold = 8000;
  {
    VectorXd v = VectorXd::Zero(2000);
    sink = sink + v.sum();
  }
  g_threshold = static_cast<size_t>(-1);
  std::printf("HEAP_SELFTEST %d (expect 1)\n", g_big);
  run_all(std::make_index_sequence<std::size(nc::type_cases)>{});
  std::printf("CHECK_SUMMARY ok %d fail %d\n", n_ok, n_fail);
  return n_fail == 0 ? 0 : 1;
}
