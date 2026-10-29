// probe_hcub.cpp -- test whether hcubature() selects a box again after it
// split it, when every remaining box has error estimate 0.  Claim of
// stan-dev/math#3442: the two hcubature defects.
//
// The integrand is the st0-only integrand of the 7-parameter wiener_lpdf:
//   f(x) = exp(log p5(y - t0 - st0 x) - shift),  x in [0, 1].
// For each point and tolerance the probe prints:
//   stan   stan::math::hcubature
//   copy   a copy of its algorithm, with counters
//   fixed  the copy, with two changes: the loop condition uses the running
//          sum and error, and a split box is never selected again
//   ref    Boost Gauss-Kronrod 61, tolerance 1e-14
#include <stan/math/prim.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

struct Stats {
  int splits = 0;
  int dead_reselect = 0;
  int first_dead_at = -1;
  int live_zero_err = 0;
  int live = 0;
};

template <typename F, typename Pars>
double hcub_copy(const F& f, const Pars& pars, double a0, double b0,
                 int max_eval, double rel, bool fixed, Stats& st) {
  using stan::math::internal::gauss_kronrod;
  double result, err;
  std::tie(result, err) = gauss_kronrod(f, a0, b0, pars);
  int numevals = 15;
  double error = err, val = result;
  if (error <= rel * std::fabs(val) || numevals >= max_eval) {
    return result;
  }
  numevals += 30;
  std::vector<double> lo{a0}, hi{b0}, I{result}, E{err};
  std::vector<bool> dead{false};
  auto cond = [&]() {
    if (fixed) {
      return err > rel * std::fabs(result);
    }
    return error > rel * std::fabs(val);
  };
  while (numevals < max_eval && cond() && std::fabs(val) < INFINITY) {
    size_t k = std::distance(E.begin(), std::max_element(E.begin(), E.end()));
    if (dead[k]) {
      if (st.first_dead_at < 0) {
        st.first_dead_at = st.splits;
      }
      ++st.dead_reselect;
    }
    const double w = (hi[k] - lo[k]) / 2;
    double r1, e1, r2, e2;
    std::tie(r1, e1) = gauss_kronrod(f, lo[k] + w, hi[k], pars);
    std::tie(r2, e2) = gauss_kronrod(f, lo[k], hi[k] - w, pars);
    result += r1 + r2 - I[k];
    err += e1 + e2 - E[k];
    I[k] = 0;
    E[k] = fixed ? -1.0 : 0.0;
    dead[k] = true;
    const double l = lo[k], h = hi[k];
    lo.push_back(l + w);
    hi.push_back(h);
    I.push_back(r1);
    E.push_back(e1);
    dead.push_back(false);
    lo.push_back(l);
    hi.push_back(h - w);
    I.push_back(r2);
    E.push_back(e2);
    dead.push_back(false);
    numevals += 30;
    ++st.splits;
  }
  double sum = 0;
  for (size_t k = 0; k < I.size(); ++k) {
    sum += I[k];
    if (!dead[k]) {
      ++st.live;
      st.live_zero_err += (E[k] == 0.0);
    }
  }
  return sum;
}

void run(const char* id, double y, double a, double t0, double w, double v,
         double sv, double st0, bool use_shift) {
  using stan::math::internal::GradientCalc;
  using stan::math::internal::wiener5_density;
  const double lerr0 = std::log(1e-12) - std::log(2.0);
  const double t_max = std::fmin(1.0, (y - t0) / st0);
  const double shift
      = use_shift ? wiener5_density<GradientCalc::OFF>(
            y - t0 - st0 * t_max / 2, a, v, w, sv, std::log(1e-12))
                  : 0.0;
  const double lerr = lerr0 + shift;
  auto f = [](const auto& x, double y, double a, double v, double w,
              double t0, double sv, double st0, double lerr, double shift) {
    stan::scalar_seq_view<std::decay_t<decltype(x)>> xv(x);
    const double yd = y - (t0 + st0 * xv[0]);
    if (yd <= 0) {
      return 0.0;
    }
    return std::exp(wiener5_density<GradientCalc::OFF>(yd, a, v, w, sv, lerr)
                    - shift);
  };
  auto pars = std::make_tuple(y, a, v, w, t0, sv, st0, lerr, shift);
  auto g = [&](double x) { return f(x, y, a, v, w, t0, sv, st0, lerr, shift); };
  double gerr = 0;
  const double ref = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(
      g, 0.0, t_max, 20, 1e-14, &gerr);
  for (double prec : {1e-4, 1e-6, 1e-8, 1e-10}) {
    const double rel = 0.9 * prec / 2;
    Eigen::VectorXd lo = Eigen::VectorXd::Zero(1);
    Eigen::VectorXd hi = Eigen::VectorXd::Constant(1, t_max);
    const double stan_val = stan::math::hcubature(f, pars, 1, lo, hi, 6000, 0.0, rel);
    Stats sc, sf;
    const double copy_val = hcub_copy(f, pars, 0.0, t_max, 6000, rel, false, sc);
    const double fixed_val = hcub_copy(f, pars, 0.0, t_max, 6000, rel, true, sf);
    std::printf(
        "%s\tshift=%d\tprec=%.0e\tref=%.17g\tstan_rel=%.3e\tcopy_rel=%.3e"
        "\tfixed_rel=%.3e\tcopy_splits=%d\tdead_reselect=%d\tfirst_dead_at=%d"
        "\tlive=%d\tlive_zero_err=%d\tfixed_splits=%d\tfixed_dead=%d\n",
        id, use_shift, prec, ref, stan_val / ref - 1, copy_val / ref - 1,
        fixed_val / ref - 1, sc.splits, sc.dead_reselect, sc.first_dead_at,
        sc.live, sc.live_zero_err, sf.splits, sf.dead_reselect);
  }
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  for (bool s : {false, true}) {
    run("d17", 0.56925218279779977, 0.76009965905648347, 0.3210161323088358,
        0.68812100253640573, 2.9357113057129585, 1.1025981486567207,
        0.13881030510921391, s);
    run("d221", 0.66043406262219984, 0.83632055804413086, 0.24616775491872767,
        0.20675767209269905, 2.1502921448431254, 0.3834153712068723,
        0.25025043319223417, s);
  }
  return 0;
}
