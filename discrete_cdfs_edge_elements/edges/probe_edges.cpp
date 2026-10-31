// Probe behind the table of the PR text: for each discrete cdf, lcdf and
// lccdf, f({ordinary, edge}) and f({edge, ordinary}) against f({ordinary})
// (edge with probability 1) or against 0 / -inf with zero gradient (edge
// with probability 0). Reverse mode, scalar parameters.
#include <stan/math/rev.hpp>
#include <climits>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using stan::math::var;
namespace sm = stan::math;

using Fn = std::function<var(const std::vector<int>&, const std::vector<var>&)>;

struct Res {
  double val;
  std::vector<double> grad;
  bool threw;
};

static Res eval(const Fn& f, const std::vector<int>& n,
                const std::vector<double>& p) {
  Res r{0, {}, false};
  std::vector<var> pv(p.begin(), p.end());
  try {
    var v = f(n, pv);
    r.val = v.val();
    v.grad();
    for (auto& x : pv) {
      r.grad.push_back(x.adj());
    }
  } catch (const std::exception& e) {
    r.threw = true;
  }
  sm::recover_memory();
  return r;
}

static bool same(double a, double b) {
  if (std::isinf(a) || std::isinf(b)) {
    return a == b;
  }
  return std::fabs(a - b) <= 1e-12 * std::fmax(1.0, std::fabs(b));
}

static int n_bad = 0;

// kind "one": the edge element has probability 1
// kind "zero": the edge element has probability 0; z is the expected value
static void check(const std::string& name, const Fn& f, int ordinary, int edge,
                  const std::vector<double>& p, const std::string& kind,
                  double z) {
  const Res alone = eval(f, {ordinary}, p);
  for (int order = 0; order < 2; ++order) {
    const std::vector<int> n
        = order == 0 ? std::vector<int>{ordinary, edge}
                     : std::vector<int>{edge, ordinary};
    const Res both = eval(f, n, p);
    bool ok = !both.threw && !alone.threw;
    double expect = kind == "one" ? alone.val : z;
    double max_dg = 0;
    bool any_nan = false;
    if (ok) {
      ok = same(both.val, expect);
      for (size_t k = 0; k < both.grad.size(); ++k) {
        const double target = kind == "one" ? alone.grad[k] : 0.0;
        const double dg = std::fabs(both.grad[k] - target);
        any_nan = any_nan || std::isnan(dg);
        max_dg = std::fmax(max_dg, dg);
        ok = ok && (kind == "one" ? same(both.grad[k], target)
                                  : both.grad[k] == 0.0);
      }
    }
    n_bad += !ok;
    std::printf("%-26s %-4s edge=%-11d %-9s %s  value %-24.17g expected %-24.17g"
                "  max|grad diff| %.3g%s%s\n",
                name.c_str(), kind.c_str(), edge,
                order == 0 ? "ord,edge" : "edge,ord", ok ? "OK " : "BAD",
                both.val, expect, max_dg, any_nan ? " NaN" : "",
                both.threw ? "  THROW" : "");
  }
}

int main() {
  const double inf = std::numeric_limits<double>::infinity();
  const int M = INT_MAX;
#define F1(fn) \
  [](const std::vector<int>& n, const std::vector<var>& p) { return sm::fn(n, p[0]); }
#define F2(fn)                                              \
  [](const std::vector<int>& n, const std::vector<var>& p) { \
    return sm::fn(n, p[0], p[1]);                           \
  }
#define F3(fn)                                              \
  [](const std::vector<int>& n, const std::vector<var>& p) { \
    return sm::fn(n, p[0], p[1], p[2]);                     \
  }
#define FBIN(fn)                                            \
  [](const std::vector<int>& n, const std::vector<var>& p) { \
    return sm::fn(n, 10, p[0]);                             \
  }
#define FBB(fn)                                             \
  [](const std::vector<int>& n, const std::vector<var>& p) { \
    return sm::fn(n, 10, p[0], p[1]);                       \
  }
#define FPB1(fn)                                                 \
  [](const std::vector<int>& n, const std::vector<var>& p) {      \
    Eigen::Matrix<var, -1, 1> th(1);                             \
    th << p[0];                                                  \
    return sm::fn(n, std::vector<Eigen::Matrix<var, -1, 1>>(n.size(), th)); \
  }
#define FPB(fn)                                                  \
  [](const std::vector<int>& n, const std::vector<var>& p) {      \
    Eigen::Matrix<var, -1, 1> th(3);                             \
    th << p[0], p[1], p[2];                                      \
    return sm::fn(n, std::vector<Eigen::Matrix<var, -1, 1>>(n.size(), th)); \
  }

  const std::vector<double> th{0.3}, ab{2.5, 1.5}, mp{4.0, 2.0},
      rab{3.0, 2.5, 1.5}, ys{1.5}, pb{0.5, 0.2, 0.7};

  check("bernoulli_cdf", F1(bernoulli_cdf), 0, 1, th, "one", 0);
  check("bernoulli_cdf", F1(bernoulli_cdf), 0, -1, th, "zero", 0);
  check("bernoulli_lcdf", F1(bernoulli_lcdf), 0, 1, th, "one", 0);
  check("bernoulli_lcdf", F1(bernoulli_lcdf), 0, -1, th, "zero", -inf);
  check("bernoulli_lccdf", F1(bernoulli_lccdf), 0, -1, th, "one", 0);
  check("bernoulli_lccdf", F1(bernoulli_lccdf), 0, 1, th, "zero", -inf);

  check("binomial_cdf", FBIN(binomial_cdf), 3, 10, th, "one", 0);
  check("binomial_cdf", FBIN(binomial_cdf), 3, -1, th, "zero", 0);
  check("binomial_lcdf", FBIN(binomial_lcdf), 3, 10, th, "one", 0);
  check("binomial_lcdf", FBIN(binomial_lcdf), 3, -1, th, "zero", -inf);
  check("binomial_lccdf", FBIN(binomial_lccdf), 3, -1, th, "one", 0);
  check("binomial_lccdf", FBIN(binomial_lccdf), 3, 10, th, "zero", -inf);

  check("beta_binomial_cdf", FBB(beta_binomial_cdf), 3, 10, ab, "one", 0);
  check("beta_binomial_cdf", FBB(beta_binomial_cdf), 3, -1, ab, "zero", 0);
  check("beta_binomial_lcdf", FBB(beta_binomial_lcdf), 3, 10, ab, "one", 0);
  check("beta_binomial_lcdf", FBB(beta_binomial_lcdf), 3, -1, ab, "zero", -inf);
  check("beta_binomial_lccdf", FBB(beta_binomial_lccdf), 3, -1, ab, "one", 0);
  check("beta_binomial_lccdf", FBB(beta_binomial_lccdf), 3, 10, ab, "zero",
        -inf);

  check("neg_binomial_cdf", F2(neg_binomial_cdf), 3, M, ab, "one", 0);
  check("neg_binomial_cdf", F2(neg_binomial_cdf), 3, -1, ab, "zero", 0);
  check("neg_binomial_lcdf", F2(neg_binomial_lcdf), 3, M, ab, "one", 0);
  check("neg_binomial_lcdf", F2(neg_binomial_lcdf), 3, -1, ab, "zero", -inf);
  check("neg_binomial_lccdf", F2(neg_binomial_lccdf), 3, -1, ab, "one", 0);
  check("neg_binomial_lccdf", F2(neg_binomial_lccdf), 3, M, ab, "zero", -inf);

  check("neg_binomial_2_cdf", F2(neg_binomial_2_cdf), 3, M, mp, "one", 0);
  check("neg_binomial_2_cdf", F2(neg_binomial_2_cdf), 3, -1, mp, "zero", 0);
  check("neg_binomial_2_lcdf", F2(neg_binomial_2_lcdf), 3, M, mp, "one", 0);
  check("neg_binomial_2_lcdf", F2(neg_binomial_2_lcdf), 3, -1, mp, "zero",
        -inf);
  check("neg_binomial_2_lccdf", F2(neg_binomial_2_lccdf), 3, -1, mp, "one", 0);
  check("neg_binomial_2_lccdf", F2(neg_binomial_2_lccdf), 3, M, mp, "zero",
        -inf);

  check("beta_neg_binomial_cdf", F3(beta_neg_binomial_cdf), 3, M, rab, "one",
        0);
  check("beta_neg_binomial_cdf", F3(beta_neg_binomial_cdf), 3, -1, rab, "zero",
        0);
  check("beta_neg_binomial_lcdf", F3(beta_neg_binomial_lcdf), 3, M, rab, "one",
        0);
  check("beta_neg_binomial_lcdf", F3(beta_neg_binomial_lcdf), 3, -1, rab,
        "zero", -inf);
  check("beta_neg_binomial_lccdf", F3(beta_neg_binomial_lccdf), 3, -1, rab,
        "one", 0);
  check("beta_neg_binomial_lccdf", F3(beta_neg_binomial_lccdf), 3, M, rab,
        "zero", -inf);

  check("yule_simon_cdf", F1(yule_simon_cdf), 3, M, ys, "one", 0);
  check("yule_simon_cdf", F1(yule_simon_cdf), 3, 0, ys, "zero", 0);
  check("yule_simon_lcdf", F1(yule_simon_lcdf), 3, M, ys, "one", 0);
  check("yule_simon_lcdf", F1(yule_simon_lcdf), 3, 0, ys, "zero", -inf);
  check("yule_simon_lccdf", F1(yule_simon_lccdf), 3, 0, ys, "one", 0);
  check("yule_simon_lccdf", F1(yule_simon_lccdf), 3, M, ys, "zero", -inf);

  check("poisson_binomial_cdf", FPB(poisson_binomial_cdf), 1, 3, pb, "one", 0);
  check("poisson_binomial_lcdf", FPB(poisson_binomial_lcdf), 1, 3, pb, "one",
        0);
  check("poisson_binomial_lccdf", FPB(poisson_binomial_lccdf), 1, 3, pb, "zero",
        -inf);
  // One trial per element: the develop code adds -inf in the main loop
  check("poisson_binomial_lccdf_1", FPB1(poisson_binomial_lccdf), 0, 1, th,
        "zero", -inf);

  // The numbers that the PR text quotes
  {
    Res r = eval(F2(neg_binomial_lccdf), {3, -1}, ab);
    Res a = eval(F2(neg_binomial_lccdf), {3}, ab);
    std::printf("CLAIM neg_binomial_lccdf({3, -1} | 2.5, 1.5) = %.6g;"
                " neg_binomial_lccdf(3 | 2.5, 1.5) = %.6g\n",
                r.val, a.val);
  }
  {
    Res r = eval(F2(neg_binomial_2_lcdf), {3, M}, mp);
    Res a = eval(F2(neg_binomial_2_lcdf), {3}, mp);
    std::printf("CLAIM neg_binomial_2_lcdf({3, INT_MAX} | 4, 2): value %.6g,"
                " d/dmu %.6g, d/dphi %.6g;"
                " neg_binomial_2_lcdf(3 | 4, 2): value %.6g, d/dmu %.6g,"
                " d/dphi %.6g\n",
                r.val, r.grad[0], r.grad[1], a.val, a.grad[0], a.grad[1]);
  }
  {
    Eigen::VectorXd th3(3);
    th3 << 0.5, 0.2, 0.7;
    const double v = sm::poisson_binomial_lccdf(3, th3);
    const double l53 = std::log(std::ldexp(1.0, -53));
    std::printf("CLAIM poisson_binomial_lccdf(3 | 0.5, 0.2, 0.7) = %.17g;"
                " log(2^-53) = %.17g; equal: %s\n",
                v, l53, v == l53 ? "yes" : "no");
  }
  {
    // discrete_range has no parameters; values only
    std::printf("CLAIM discrete_range_lccdf({2, 0} | 1, 5) = %.6g,"
                " discrete_range_lccdf(2 | 1, 5) = %.6g;"
                " discrete_range_lcdf({2, 9} | 1, 5) = %.6g,"
                " discrete_range_lcdf(2 | 1, 5) = %.6g\n",
                sm::discrete_range_lccdf(std::vector<int>{2, 0}, 1, 5),
                sm::discrete_range_lccdf(2, 1, 5),
                sm::discrete_range_lcdf(std::vector<int>{2, 9}, 1, 5),
                sm::discrete_range_lcdf(2, 1, 5));
  }

  std::printf("BAD rows: %d\n", n_bad);
  return 0;
}
