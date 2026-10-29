// probe_sib.cpp --- evaluate Stan Math distribution functions and their
// gradients with stan::math::var at the points of a hex point file.
//
// Large-shape spot check on develop (6110e533f3). Copied from
// ../../density_sweep_work/probe_sweep.cpp.
//
// Usage: probe_sib <points file>
//
// Input line:   <TAG> <LABEL> <hex> <hex> ...
// Output line:  <TAG> <LABEL> <hex args...> => <hex value> <hex grads...>
//               <status> <milliseconds>
//
// status: '.' finished, 'T' the call threw (every number is nan then).
// Every call uses the default (propto = false) overload.
// Argument and gradient order per tag: the SHAPES table below and the
// ARITY table in ref_sib.py. Do not change one without the other.

#include <stan/math/rev.hpp>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::var;

namespace {

struct Shape {
  const char* tag;
  int n_args;
  int n_grads;
};

const Shape SHAPES[] = {
    {"NB", 3, 2},    {"NB2", 3, 2},   {"NB2L", 3, 2},  {"GLM1", 5, 3},
    {"GLM5", 13, 3}, {"BNB", 4, 3},   {"DM", 6, 3},    {"LKJ", 4, 1},
    {"LKJC", 6, 1},  {"BBC", 4, 2},   {"BBLC", 4, 2},  {"BBLCC", 4, 2},
    {"YSC", 2, 1},   {"YSLC", 2, 1},  {"YSLCC", 2, 1}, {"BNBC", 4, 3},
    {"BNBLC", 4, 3}, {"BNBLCC", 4, 3},
};

const Shape* shape_of(const std::string& tag) {
  for (const Shape& s : SHAPES) {
    if (tag == s.tag) {
      return &s;
    }
  }
  return nullptr;
}

int to_int(double x) { return static_cast<int>(x); }

// ── One evaluation per tag ───────────────────────────────────────────────

void eval_one(const std::string& tag, const std::vector<double>& a,
              double& value, std::vector<double>& grads) {
  using namespace stan::math;
  var lp;
  std::vector<var> wrt;
  if (tag == "NB" || tag == "NB2" || tag == "NB2L") {
    var p1 = a[1], p2 = a[2];
    if (tag == "NB") {
      lp = neg_binomial_lpmf(to_int(a[0]), p1, p2);
    } else if (tag == "NB2") {
      lp = neg_binomial_2_lpmf(to_int(a[0]), p1, p2);
    } else {
      lp = neg_binomial_2_log_lpmf(to_int(a[0]), p1, p2);
    }
    wrt = {p1, p2};
  } else if (tag == "GLM1" || tag == "GLM5") {
    // GLM1: y, x, alpha, beta, phi. GLM5: y1..y5, x1..x5, alpha, beta, phi.
    const int n_obs = (tag == "GLM1") ? 1 : 5;
    std::vector<int> y(n_obs);
    Eigen::MatrixXd x(n_obs, 1);
    for (int i = 0; i < n_obs; ++i) {
      y[i] = to_int(a[i]);
      x(i, 0) = a[n_obs + i];
    }
    var alpha = a[2 * n_obs];
    Matrix<var, Dynamic, 1> beta(1);
    beta(0) = a[2 * n_obs + 1];
    var phi = a[2 * n_obs + 2];
    lp = neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi);
    wrt = {alpha, beta(0), phi};
  } else if (tag == "BNB" || tag == "BNBC" || tag == "BNBLC"
             || tag == "BNBLCC") {
    var r = a[1], alpha = a[2], beta = a[3];
    const int n = to_int(a[0]);
    if (tag == "BNB") {
      lp = beta_neg_binomial_lpmf(n, r, alpha, beta);
    } else if (tag == "BNBC") {
      lp = beta_neg_binomial_cdf(n, r, alpha, beta);
    } else if (tag == "BNBLC") {
      lp = beta_neg_binomial_lcdf(n, r, alpha, beta);
    } else {
      lp = beta_neg_binomial_lccdf(n, r, alpha, beta);
    }
    wrt = {r, alpha, beta};
  } else if (tag == "DM") {
    std::vector<int> ns{to_int(a[0]), to_int(a[1]), to_int(a[2])};
    Matrix<var, Dynamic, 1> alpha(3);
    alpha(0) = a[3];
    alpha(1) = a[4];
    alpha(2) = a[5];
    lp = dirichlet_multinomial_lpmf(ns, alpha);
    wrt = {alpha(0), alpha(1), alpha(2)};
  } else if (tag == "LKJ") {
    // r12, r13, r23, eta
    Eigen::MatrixXd omega(3, 3);
    omega << 1.0, a[0], a[1], a[0], 1.0, a[2], a[1], a[2], 1.0;
    var eta = a[3];
    lp = lkj_corr_lpdf(omega, eta);
    wrt = {eta};
  } else if (tag == "LKJC") {
    // L10, L11, L20, L21, L22, eta; L00 = 1
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(3, 3);
    L(0, 0) = 1.0;
    L(1, 0) = a[0];
    L(1, 1) = a[1];
    L(2, 0) = a[2];
    L(2, 1) = a[3];
    L(2, 2) = a[4];
    var eta = a[5];
    lp = lkj_corr_cholesky_lpdf(L, eta);
    wrt = {eta};
  } else if (tag == "BBC" || tag == "BBLC" || tag == "BBLCC") {
    var alpha = a[2], beta = a[3];
    const int n = to_int(a[0]), N = to_int(a[1]);
    if (tag == "BBC") {
      lp = beta_binomial_cdf(n, N, alpha, beta);
    } else if (tag == "BBLC") {
      lp = beta_binomial_lcdf(n, N, alpha, beta);
    } else {
      lp = beta_binomial_lccdf(n, N, alpha, beta);
    }
    wrt = {alpha, beta};
  } else if (tag == "YSC" || tag == "YSLC" || tag == "YSLCC") {
    var alpha = a[1];
    const int n = to_int(a[0]);
    if (tag == "YSC") {
      lp = yule_simon_cdf(n, alpha);
    } else if (tag == "YSLC") {
      lp = yule_simon_lcdf(n, alpha);
    } else {
      lp = yule_simon_lccdf(n, alpha);
    }
    wrt = {alpha};
  } else {
    std::fprintf(stderr, "unknown tag %s\n", tag.c_str());
    std::exit(2);
  }
  lp.grad();
  value = lp.val();
  grads.clear();
  for (const var& v : wrt) {
    grads.push_back(v.adj());
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: probe_sib <points file>\n");
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream ss(line);
    std::string tag, label;
    ss >> tag >> label;
    const Shape* shape = shape_of(tag);
    if (shape == nullptr) {
      std::fprintf(stderr, "unknown tag %s\n", tag.c_str());
      return 2;
    }
    std::vector<double> a;
    std::string token;
    while (ss >> token) {
      a.push_back(std::strtod(token.c_str(), nullptr));
    }
    if (static_cast<int>(a.size()) != shape->n_args) {
      std::fprintf(stderr, "%s: %zu args, want %d\n", tag.c_str(), a.size(),
                   shape->n_args);
      return 2;
    }
    double value = std::nan("");
    std::vector<double> grads(shape->n_grads, std::nan(""));
    bool threw = false;
    const auto t0 = std::chrono::steady_clock::now();
    try {
      eval_one(tag, a, value, grads);
    } catch (const std::exception& e) {
      threw = true;
      value = std::nan("");
      grads.assign(shape->n_grads, std::nan(""));
    }
    const auto t1 = std::chrono::steady_clock::now();
    stan::math::recover_memory();
    const double ms
        = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("%s %s", tag.c_str(), label.c_str());
    for (double x : a) {
      std::printf(" %a", x);
    }
    std::printf(" => %a", value);
    for (double g : grads) {
      std::printf(" %a", g);
    }
    std::printf(" %c %.3f\n", threw ? 'T' : '.', ms);
    std::fflush(stdout);
  }
  return 0;
}
