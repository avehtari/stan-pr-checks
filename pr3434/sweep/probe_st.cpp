// probe_st.cpp --- student_t_lpdf, multi_student_t_cholesky_lpdf,
// multi_student_t_lpdf, yule_simon_lpmf and lbeta with stan::math::var at
// the points of a hex point file (rows 10 to 13 of the worst-error table).
//
// Usage: probe_st <points file>
//
// Input line:   <TAG> <hex> <hex> ...      ('#' lines are copied: regions)
// Output line:  <TAG> <hex args...> => <hex value> <hex grads...> <status>
// status: '.' finished, 'T' the call threw (every number is nan then).
//
// Tags and argument order (the gradients follow the same order):
//   ST    y nu mu sigma       student_t_lpdf
//   MST   y1 y2 nu            multi_student_t_cholesky_lpdf, mu = 0, L = I
//   MSTF  y1 y2 nu            multi_student_t_lpdf, mu = 0, Sigma = I
//   YS    n alpha             yule_simon_lpmf (gradient in alpha only)
//   LBETA d x                 lbeta (not a density)
// MST and MSTF print d/dy1, d/dy2, d/dnu. Every call uses propto = false.

#include <stan/math/rev.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::var;

namespace {

int n_args_of(const std::string& tag) {
  if (tag == "ST") {
    return 4;
  }
  if (tag == "MST" || tag == "MSTF") {
    return 3;
  }
  if (tag == "YS" || tag == "LBETA") {
    return 2;
  }
  return -1;
}

int n_grads_of(const std::string& tag) {
  return tag == "YS" ? 1 : n_args_of(tag);
}

void eval_one(const std::string& tag, const std::vector<double>& a,
              double& value, std::vector<double>& grads) {
  if (tag == "ST") {
    var y = a[0], nu = a[1], mu = a[2], sigma = a[3];
    var lp = stan::math::student_t_lpdf(y, nu, mu, sigma);
    lp.grad();
    value = lp.val();
    grads = {y.adj(), nu.adj(), mu.adj(), sigma.adj()};
  } else if (tag == "MST" || tag == "MSTF") {
    Matrix<var, Dynamic, 1> y(2);
    y(0) = a[0];
    y(1) = a[1];
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(2);
    Eigen::MatrixXd ident = Eigen::MatrixXd::Identity(2, 2);
    var nu = a[2];
    var lp = (tag == "MST")
                 ? stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, ident)
                 : stan::math::multi_student_t_lpdf(y, nu, mu, ident);
    lp.grad();
    value = lp.val();
    grads = {y(0).adj(), y(1).adj(), nu.adj()};
  } else if (tag == "YS") {
    var alpha = a[1];
    var lp = stan::math::yule_simon_lpmf(static_cast<int>(a[0]), alpha);
    lp.grad();
    value = lp.val();
    grads = {alpha.adj()};
  } else {
    var d = a[0], x = a[1];
    var lb = stan::math::lbeta(d, x);
    lb.grad();
    value = lb.val();
    grads = {d.adj(), x.adj()};
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: probe_st <points file>\n");
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '#') {
      std::printf("%s\n", line.c_str());
      continue;
    }
    std::istringstream ss(line);
    std::string tag;
    ss >> tag;
    const int n_args = n_args_of(tag);
    if (n_args < 0) {
      std::fprintf(stderr, "unknown tag %s\n", tag.c_str());
      return 2;
    }
    std::vector<double> a;
    std::string token;
    while (ss >> token) {
      a.push_back(std::strtod(token.c_str(), nullptr));
    }
    if (static_cast<int>(a.size()) != n_args) {
      std::fprintf(stderr, "%s: %zu args, want %d\n", tag.c_str(), a.size(),
                   n_args);
      return 2;
    }
    double value = std::nan("");
    std::vector<double> grads(n_grads_of(tag), std::nan(""));
    bool threw = false;
    try {
      eval_one(tag, a, value, grads);
    } catch (const std::exception&) {
      threw = true;
      value = std::nan("");
      grads.assign(n_grads_of(tag), std::nan(""));
    }
    stan::math::recover_memory();
    std::printf("%s", tag.c_str());
    for (double x : a) {
      std::printf(" %a", x);
    }
    std::printf(" => %a", value);
    for (double g : grads) {
      std::printf(" %a", g);
    }
    std::printf(" %c\n", threw ? 'T' : '.');
  }
  return 0;
}
