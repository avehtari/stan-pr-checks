// Time of one gradient evaluation of normal_lpdf variants.  One binary per
// variant: -DNV=<0..7> (see normal_variants.hpp; 6 glm_k1 = normal_id_glm
// with K = 1, 7 matvec_k1 = its two double products).
// Output: RESULT <variant> <N> <median_ns> <min_ns> <reps>

#include "normal_variants.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

#ifndef NV
#error "define NV"
#endif

constexpr const char* names[]
    = {"dev", "f1", "f2", "f3", "f4", "compute_only", "glm_k1", "matvec_k1"};

struct data_t {
  Eigen::VectorXd y, eta;
  Eigen::MatrixXd x;
  explicit data_t(int N) : y(N), eta(N), x(N, 1) {
    std::mt19937 rng(7919u * N);
    std::normal_distribution<double> z(0, 1);
    for (int i = 0; i < N; ++i) {
      eta[i] = 0.3 + z(rng);
      y[i] = eta[i] + 0.7 * z(rng);
      x(i, 0) = eta[i];
    }
  }
};

volatile double sink = 0;

static double one(const data_t& d) {
  double v;
  if constexpr (NV == 7) {
    Eigen::VectorXd e = d.x * Eigen::VectorXd::Constant(1, 1.0);
    Eigen::VectorXd g = d.x.transpose() * d.y;
    return e[0] + g[0];
  } else {
    {
      stan::math::var sigma(0.7);
      if constexpr (NV == 6) {
        Eigen::Matrix<stan::math::var, -1, 1> beta(1);
        beta[0] = 1.0;
        stan::math::var lp
            = stan::math::normal_id_glm_lpdf<true>(d.y, d.x, 0.0, beta, sigma);
        lp.grad();
        v = lp.val() + beta[0].adj();
      } else {
        nv::vm mu(d.eta);
        stan::math::var lp;
        if constexpr (NV == 0) lp = nv::dev(d.y, mu, sigma);
        if constexpr (NV == 1) lp = nv::f1(d.y, mu, sigma);
        if constexpr (NV == 2) lp = nv::f2(d.y, mu, sigma);
        if constexpr (NV == 3) lp = nv::f3(d.y, mu, sigma);
        if constexpr (NV == 4) lp = nv::f4(d.y, mu, sigma);
        if constexpr (NV == 5) lp = nv::compute_only(d.y, mu, sigma);
        lp.grad();
        v = lp.val() + mu.adj()[0];
      }
      v += sigma.adj();
    }
    stan::math::recover_memory();
    return v;
  }
}

static double time_batch(const data_t& d, int reps) {
  const auto t0 = std::chrono::steady_clock::now();
  double acc = 0;
  for (int r = 0; r < reps; ++r) {
    acc += one(d);
  }
  const auto t1 = std::chrono::steady_clock::now();
  sink = sink + acc;
  return std::chrono::duration<double, std::nano>(t1 - t0).count() / reps;
}

int main() {
  for (int N : {100, 1000, 10000, 100000}) {
    data_t d(N);
    double t = time_batch(d, 3);
    int reps = std::max(1, static_cast<int>(2e7 / t));
    t = time_batch(d, reps);
    reps = std::max(1, static_cast<int>(2e7 / t));
    std::vector<double> ts;
    for (int b = 0; b < 7; ++b) {
      ts.push_back(time_batch(d, reps));
    }
    std::sort(ts.begin(), ts.end());
    std::printf("RESULT %s %d %.1f %.1f %d\n", names[NV], N, ts[3], ts[0],
                reps);
    std::fflush(stdout);
  }
  std::printf("SINK %g\n", static_cast<double>(sink));
  return 0;
}
