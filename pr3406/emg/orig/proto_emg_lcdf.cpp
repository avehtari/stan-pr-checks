// Prototype (PLAN 2026-10-20, step 1) of exp_mod_normal_lcdf through erfcx
// differences, value and the four partials, mu = 0, sigma = 1, lambda = v.
// Output in the format of probe_emg_phi.cpp ("E lcdf z v value dy dmu dsigma
// dlambda"), so score_emg_phi.py and score_emg_detail.py score it; the
// branch per point goes to stderr.
//
// z = (y - mu) / sigma, v = lambda sigma, u = -z / sqrt 2, d = v / sqrt 2,
// t = u + d, E = erfcx, g = exp(-z^2 / 2) / 2, phi = phi(z) = g sqrt(2 / pi).
//   Phi(z) = g E(u), b = exp(v^2/2 - v z) Phi(z - v) = g E(t),
//   F = g Delta(u, d), Delta(a, d) = E(a) - E(a + d);  P = 1 - F.
//   F_z = v b,  F_v = phi - (v - z) b = -g E'(t) / sqrt 2,
//   sigma F_sigma = -z F_z + v F_v = v (phi - v b).
#include <stan/math/prim.hpp>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace sm = stan::math;
namespace smi = stan::math::internal;

constexpr double SQRT2 = 1.41421356237309504880;
constexpr double SQRT_PI = 1.77245385090551602730;
constexpr double SQRT_2_OVER_PI = 0.79788456080286535588;

double E(double x) { return sm::erfcx(x); }
double dE(double x, double ex) { return smi::erfcx_derivative(x, ex); }

// Cody tail: E(x) = (1 + h(x)) / (sqrt(pi) x), h = w T(w), w = 1 / x^2,
// T = sqrt(pi) P(w) / Q(w), P and Q ascending in w (erfcx_tail_p, _q).
// For x, y >= 4: h(y) - h(x) without cancellation, through the divided
// differences of A(w) = w P(w) and Q(w):
//   h(y) - h(x) = sqrt(pi) (w1 - w0) [DD_A Q(w0) - A(w0) DD_Q] / (Q(w1) Q(w0)).
double poly(const double* c, int n, double w) {
  double r = c[n - 1];
  for (int k = n - 2; k >= 0; --k) {
    r = r * w + c[k];
  }
  return r;
}
// (A(w1) - A(w0)) / (w1 - w0) for A(w) = sum c[k] w^k.
double divdiff(const double* c, int n, double w1, double w0) {
  double s = 0.0;
  for (int k = 1; k < n; ++k) {
    double t = 0.0, p1 = 1.0;
    for (int j = 0; j < k; ++j) {  // sum_j w1^j w0^(k-1-j)
      t += p1 * std::pow(w0, k - 1 - j);
      p1 *= w1;
    }
    s += c[k] * t;
  }
  return s;
}
double h_tail(double x) {
  const double w = 1.0 / (x * x);
  const double* p = smi::erfcx_tail_p.data();
  const double* q = smi::erfcx_tail_q.data();
  return SQRT_PI * w * poly(p, 6, w) / poly(q, 6, w);
}
// h(y) - h(x) with y - x = delta given exactly: a rounded y - x would carry
// the rounding error of y into a difference of size delta.
double h_diff(double x, double y, double delta) {
  const double w0 = 1.0 / (x * x), w1 = 1.0 / (y * y);
  const double dw = -delta * (y + x) * w0 * w1;
  const double* p = smi::erfcx_tail_p.data();
  const double* q = smi::erfcx_tail_q.data();
  double a[7];  // A(w) = w P(w)
  a[0] = 0.0;
  for (int k = 0; k < 6; ++k) {
    a[k + 1] = p[k];
  }
  const double a0 = poly(a, 7, w0), q0 = poly(q, 6, w0), q1 = poly(q, 6, w1);
  const double num = divdiff(a, 7, w1, w0) * q0 - a0 * divdiff(q, 6, w1, w0);
  return SQRT_PI * dw * num / (q1 * q0);
}

// E^(1), E^(3), E^(5), E^(7) at m by E^(k+1) = 2 m E^(k) + 2 k E^(k-1).
void odd_derivs(double m, double out[4]) {
  double e[9];
  e[0] = E(m);
  e[1] = dE(m, e[0]);
  for (int k = 1; k < 8; ++k) {
    e[k + 1] = 2.0 * m * e[k] + 2.0 * k * e[k - 1];
  }
  out[0] = e[1];
  out[1] = e[3];
  out[2] = e[5];
  out[3] = e[7];
}

// log Delta = log(E(lo) - E(hi)), hi - lo = d > 0 given exactly.
// upper_exact: hi is the exact argument (lo = hi - d is derived); else lo is.
// The differences use d itself, never a rounded hi - lo.
// which: 0 direct, 1 midpoint series, 2 tail form.
double log_delta_core(double lo, double hi, double d, bool upper_exact,
                      int& which) {
  const double elo = E(lo);
  if (lo >= 4.0) {
    which = 2;
    // R = E(hi) / E(lo) = (lo / hi) (1 + h(hi)) / (1 + h(lo)).
    const double log_ratio_x
        = upper_exact ? std::log1p(-d / hi) : -std::log1p(d / lo);
    const double log_r
        = log_ratio_x + std::log1p(h_diff(lo, hi, d) / (1.0 + h_tail(lo)));
    return std::log(elo) + sm::log1m_exp(log_r);
  }
  const double ehi = E(hi);
  if (ehi / elo < 0.99) {
    which = 0;
    return std::log(elo - ehi);
  }
  which = 1;
  const double m = upper_exact ? hi - 0.5 * d : lo + 0.5 * d;
  double o[4];
  odd_derivs(m, o);
  const double h = 0.5 * d;
  const double h2 = h * h;
  const double s
      = o[0] + h2 * (o[1] / 6.0 + h2 * (o[2] / 120.0 + h2 * o[3] / 5040.0));
  return std::log(-2.0 * h * s);
}
// E(a) - E(a + d), a exact.
double log_delta(double a, double d, int& which) {
  return log_delta_core(a, a + d, d, false, which);
}
// E(c - d) - E(c), c exact.
double log_delta_upper(double c, double d, int& which) {
  return log_delta_core(c - d, c, d, true, which);
}

struct Out {
  double value, dz, dv, ds;      // lcdf; ds = sigma * d/dsigma
  double cvalue, cdz, cdv, cds;  // lccdf
  int branch;
};

// lccdf from the lcdf where F <= 1/2: P = 1 - F >= 1/2, and each lccdf
// partial is minus the lcdf partial times F / P.
void complement_small_f(Out& o) {
  const double f = std::exp(o.value);
  o.cvalue = std::log1p(-f);
  const double ratio = std::exp(o.value - o.cvalue);
  o.cdz = -o.dz * ratio;
  o.cdv = -o.dv * ratio;
  o.cds = -o.ds * ratio;
}

Out proto(double z, double v) {
  const double u = -z / SQRT2, d = v / SQRT2, t = u + d;
  const double log_g = -0.5 * z * z - M_LN2;
  const double g = std::exp(log_g);
  const double phi = g * SQRT_2_OVER_PI;
  // log b: for t >= 0 as log g + log E(t), which has no cancellation; the
  // form v^2/2 - v z + log Phi(z - v) adds terms of size v^2/2.
  const double log_b = t >= 0.0 ? log_g + std::log(E(t))
                                 : v * (0.5 * v - z) + sm::std_normal_lcdf(z - v);
  const double f_est = std::exp(sm::std_normal_lcdf(z)) - std::exp(log_b);
  Out o;
  if (z <= 0.0 || (t >= 0.0 && f_est <= 0.5 && z < 37.0)) {
    // F = g Delta(u, d)
    int which;
    const double ld = log_delta(u, d, which);
    o.branch = which;
    o.value = log_g + ld;
    const double inv_delta = std::exp(-ld);
    const double et = E(t);
    const double det = dE(t, et);
    o.dz = v * et * inv_delta;
    o.dv = -det / SQRT2 * inv_delta;
    // phi - v b = g sqrt 2 (u E(t) - E'(t) / 2): the accurate E' replaces
    // 1 / sqrt(pi) - t E(t), which cancels for large t.
    o.ds = v * SQRT2 * (u * et - 0.5 * det) * inv_delta;
    complement_small_f(o);
  } else if (t < 0.0 && f_est <= 0.5) {
    // z > 0, v < z: F = -expm1(d (2u + d)) + g Delta(-t, d), both > 0.
    int which;
    const double ld = log_delta_upper(-u, d, which);
    o.branch = 4 + which;
    const double first = -std::expm1(d * (2.0 * u + d));
    const double second = std::exp(log_g + ld);
    const double f = first + second;
    o.value = std::log(f);
    // b = g E(t) = exp(t^2 - u^2) - g E(-t), E(t) = 2 exp(t^2) - E(-t).
    const double b = std::exp(d * (2.0 * u + d)) - g * E(-t);
    o.dz = v * b / f;
    o.dv = (phi - (v - z) * b) / f;
    o.ds = v * (phi - v * b) / f;
    complement_small_f(o);
  } else {
    // F > 1/2: log F = log1m_exp(log P), P = Phi(-z) + b.
    o.branch = 3;
    const double log_q = sm::std_normal_lcdf(-z);
    const double log_p = sm::log_sum_exp(log_q, log_b);
    o.value = sm::log1m_exp(log_p);
    o.cvalue = log_p;
    // X / F for the lcdf and -X / P for the lccdf, with the tiny factors
    // (b, g, phi) divided in log space: b and P can both be below DBL_MIN.
    const double b_f = std::exp(log_b - o.value), b_p = std::exp(log_b - log_p);
    o.dz = v * b_f;
    o.cdz = -v * b_p;
    if (t >= 0.0) {
      const double et = E(t), det = dE(t, et);
      const double g_f = std::exp(log_g - o.value), g_p = std::exp(log_g - log_p);
      o.dv = -g_f * det / SQRT2;
      o.cdv = g_p * det / SQRT2;
      o.ds = v * g_f * SQRT2 * (u * et - 0.5 * det);
      o.cds = -v * g_p * SQRT2 * (u * et - 0.5 * det);
    } else {
      const double log_phi = log_g + std::log(SQRT_2_OVER_PI);
      const double phi_f = std::exp(log_phi - o.value), phi_p = std::exp(log_phi - log_p);
      o.dv = phi_f - (v - z) * b_f;
      o.cdv = -(phi_p - (v - z) * b_p);
      o.ds = v * (phi_f - v * b_f);
      o.cds = -v * (phi_p - v * b_p);
    }
  }
  return o;
}

void emit(double z, double v) {
  const Out o = proto(z, v);
  // mu = 0, sigma = 1, lambda = v: dy = dz, dmu = -dz, dsigma = ds,
  // dlambda = sigma dv.
  std::printf("E lcdf %a %a %a %a %a %a %a\n", z, v, o.value, o.dz, -o.dz,
              o.ds, o.dv);
  std::printf("E lccdf %a %a %a %a %a %a %a\n", z, v, o.cvalue, o.cdz,
              -o.cdz, o.cds, o.cdv);
  std::fprintf(stderr, "%g %g branch %d\n", z, v, o.branch);
}

int main(int argc, char** argv) {
  const bool dense = argc > 1 && std::string(argv[1]) == "dense";
  std::vector<double> zs, vs;
  if (!dense) {
    zs = {-1e4, -1e3, -300, -100, -30, -10, -5, -3, -1, -0.3, 0, 0.3, 1, 3, 10, 30, 100};
    vs = {1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1, 3, 10, 100};
  } else {
    for (double lz = -4; lz <= 2.01; lz += 0.25) {
      zs.push_back(-std::pow(10.0, lz));
      zs.push_back(std::pow(10.0, lz));
    }
    for (double lv = -10; lv <= 2.01; lv += 0.5) {
      vs.push_back(std::pow(10.0, lv));
    }
  }
  for (double z : zs) {
    for (double v : vs) {
      emit(z, v);
    }
  }
  return 0;
}
