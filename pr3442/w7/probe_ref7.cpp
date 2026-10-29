// probe_ref7.cpp -- C++ references for the 7-parameter wiener_lpdf, lcdf and
// lccdf (own adaptive Gauss-Kronrod; no hcubature).  Claims of
// stan-dev/math#3442: all 7-parameter accuracy numbers and test references.
//
// Point sets (environment): REF_D=1 the 300 D points; REF_P5=1 the p5 points
// (cdf only); REF_SWST0=1 the 240 points of swst0_points.hpp;
// REF_LIST="y/a/t0/w/v/sv/sw/st0:..." given points (ids l0, l1, ...); none
// of these: 4 D points and the extreme grid from index REF_START (default 0).
// Modifiers: REF_RANGE="a:b" (D and SWST0 subsets), REF_CDF=1 (lcdf and
// lccdf values), REF_CDF_FD=1 (lcdf and lccdf with finite-difference
// gradients), REF_SCORE=1 (d/dsw and d/dst0 under the integral).
//
// The 7-parameter density is the 5-parameter density averaged over a box:
//   p7 = 1 / (sw st0) int_{w - sw/2}^{w + sw/2} int_0^T p5(y - t0 - tau | om)
//        dtau dom,   T = min(st0, y - t0),
// with the integral over om (tau) left out if sw (st0) is 0. The probe
// integrates exp(log p5 - M) (M = largest log p5 on a grid) with its own
// adaptive Gauss-Kronrod 15 rule for vector integrands (global error
// queue, geometric start partition), nested for two dimensions. The log p5
// values and gradients come from the internal 5-parameter functions with a
// relative series tolerance (see eval5); the 7-parameter code and
// hcubature() are not used. Gradients:
//   d/dtheta log p7 = int p5 d/dtheta log p5 / int p5, theta in y a w v sv
//   d/dt0 = -d/dy
//   d/dsw = -1/sw + (G(w + sw/2) + G(w - sw/2)) / 2 / J,  G = int_0^T p5 dtau
//   d/dst0 = -1/st0 + int p5(y - t0 - st0 | om) dom / J   if st0 < y - t0
//            -1/st0                                      otherwise
// Output: REF <id> <y a t0 w v sv sw st0> <value> <g_y g_a g_t0 g_w g_v g_sv
// g_sw g_st0> <boxes outer> <max boxes inner> <rel err est>, or FAIL <msg>.
#include <stan/math/rev.hpp>
#include "swst0_points.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <random>
#include <string>
#include <vector>

using stan::math::var;

namespace {

using Vec = std::array<double, 6>;

// Gauss-Kronrod 15 nodes and weights (QUADPACK qk15)
const double XGK[8]
    = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
       0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
       0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
       0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
const double WGK[8]
    = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
       0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
       0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
       0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
const double WG[4]
    = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
       0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

struct Box {
  double a, b;
  Vec I, E;
};

Box gk15(const std::function<Vec(double)>& f, double a, double b) {
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  Box bx{a, b, {}, {}};
  Vec K{}, G{};
  const Vec fc = f(c);
  for (int m = 0; m < 6; ++m) {
    K[m] = WGK[7] * fc[m];
    G[m] = WG[3] * fc[m];
  }
  for (int k = 0; k < 7; ++k) {
    const Vec f1 = f(c - h * XGK[k]);
    const Vec f2 = f(c + h * XGK[k]);
    for (int m = 0; m < 6; ++m) {
      K[m] += WGK[k] * (f1[m] + f2[m]);
      if (k % 2 == 1) {
        G[m] += WG[k / 2] * (f1[m] + f2[m]);
      }
    }
  }
  for (int m = 0; m < 6; ++m) {
    bx.I[m] = h * K[m];
    bx.E[m] = std::fabs(h * (K[m] - G[m]));
  }
  return bx;
}

// adaptive: split the box with the largest scaled error until the summed
// error of every component is below tol * scale, scale = max(|I_m|, 1e-6
// I_0 |typical gradient|) -- here max(|I_m|, 1e-6 |I_0|)
Vec integrate(const std::function<Vec(double)>& f, double a, double b,
              double tol, int max_boxes, int* n_boxes, double* rel_err) {
  // start from a partition that is geometric toward both ends, so that a
  // peak of width down to 2^-40 (b - a) at an end is resolved
  std::vector<double> cuts{a, b};
  for (int k = 2; k <= 40; ++k) {
    const double d = (b - a) * std::ldexp(1.0, -k);
    cuts.push_back(a + d);
    cuts.push_back(b - d);
  }
  cuts.push_back(0.5 * (a + b));
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  std::vector<Box> boxes;
  for (size_t i = 0; i + 1 < cuts.size(); ++i) {
    boxes.push_back(gk15(f, cuts[i], cuts[i + 1]));
  }
  max_boxes += static_cast<int>(boxes.size());
  auto totals = [&](Vec& I, Vec& E) {
    I.fill(0.0);
    E.fill(0.0);
    for (const auto& bx : boxes) {
      for (int m = 0; m < 6; ++m) {
        I[m] += bx.I[m];
        E[m] += bx.E[m];
      }
    }
  };
  Vec I, E;
  while (true) {
    totals(I, E);
    Vec scale;
    for (int m = 0; m < 6; ++m) {
      scale[m] = std::fmax(std::fabs(I[m]), 1e-6 * std::fabs(I[0]));
      scale[m] = std::fmax(scale[m], 1e-300);
    }
    double worst = 0;
    for (int m = 0; m < 6; ++m) {
      worst = std::fmax(worst, E[m] / scale[m]);
    }
    *rel_err = worst;
    if (worst <= tol || static_cast<int>(boxes.size()) >= max_boxes) {
      break;
    }
    size_t k = 0;
    double kmax = -1;
    for (size_t i = 0; i < boxes.size(); ++i) {
      double s = 0;
      for (int m = 0; m < 6; ++m) {
        s = std::fmax(s, boxes[i].E[m] / scale[m]);
      }
      if (s > kmax) {
        kmax = s;
        k = i;
      }
    }
    const double lo = boxes[k].a, hi = boxes[k].b, mid = 0.5 * (lo + hi);
    boxes[k] = gk15(f, lo, mid);
    boxes.push_back(gk15(f, mid, hi));
  }
  *n_boxes = static_cast<int>(boxes.size());
  return I;
}

// log p5 and its gradients in y, a, w, v, sv. The internal functions get a
// tolerance relative to the factor in front of the series: the public
// 5-parameter wiener_lpdf uses the absolute tolerance 1e-12, which truncates
// the series too early where the density is small.
Vec eval5(double t, double a, double w, double v, double sv) {
  using stan::math::internal::GradientCalc;
  namespace in = stan::math::internal;
  Vec r;
  if (!(t > 0)) {
    r.fill(0.0);
    r[0] = -INFINITY;
    return r;
  }
  const double let = in::wiener5_compute_log_error_term(t, a, v, w, sv);
  const double err = std::log(1e-15) + std::fmin(let, 0.0);
  r[0] = in::wiener5_density<GradientCalc::OFF>(t, a, v, w, sv, err);
  r[1] = in::wiener5_grad_t<GradientCalc::OFF>(t, a, v, w, sv, err);
  r[2] = in::wiener5_grad_a<GradientCalc::OFF>(t, a, v, w, sv, err);
  r[3] = in::wiener5_grad_w<GradientCalc::OFF>(t, a, v, w, sv, err);
  r[4] = in::wiener5_grad_v<GradientCalc::OFF>(t, a, v, w, sv, err);
  r[5] = in::wiener5_grad_sv<GradientCalc::OFF>(t, a, v, w, sv, err);
  return r;
}

void ref(const std::string& id, const std::vector<double>& p) {
  const double y = p[0], a = p[1], t0 = p[2], w = p[3], v = p[4], sv = p[5],
               sw = p[6], st0 = p[7];
  std::printf("REF\t%s", id.c_str());
  for (double q : p) {
    std::printf("\t%.17g", q);
  }
  try {
    const double T = std::fmin(st0, y - t0);
    const double wlo = w - sw / 2, whi = w + sw / 2;
    double M = -INFINITY;
    for (int i = 0; i <= 40; ++i) {
      const double om = sw > 0 ? wlo + sw * i / 40.0 : w;
      for (int j = 0; j <= 40; ++j) {
        const double tau = st0 > 0 ? T * j / 40.0 : 0.0;
        M = std::fmax(M, eval5(y - t0 - tau, a, om, v, sv)[0]);
        if (st0 == 0) {
          break;
        }
      }
      if (sw == 0) {
        break;
      }
    }
    if (!std::isfinite(M)) {
      throw std::runtime_error("no finite log p5 on the grid");
    }
    // integrand at (tau, om): [weight, weighted gradients y a w v sv]
    auto f = [&](double tau, double om) {
      Vec r = eval5(y - t0 - tau, a, om, v, sv);
      const double e = std::isfinite(r[0]) ? std::exp(r[0] - M) : 0.0;
      Vec out;
      out[0] = e;
      for (int m = 1; m < 6; ++m) {
        out[m] = e == 0.0 ? 0.0 : e * r[m];
      }
      return out;
    };
    int inner_max = 0, nb = 0;
    double err_in = 0, err_out = 0;
    auto G = [&](double om) {
      if (st0 == 0) {
        return f(0.0, om);
      }
      int n = 0;
      double e = 0;
      Vec r = integrate([&](double tau) { return f(tau, om); }, 0.0, T, 1e-11,
                        100, &n, &e);
      inner_max = std::max(inner_max, n);
      err_in = std::fmax(err_in, e);
      return r;
    };
    Vec J = sw == 0 ? G(w)
                    : integrate(G, wlo, whi, 1e-10, 100, &nb, &err_out);
    const double J0 = J[0];
    double value = M + std::log(J0);
    if (sw > 0) {
      value -= std::log(sw);
    }
    if (st0 > 0) {
      value -= std::log(st0);
    }
    std::array<double, 8> g;
    g[0] = J[1] / J0;
    g[1] = J[2] / J0;
    g[2] = -g[0];
    g[3] = J[3] / J0;
    g[4] = J[4] / J0;
    g[5] = J[5] / J0;
    g[6] = sw > 0 ? -1 / sw + 0.5 * (G(whi)[0] + G(wlo)[0]) / J0 : 0.0;
    if (st0 == 0) {
      g[7] = 0.0;
    } else if (st0 < y - t0) {
      double edge;
      if (sw == 0) {
        edge = f(st0, w)[0];
      } else {
        int n = 0;
        double e = 0;
        edge = integrate([&](double om) { return f(st0, om); }, wlo, whi,
                         1e-10, 100, &n, &e)[0];
      }
      g[7] = -1 / st0 + edge / J0;
    } else {
      g[7] = -1 / st0;
    }
    // REF_SCORE = "1": d/dsw and d/dst0 in the score form (differentiation
    // under the integral), which has no -1/sw or -1/st0 cancellation:
    //   d/dsw  = int (om - w) / sw d/dw p5 / J
    //   d/dst0 = -int tau / st0 d/dt p5 / J
    // (the second is the endpoint formula after integration by parts; p5 = 0
    // at t = 0). The endpoint values are appended after the error columns.
    const double g6_end = g[6], g7_end = g[7];
    const bool score = std::getenv("REF_SCORE") != nullptr;
    if (score && sw > 0) {
      auto h = [&](double om) {
        Vec r = G(om);
        Vec o{};
        o[0] = (om - w) / sw * r[3];
        return o;
      };
      int n = 0;
      double e = 0;
      g[6] = integrate(h, wlo, whi, 1e-10, 100, &n, &e)[0] / J0;
      err_out = std::fmax(err_out, e);
    }
    if (score && st0 > 0) {
      auto inner = [&](double om) {
        int n = 0;
        double e = 0;
        Vec r = integrate(
            [&](double tau) {
              Vec q = f(tau, om);
              Vec o{};
              o[0] = tau / st0 * q[1];
              return o;
            },
            0.0, T, 1e-11, 100, &n, &e);
        err_in = std::fmax(err_in, e);
        return r;
      };
      double num;
      if (sw == 0) {
        num = inner(w)[0];
      } else {
        int n = 0;
        double e = 0;
        num = integrate(inner, wlo, whi, 1e-10, 100, &n, &e)[0];
        err_out = std::fmax(err_out, e);
      }
      g[7] = -num / J0;
    }
    std::printf("\t%.17g", value);
    for (double q : g) {
      std::printf("\t%.17g", q);
    }
    std::printf("\t%d\t%d\t%.2g", nb, inner_max, std::fmax(err_in, err_out));
    if (score) {
      std::printf("\t%.17g\t%.17g", g6_end, g7_end);
    }
    std::printf("\n");
  } catch (const std::exception& e) {
    std::printf("\tFAIL\t%s\n", e.what());
  }
}

// Reference values of the 7-parameter lcdf and lccdf (unnormalized, upper
// boundary). With d = y - t0 and T = min(st0, d):
//   F7 = 1/(sw st0) int dom int_0^d p5(s | om) min(T, d - s) ds
//   P7 - F7 = 1/(sw st0) int dom int_{d - st0}^inf p5(s | om)
//             clamp(s - d + st0, 0, st0) ds
// (for st0 = 0 the weights are 1 on [0, d] and on [d, inf)). The second
// form has no cancellation. The upper limit is d + 160 a^2 / pi^2, where
// p5 has decayed by more than exp(-80). Output: RCDF <id> <lcdf> <lccdf>.
std::pair<double, double> ref_cdf_values(const std::vector<double>& p);

void ref_cdf(const std::string& id, const std::vector<double>& p) {
  const auto r = ref_cdf_values(p);
  std::printf("RCDF\t%s\t%.17g\t%.17g\n", id.c_str(), r.first, r.second);
}

// Central finite differences of the reference lcdf and lccdf in each of the
// 8 arguments, step h = 1e-5 max(|x|, 0.05). Output: RFD <id> <lcdf>
// <lccdf> <8 lcdf gradients> <8 lccdf gradients>.
void ref_cdf_fd(const std::string& id, const std::vector<double>& p) {
  const auto r0 = ref_cdf_values(p);
  std::vector<double> g_cdf(8), g_ccdf(8);
  for (int k = 0; k < 8; ++k) {
    if ((k == 6 || k == 7) && p[k] == 0) {
      g_cdf[k] = g_ccdf[k] = NAN;
      continue;
    }
    const double h = 1e-5 * std::fmax(std::fabs(p[k]), 0.05);
    std::vector<double> pp = p, pm = p;
    pp[k] += h;
    pm[k] -= h;
    const auto rp = ref_cdf_values(pp);
    const auto rm = ref_cdf_values(pm);
    g_cdf[k] = (rp.first - rm.first) / (2 * h);
    g_ccdf[k] = (rp.second - rm.second) / (2 * h);
  }
  std::printf("RFD\t%s\t%.17g\t%.17g", id.c_str(), r0.first, r0.second);
  for (double q : g_cdf) {
    std::printf("\t%.12g", q);
  }
  for (double q : g_ccdf) {
    std::printf("\t%.12g", q);
  }
  std::printf("\n");
}

std::pair<double, double> ref_cdf_values(const std::vector<double>& p) {
  using stan::math::internal::GradientCalc;
  namespace in = stan::math::internal;
  const double y = p[0], a = p[1], t0 = p[2], w = p[3], v = p[4], sv = p[5],
               sw = p[6], st0 = p[7];
  const double d = y - t0;
  const double T = std::fmin(st0, d);
  const double wlo = w - sw / 2, whi = w + sw / 2;
  const double s_max = d + 160.0 * a * a / (M_PI * M_PI);
  auto lp5 = [&](double s, double om) -> double {
    if (!(s > 0)) {
      return -INFINITY;
    }
    const double let = in::wiener5_compute_log_error_term(s, a, v, om, sv);
    return in::wiener5_density<GradientCalc::OFF>(
        s, a, v, om, sv, std::log(1e-15) + std::fmin(let, 0.0));
  };
  // shift: largest log p5 on a grid over [0, s_max] x [wlo, whi]
  double M = -INFINITY;
  for (int i = 0; i <= 20; ++i) {
    const double om = sw > 0 ? wlo + sw * i / 20.0 : w;
    for (int j = 1; j <= 400; ++j) {
      M = std::fmax(M, lp5(s_max * j / 400.0, om));
    }
    if (sw == 0) {
      break;
    }
  }
  auto weight_cdf = [&](double s) {
    return st0 > 0 ? std::fmin(T, d - s) : 1.0;
  };
  auto weight_ccdf = [&](double s) {
    return st0 > 0 ? std::fmin(std::fmax(s - d + st0, 0.0), st0) : 1.0;
  };
  // The weights have a kink (at s = d - T for the cdf, at s = d for the
  // ccdf); the integral is split there, as the adaptive rule converges
  // slowly across a kink.
  auto inner_piece = [&](double om, double lo, double hi, bool cdf) {
    int n = 0;
    double e = 0;
    return integrate(
        [&](double s) {
          Vec r{};
          const double l = lp5(s, om);
          r[0] = std::isfinite(l)
                     ? std::exp(l - M) * (cdf ? weight_cdf(s) : weight_ccdf(s))
                     : 0.0;
          return r;
        },
        lo, hi, 1e-11, 100, &n, &e)[0];
  };
  auto inner = [&](double om, double lo, double hi, bool cdf) {
    const double kink = cdf ? d - T : d;
    if (st0 > 0 && kink > lo && kink < hi) {
      return inner_piece(om, lo, kink, cdf) + inner_piece(om, kink, hi, cdf);
    }
    return inner_piece(om, lo, hi, cdf);
  };
  auto outer = [&](bool cdf) {
    const double lo = cdf ? 0.0 : std::fmax(0.0, d - st0);
    const double hi = cdf ? d : s_max;
    if (sw == 0) {
      return inner(w, lo, hi, cdf);
    }
    int n = 0;
    double e = 0;
    return integrate(
        [&](double om) {
          Vec r{};
          r[0] = inner(om, lo, hi, cdf);
          return r;
        },
        wlo, whi, 1e-10, 100, &n, &e)[0];
  };
  double norm = 0;
  if (sw > 0) {
    norm += std::log(sw);
  }
  if (st0 > 0) {
    norm += std::log(st0);
  }
  return {M + std::log(outer(true)) - norm, M + std::log(outer(false)) - norm};
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  // REF_SWST0 = "1": the points of swst0_points.hpp (small sw and st0);
  // REF_RANGE = "a:b" selects points a .. b-1 of that list
  if (std::getenv("REF_SWST0")) {
    const auto pts = swst0_points();
    size_t lo = 0, hi = pts.size();
    if (const char* range = std::getenv("REF_RANGE")) {
      const std::string s(range);
      const size_t c = s.find(':');
      lo = std::atoi(s.substr(0, c).c_str());
      hi = std::min(hi, static_cast<size_t>(std::atoi(s.substr(c + 1).c_str())));
    }
    for (size_t i = lo; i < hi; ++i) {
      ref(pts[i].first, pts[i].second);
    }
    return 0;
  }
  // REF_LIST = "y/a/t0/w/v/sv/sw/st0:y/a/...": compute only these points
  // (ids l0, l1, ...)
  if (const char* list = std::getenv("REF_LIST")) {
    std::string s(list);
    size_t pos = 0;
    int id = 0;
    while (pos < s.size()) {
      size_t end = s.find(':', pos);
      if (end == std::string::npos) {
        end = s.size();
      }
      std::vector<double> p;
      size_t q = pos;
      while (q < end) {
        size_t r = s.find('/', q);
        if (r == std::string::npos || r > end) {
          r = end;
        }
        p.push_back(std::strtod(s.substr(q, r - q).c_str(), nullptr));
        q = r + 1;
      }
      if (p.size() == 8) {
        // REF_CDF_FD = "1": cdf references with finite-difference gradients
        if (std::getenv("REF_CDF_FD")) {
          ref_cdf_fd("l" + std::to_string(id), p);
        } else {
          ref("l" + std::to_string(id), p);
        }
      }
      ++id;
      pos = end + 1;
    }
    return 0;
  }
  // REF_D = "1": the 300 D points of probe_shift.cpp (ids d0 .. d299), only
  if (std::getenv("REF_D")) {
    std::mt19937 rng(20261009);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int k = 0; k < 300; ++k) {
      const double aa = 0.5 + 2.5 * U(rng);
      const double vv = -3.0 + 6.0 * U(rng);
      const double ww = 0.2 + 0.6 * U(rng);
      const double tt0 = 0.1 + 0.3 * U(rng);
      const double svv = 1.5 * U(rng);
      const double sww
          = (k % 3 == 2) ? 0.0 : 0.8 * U(rng) * 2.0 * std::fmin(ww, 1 - ww);
      const double stt0 = (k % 3 == 1) ? 0.0 : 0.3 * U(rng);
      const double yy = tt0 + 0.05 + 2.0 * U(rng);
      // REF_RANGE = "a:b": only points a .. b-1
      if (const char* range = std::getenv("REF_RANGE")) {
        const std::vector<double> r = [&] {
          std::vector<double> out;
          std::string s(range);
          const size_t c = s.find(':');
          out.push_back(std::atof(s.substr(0, c).c_str()));
          out.push_back(std::atof(s.substr(c + 1).c_str()));
          return out;
        }();
        if (k < r[0] || k >= r[1]) {
          continue;
        }
      }
      const std::vector<double> p{yy, aa, tt0, ww, vv, svv, sww, stt0};
      // REF_CDF = "1": the lcdf and lccdf references instead; with
      // REF_CDF_FD = "1" also their finite-difference gradients
      if (std::getenv("REF_CDF_FD")) {
        ref_cdf_fd("d" + std::to_string(k), p);
      } else if (std::getenv("REF_CDF")) {
        ref_cdf("d" + std::to_string(k), p);
      } else {
        ref("d" + std::to_string(k), p);
      }
    }
    return 0;
  }
  // REF_P5 = "1": cdf references with finite-difference gradients at the p5
  // points of probe_w4.cpp (ids p5<k>)
  if (std::getenv("REF_P5")) {
    std::mt19937 rng(20261010);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int s = 0; s < 2; ++s) {
      for (int k = 0; k < 300; ++k) {
        const double a = 0.5 + 2.5 * U(rng);
        const double v = -5.0 + 10.0 * U(rng);
        const double w = 0.2 + 0.6 * U(rng);
        const double t = 0.05 + 3.0 * U(rng);
        const double sv = s == 1 ? 1.5 * U(rng) : 0.0;
        if (s == 1) {
          ref_cdf_fd("p5" + std::to_string(k),
                     {0.1 + t, a, 0.1, w, v, sv, 0.0, 0.0});
        }
      }
    }
    return 0;
  }
  // REF_START = first E index to compute; the D points only if it is 0
  const char* start_env = std::getenv("REF_START");
  const int start = start_env ? std::atoi(start_env) : 0;
  if (start == 0) {
  // D points where develop or the branch was wrong at precision 1e-10
  ref("d17", {0.56925218279779977, 0.76009965905648347, 0.3210161323088358,
              0.68812100253640573, 2.9357113057129585, 1.1025981486567207, 0,
              0.13881030510921391});
  ref("d55", {0.80675126136439368, 2.5310625226778387, 0.36755347829230678,
              0.46748130089354978, -1.4643994108608576, 0.75455166358992853,
              0.58154398269106056, 0});
  ref("d221", {0.66043406262219984, 0.83632055804413086, 0.24616775491872767,
               0.20675767209269905, 2.1502921448431254, 0.3834153712068723, 0,
               0.25025043319223417});
  ref("d246", {1.6075726778128945, 2.7699834918376043, 0.27332973706351155,
               0.38629891015529771, 2.5628925485873895, 0.99808830122225589,
               0.55409674526881225, 0.26327625131261939});
  }
  // the E grid of probe_shift.cpp
  int k = 0;
  for (double yd : {1e-6, 1e-4, 10.0, 50.0}) {
    for (double aa : {0.1, 8.0}) {
      for (double vv : {-30.0, -10.0, 10.0, 30.0}) {
        for (double ww : {1e-4, 0.5}) {
          for (double stt0 : {0.0, 1e-3}) {
            const double sww = 0.5 * std::fmin(ww, 1 - ww);
            if (k >= start) {
              ref("e" + std::to_string(k),
                  {0.2 + yd, aa, 0.2, ww, vv, 0.1, sww, stt0});
            }
            ++k;
          }
        }
      }
    }
  }
  return 0;
}
