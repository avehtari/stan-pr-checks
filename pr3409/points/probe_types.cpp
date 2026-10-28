// probe_types.cpp --- shape gradients of inc_beta(a, b, z) in each autodiff
// type, with both shapes or one shape as autodiff arguments, against the
// double root grad_reg_inc_beta (9-argument form).  Uses only functions
// that exist on develop.
//
// Usage: probe_types <points file>
// Output: one line per region with the largest relative difference to the
// double root and the number of points that differ, per configuration.

#include <stan/math/mix.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

namespace {

using stan::math::fvar;
using stan::math::var;
using ffvar = fvar<fvar<double>>;

const char* kNames[] = {"var2",   "var_a",  "var_b",   "fvar2",
                        "fvar_a", "fvar_b", "ffvar2",  "ffvar_a",
                        "ffvar_b"};
constexpr int kConf = 9;

// out[c][0] = d/da, out[c][1] = d/db; NaN where the configuration does not
// give that derivative
void eval(double a, double b, double z, double out[kConf][2]) {
  for (int c = 0; c < kConf; ++c) {
    out[c][0] = out[c][1] = std::nan("");
  }
  {
    stan::math::nested_rev_autodiff nested;
    var av = a, bv = b;
    var r = stan::math::inc_beta(av, bv, z);
    r.grad();
    out[0][0] = av.adj();
    out[0][1] = bv.adj();
  }
  {
    stan::math::nested_rev_autodiff nested;
    var av = a;
    var r = stan::math::inc_beta(av, b, z);
    r.grad();
    out[1][0] = av.adj();
  }
  {
    stan::math::nested_rev_autodiff nested;
    var bv = b;
    var r = stan::math::inc_beta(a, bv, z);
    r.grad();
    out[2][1] = bv.adj();
  }
  for (int k = 0; k < 2; ++k) {
    fvar<double> af(a, k == 0 ? 1.0 : 0.0), bf(b, k == 1 ? 1.0 : 0.0);
    out[3][k] = stan::math::inc_beta(af, bf, fvar<double>(z, 0.0)).d_;
  }
  out[4][0] = stan::math::inc_beta(fvar<double>(a, 1.0), b, z).d_;
  out[5][1] = stan::math::inc_beta(a, fvar<double>(b, 1.0), z).d_;
  for (int k = 0; k < 2; ++k) {
    ffvar af = a, bf = b;
    (k == 0 ? af : bf).d_.val_ = 1.0;
    out[6][k] = stan::math::inc_beta(af, bf, ffvar(z)).d_.val_;
  }
  {
    ffvar af = a;
    af.d_.val_ = 1.0;
    out[7][0] = stan::math::inc_beta(af, b, z).d_.val_;
  }
  {
    ffvar bf = b;
    bf.d_.val_ = 1.0;
    out[8][1] = stan::math::inc_beta(a, bf, z).d_.val_;
  }
}

struct Stat {
  double worst = 0;
  long differ = 0;
  long nonfinite = 0;
};

}  // namespace

int main(int argc, char** argv) {
  const char* path = (argc > 1) ? argv[1] : "points_beta_sub.txt";
  FILE* fh = std::fopen(path, "r");
  if (fh == nullptr) {
    std::printf("cannot open %s\n", path);
    return 1;
  }
  std::vector<std::string> order;
  std::map<std::string, std::vector<Stat>> stats;
  std::map<std::string, long> count;
  std::string region = "?";
  char line[512];
  while (std::fgets(line, sizeof(line), fh)) {
    if (line[0] == '#') {
      region = std::string(line + 2);
      region.erase(region.find_last_not_of("\n") + 1);
      order.push_back(region);
      stats[region].assign(kConf, Stat());
      continue;
    }
    char ta[128], tb[128], tz[128];
    if (std::sscanf(line, "%127s %127s %127s", ta, tb, tz) != 3) {
      continue;
    }
    const double a = std::strtod(ta, nullptr);
    const double b = std::strtod(tb, nullptr);
    const double z = std::strtod(tz, nullptr);
    double ref[2] = {0, 0};
    double out[kConf][2];
    try {
      using stan::math::digamma;
      stan::math::grad_reg_inc_beta(ref[0], ref[1], a, b, z, digamma(a),
                                    digamma(b), digamma(a + b),
                                    stan::math::beta(a, b));
      eval(a, b, z, out);
    } catch (const std::exception&) {
      continue;
    }
    ++count[region];
    for (int c = 0; c < kConf; ++c) {
      Stat& s = stats[region][c];
      for (int k = 0; k < 2; ++k) {
        if ((c % 3 == 1 && k == 1) || (c % 3 == 2 && k == 0)) {
          continue;  // the configuration does not give this derivative
        }
        const double v = out[c][k];
        if (!std::isfinite(v) || !std::isfinite(ref[k])) {
          if (!(std::isnan(v) && std::isnan(ref[k]))) {
            ++s.nonfinite;
          }
          continue;
        }
        if (v != ref[k]) {
          ++s.differ;
          const double m = std::max(std::fabs(v), std::fabs(ref[k]));
          s.worst = std::max(s.worst, std::fabs(v - ref[k]) / m);
        }
      }
    }
  }
  std::fclose(fh);
  std::printf("%-26s %5s", "region", "n");
  for (int c = 0; c < kConf; ++c) {
    std::printf(" %15s", kNames[c]);
  }
  std::printf("\n");
  for (const auto& r : order) {
    std::printf("%-26.26s %5ld", r.c_str(), count[r]);
    for (int c = 0; c < kConf; ++c) {
      const Stat& s = stats[r][c];
      std::printf(" %7.1e/%3ld/%2ld", s.worst, s.differ, s.nonfinite);
    }
    std::printf("\n");
  }
  std::printf("cell: largest relative difference / points that differ / "
              "non-finite where the root is finite or the reverse\n");
  return 0;
}
