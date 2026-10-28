// probe_throws.cpp --- for each point of a point file, which of the gradient
// roots throws: inc_beta_dda, inc_beta_ddb, the 9-argument
// grad_reg_inc_beta, and beta(a, b).  Each call is made on its own.
//
// Usage: probe_throws <points file>
// Output: one line per point with at least one throw, then counts.

#include <stan/math/prim.hpp>

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>

namespace {

bool throws(const std::function<void()>& f, std::string* msg) {
  try {
    f();
    return false;
  } catch (const std::exception& e) {
    *msg = e.what();
    if (msg->size() > 60) {
      *msg = msg->substr(0, 60);
    }
    return true;
  }
}

}  // namespace

int main(int argc, char** argv) {
  FILE* fh = std::fopen(argc > 1 ? argv[1] : "points_beta.txt", "r");
  if (fh == nullptr) {
    return 1;
  }
  using stan::math::digamma;
  long n_points = 0;
  long count[4] = {0, 0, 0, 0};
  char line[512];
  while (std::fgets(line, sizeof(line), fh)) {
    char ta[128], tb[128], tz[128];
    if (line[0] == '#'
        || std::sscanf(line, "%127s %127s %127s", ta, tb, tz) != 3) {
      continue;
    }
    ++n_points;
    const double a = std::strtod(ta, nullptr);
    const double b = std::strtod(tb, nullptr);
    const double z = std::strtod(tz, nullptr);
    std::string m[4];
    const bool t[4] = {
        throws(
            [&] {
              stan::math::inc_beta_dda(a, b, z, digamma(a), digamma(a + b));
            },
            &m[0]),
        throws(
            [&] {
              stan::math::inc_beta_ddb(a, b, z, digamma(b), digamma(a + b));
            },
            &m[1]),
        throws(
            [&] {
              double g1 = 0;
              double g2 = 0;
              stan::math::grad_reg_inc_beta(g1, g2, a, b, z, digamma(a),
                                            digamma(b), digamma(a + b),
                                            stan::math::beta(a, b));
            },
            &m[2]),
        throws([&] { stan::math::beta(a, b); }, &m[3])};
    if (t[0] || t[1] || t[2] || t[3]) {
      std::printf("%.17g %.17g %.17g  dda=%d ddb=%d root9=%d beta=%d  %s\n",
                  a, b, z, t[0], t[1], t[2], t[3],
                  (t[0] ? m[0] : t[1] ? m[1] : t[2] ? m[2] : m[3]).c_str());
    }
    for (int i = 0; i < 4; ++i) {
      count[i] += t[i];
    }
  }
  std::fclose(fh);
  std::printf("# points %ld  throws: dda %ld  ddb %ld  root9 %ld  beta %ld\n",
              n_points, count[0], count[1], count[2], count[3]);
  return 0;
}
