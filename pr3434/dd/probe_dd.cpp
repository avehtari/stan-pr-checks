// probe_dd.cpp --- digamma_diff(x, d) against the plain digamma difference.
//
// Usage: probe_dd <points file>
// Input line:   <x hex> <d hex>
// Output line:  <x hex> <d hex> => <digamma_diff> <plain difference>
//               <d/dx of digamma_diff by reverse mode> <status>
//
// Built against the branch tree only (digamma_diff does not exist on
// develop).

#include <stan/math/rev.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: probe_dd <points file>\n");
    return 2;
  }
  std::ifstream in(argv[1]);
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
    std::string tx, td;
    ss >> tx >> td;
    const double x = std::strtod(tx.c_str(), nullptr);
    const double d = std::strtod(td.c_str(), nullptr);
    double v = std::nan(""), plain = std::nan(""), dx = std::nan("");
    char st = '.';
    try {
      v = stan::math::digamma_diff(x, d);
      plain = stan::math::digamma(x + d) - stan::math::digamma(x);
      stan::math::var xv = x;
      stan::math::var r = stan::math::digamma_diff(xv, d);
      r.grad();
      dx = xv.adj();
    } catch (const std::exception&) {
      st = 'T';
    }
    stan::math::recover_memory();
    std::printf("%a %a => %a %a %a %c\n", x, d, v, plain, dx, st);
  }
  return 0;
}
