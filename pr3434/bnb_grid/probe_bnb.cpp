// probe_bnb.cpp --- evaluate beta_neg_binomial_{cdf,lcdf,lccdf} with
// stan::math::var, and hypergeometric_3F2 and grad_F32 in double with the
// arguments that the cdfs pass, at the points of a hex point file. Every
// call runs in a forked child with a time limit.
//
// Usage: probe_bnb <points file> [time limit in seconds, default 30]
//
// Input line:  BNB <label> <hex n> <hex r> <hex alpha> <hex beta>
// Output line: BNB <label> <hex args...> => then five groups
//              <fn> <hex value> <hex g1> <hex g2> <hex g3> <status> <ms>
//   fn C, LC, LCC: cdf, lcdf, lccdf; value and d/d(r, alpha, beta)
//   fn F3: hypergeometric_3F2({1, beta + n + 1, r + n + 1},
//          {n + 2, alpha + r + beta + n + 1}, 1); value, nan, nan, nan
//   fn G3: grad_F32<false, true, true, false, true, false> with the same
//          arguments, precision 1e-8 and max_steps 1e8 (the defaults of
//          the cdfs); nan, dF/da2, dF/da3, dF/db2
//   status: '.' finished, 'T' threw, 'X' time limit, 'C' child died
//   ms: time of the call in the child (the time limit for 'X')

#include <stan/math/rev.hpp>
#ifdef BNB_PROTO
#include "proto_bnb.hpp"
#endif

#include <poll.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Result {
  double v[4];
  char status;
  double ms;
};

void set_nan(double* v) {
  for (int i = 0; i < 4; ++i) {
    v[i] = std::nan("");
  }
}

// ── Run one call in a forked child with a time limit ─────────────────────

template <typename Fn>
Result run_limited(Fn&& fn, int limit_ms) {
  Result res;
  set_nan(res.v);
  res.status = 'C';
  res.ms = std::nan("");
  int fd[2];
  if (pipe(fd) != 0) {
    std::perror("pipe");
    std::exit(3);
  }
  std::fflush(stdout);
  const pid_t pid = fork();
  if (pid < 0) {
    std::perror("fork");
    std::exit(3);
  }
  if (pid == 0) {
    close(fd[0]);
    Result out;
    set_nan(out.v);
    const auto t0 = std::chrono::steady_clock::now();
    try {
      fn(out.v);
      out.status = '.';
    } catch (const std::exception&) {
      set_nan(out.v);
      out.status = 'T';
    }
    const auto t1 = std::chrono::steady_clock::now();
    out.ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const ssize_t w = write(fd[1], &out, sizeof(out));
    _exit(w == static_cast<ssize_t>(sizeof(out)) ? 0 : 1);
  }
  close(fd[1]);
  struct pollfd p = {fd[0], POLLIN, 0};
  const int ready = poll(&p, 1, limit_ms);
  if (ready == 0) {
    kill(pid, SIGKILL);
    waitpid(pid, nullptr, 0);
    close(fd[0]);
    res.status = 'X';
    res.ms = limit_ms;
    return res;
  }
  Result got;
  const ssize_t n_read = read(fd[0], &got, sizeof(got));
  close(fd[0]);
  waitpid(pid, nullptr, 0);
  if (n_read == static_cast<ssize_t>(sizeof(got))) {
    res = got;
  }
  return res;
}

// ── The five calls ────────────────────────────────────────────────────────

#ifdef BNB_PROTO
// Prototype build (-DBNB_PROTO, proto_bnb.hpp): C, LC, LCC come from
// proto::bnb_cdfs in double; F3 holds the route of the upper tail and G3
// the number of tail steps (in the slot of dF/da2).
void eval_cdf(int which, int n, double r, double alpha, double beta,
              double* v) {
  const proto::BnbCdfs res = proto::bnb_cdfs(n, r, alpha, beta);
  if (which == 0) {
    const double cdf = std::exp(res.lcdf);
    v[0] = cdf;
    for (int i = 0; i < 3; ++i) {
      v[i + 1] = cdf * res.dlcdf[i];
    }
  } else {
    v[0] = (which == 1) ? res.lcdf : res.lccdf;
    for (int i = 0; i < 3; ++i) {
      v[i + 1] = (which == 1) ? res.dlcdf[i] : res.dlccdf[i];
    }
  }
}

void eval_f3(int n, double r, double alpha, double beta, double* v) {
  v[0] = proto::bnb_cdfs(n, r, alpha, beta).route;
}

void eval_g3(int n, double r, double alpha, double beta, double* v) {
  v[1] = static_cast<double>(proto::bnb_cdfs(n, r, alpha, beta).steps);
}
#else
void eval_cdf(int which, int n, double r, double alpha, double beta,
              double* v) {
  using stan::math::var;
  var rv = r;
  var av = alpha;
  var bv = beta;
  var lp;
  if (which == 0) {
    lp = stan::math::beta_neg_binomial_cdf(n, rv, av, bv);
  } else if (which == 1) {
    lp = stan::math::beta_neg_binomial_lcdf(n, rv, av, bv);
  } else {
    lp = stan::math::beta_neg_binomial_lccdf(n, rv, av, bv);
  }
  lp.grad();
  v[0] = lp.val();
  v[1] = rv.adj();
  v[2] = av.adj();
  v[3] = bv.adj();
}

#ifdef BNB_NO_PARTS
// -DBNB_NO_PARTS: skip the 3F2 columns (for trees whose cdfs do not use the
// 3F2 series); F3 and G3 then hold nan
void eval_f3(int, double, double, double, double*) {}
void eval_g3(int, double, double, double, double*) {}
#else
// the arguments are formed in the order of the cdf code
void eval_f3(int n, double r, double alpha, double beta, double* v) {
  const double n_dbl = n;
  const double b_plus_n = beta + n_dbl;
  const double r_plus_n = r + n_dbl;
  const double a_plus_r = alpha + r;
  const std::vector<double> a{1.0, b_plus_n + 1.0, r_plus_n + 1.0};
  const std::vector<double> b{n_dbl + 2.0, a_plus_r + b_plus_n + 1.0};
  v[0] = stan::math::hypergeometric_3F2(a, b, 1.0);
}

void eval_g3(int n, double r, double alpha, double beta, double* v) {
  const double n_dbl = n;
  const double b_plus_n = beta + n_dbl;
  const double r_plus_n = r + n_dbl;
  const double a_plus_r = alpha + r;
  double dF[6];
  stan::math::grad_F32<false, true, true, false, true, false>(
      dF, 1.0, b_plus_n + 1.0, r_plus_n + 1.0, n_dbl + 2.0,
      a_plus_r + b_plus_n + 1.0, 1.0, 1e-8, static_cast<int>(1e8));
  v[1] = dF[1];
  v[2] = dF[2];
  v[3] = dF[4];
}
#endif  // BNB_NO_PARTS
#endif  // BNB_PROTO

void print_result(const char* fn, const Result& res) {
  std::printf(" %s", fn);
  for (double x : res.v) {
    std::printf(" %a", x);
  }
  std::printf(" %c %.3f", res.status, res.ms);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argc > 4) {
    std::fprintf(stderr,
                 "usage: probe_bnb <points file> [limit in s] [functions]\n");
    return 2;
  }
  const int limit_ms
      = static_cast<int>(1000 * (argc >= 3 ? std::atof(argv[2]) : 30.0));
  // functions: 1 = cdf only, 2 = cdf and lcdf, 3 (default) = all three
  const int n_fns = argc == 4 ? std::atoi(argv[3]) : 3;
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
    std::string tag, label, token;
    ss >> tag >> label;
    std::vector<double> a;
    while (ss >> token) {
      a.push_back(std::strtod(token.c_str(), nullptr));
    }
    if (tag != "BNB" || a.size() != 4) {
      std::fprintf(stderr, "bad line: %s\n", line.c_str());
      return 2;
    }
    const int n = static_cast<int>(a[0]);
    const double r = a[1];
    const double alpha = a[2];
    const double beta = a[3];
    std::printf("%s %s", tag.c_str(), label.c_str());
    for (double x : a) {
      std::printf(" %a", x);
    }
    std::printf(" =>");
    const char* names[3] = {"C", "LC", "LCC"};
    for (int which = 0; which < n_fns; ++which) {
      print_result(names[which],
                   run_limited(
                       [&](double* v) { eval_cdf(which, n, r, alpha, beta, v); },
                       limit_ms));
    }
    if (n_fns < 3) {
      std::printf("\n");
      std::fflush(stdout);
      continue;
    }
    print_result("F3", run_limited(
                           [&](double* v) { eval_f3(n, r, alpha, beta, v); },
                           limit_ms));
    print_result("G3", run_limited(
                           [&](double* v) { eval_g3(n, r, alpha, beta, v); },
                           limit_ms));
    std::printf("\n");
    std::fflush(stdout);
  }
  return 0;
}
