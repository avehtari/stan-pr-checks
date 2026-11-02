// Probe of stan::io::random_var_context for an array of tuples.
// write_array stores at element by element (at[1].1, at[1].2, at[2].1,
// at[2].2).  The var_context layout is slot by slot ("at.1" holds at[1].1
// and at[2].1).  random_var_context cuts the write_array output by
// get_dims, so vals_r("at.1") gets at[1].1 and at[1].2.  A partial init
// then reads a value <= 0 into the slot with lower=10.
#include <stan/services/util/initialize.hpp>
#include <stan/services/util/create_rng.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/array_var_context.hpp>
#include <stan/io/random_var_context.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include "rvc_tuples.hpp"
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static void print_vec(const char* label, const std::vector<double>& v) {
  std::printf("%s", label);
  for (double x : v) {
    std::printf(" %.6f", x);
  }
  std::printf("\n");
}

static size_t count(const std::string& s, const std::string& pat) {
  size_t n = 0;
  for (size_t pos = s.find(pat); pos != std::string::npos;
       pos = s.find(pat, pos + 1)) {
    ++n;
  }
  return n;
}

// Runs initialize and prints the outcome and the number of rejections.
static void run_init(stan_model& model, bool partial, double radius,
                     unsigned int seed) {
  std::vector<std::string> names{"a"};
  std::vector<double> vals{0.5};
  std::vector<std::vector<size_t>> dims{{}};
  stan::io::array_var_context partial_context(names, vals, dims);
  stan::io::empty_var_context empty_context;
  std::stringstream out, warn;
  stan::callbacks::stream_logger logger(out, out, warn, warn, warn);
  std::stringstream init_ss;
  stan::callbacks::stream_writer init_writer(init_ss);
  stan::rng_t rng = stan::services::util::create_rng(seed, 1);
  std::string result = "OK";
  try {
    if (partial) {
      stan::services::util::initialize(model, partial_context, rng, radius,
                                       false, logger, init_writer);
    } else {
      stan::services::util::initialize(model, empty_context, rng, radius,
                                       false, logger, init_writer);
    }
  } catch (const std::exception& e) {
    result = "FAILED";
  }
  std::printf("RESULT init partial=%d radius=%g seed=%u result=%s "
              "rejections=%zu\n",
              partial ? 1 : 0, radius, seed, result.c_str(),
              count(warn.str(), "Rejecting initial value"));
}

int main() {
  stan::io::empty_var_context empty_context;
  std::stringstream model_ss;
  stan_model model(empty_context, 0, &model_ss);

  stan::rng_t rng = stan::services::util::create_rng(1, 1);
  stan::io::random_var_context random_context(model, rng, 2, false);
  std::vector<double> u = random_context.get_unconstrained();
  std::vector<int> params_i;
  std::vector<double> constrained;
  model.write_array(rng, u, params_i, constrained, false, false);
  std::vector<std::string> flat_names;
  model.constrained_param_names(flat_names, false, false);
  std::printf("write_array order:");
  for (size_t n = 0; n < flat_names.size(); ++n) {
    std::printf(" %s=%.6f", flat_names[n].c_str(), constrained[n]);
  }
  std::printf("\n");
  print_vec("random_var_context vals_r(\"at.1\"):",
            random_context.vals_r("at.1"));
  print_vec("random_var_context vals_r(\"at.2\"):",
            random_context.vals_r("at.2"));

  run_init(model, false, 2, 1);
  run_init(model, true, 0, 1);
  for (unsigned int seed = 1; seed <= 3; ++seed) {
    run_init(model, true, 2, seed);
  }
  return 0;
}
