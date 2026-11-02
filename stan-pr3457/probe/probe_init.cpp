// Probe of stan::services::util::initialize with partial inits, for the
// claims of the PR text (missing entries equal the random draws, supplied
// values come back, runs without partial inits and the RNG stream do not
// change).  Built once with the develop initialize.hpp and once with the
// PR version; compare_outputs.py compares the two outputs.
#include <stan/services/util/initialize.hpp>
#include <stan/services/util/create_rng.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/io/array_var_context.hpp>
#include <stan/io/random_var_context.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include "partial_init_dependent.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct Var {
  std::string name;
  std::vector<double> values;
  std::vector<size_t> dims;
};

static stan::io::array_var_context make_context(const std::vector<Var>& vars) {
  std::vector<std::string> names;
  std::vector<double> values;
  std::vector<std::vector<size_t>> dims;
  for (const auto& v : vars) {
    names.push_back(v.name);
    values.insert(values.end(), v.values.begin(), v.values.end());
    dims.push_back(v.dims);
  }
  return stan::io::array_var_context(names, values, dims);
}

// var_context name of a flat name: keep each ":k" (tuple slot) as ".k" and
// drop each ".i" (index), so "at.1:2" -> "at.2" and "w:2:1" -> "w.2.1".
static std::string slot_name(const std::string& flat) {
  size_t pos = flat.find_first_of(".:");
  std::string slot = flat.substr(0, pos);
  while (pos != std::string::npos) {
    size_t next = flat.find_first_of(".:", pos + 1);
    if (flat[pos] == ':') {
      slot += "." + flat.substr(pos + 1, next - pos - 1);
    }
    pos = next;
  }
  return slot;
}

// One initialization; prints one line per seed.
static void run_case(stan_model& model, const std::string& label,
                     const std::vector<Var>& supplied, double radius,
                     int n_seeds, bool print_params) {
  stan::io::array_var_context init_context = make_context(supplied);
  std::set<std::string> supplied_slots;
  for (const auto& v : supplied) {
    supplied_slots.insert(v.name);
  }
  std::vector<std::string> u_names, c_names;
  model.unconstrained_param_names(u_names, false, false);
  model.constrained_param_names(c_names, false, false);
  for (int seed = 0; seed < n_seeds; ++seed) {
    stan::rng_t rng = stan::services::util::create_rng(seed, 1);
    stan::rng_t rng_copy = rng;
    stan::io::random_var_context random_context(model, rng_copy, radius,
                                                radius == 0);
    const std::vector<double> u_random = random_context.get_unconstrained();
    std::stringstream log_ss, init_ss;
    stan::callbacks::stream_logger logger(log_ss, log_ss, log_ss, log_ss,
                                          log_ss);
    stan::callbacks::stream_writer init_writer(init_ss);
    std::vector<double> params;
    try {
      params = stan::services::util::initialize(model, init_context, rng,
                                                radius, false, logger,
                                                init_writer);
    } catch (const std::exception& e) {
      std::string text = log_ss.str();
      size_t rejections = 0;
      for (size_t pos = text.find("Rejecting initial value");
           pos != std::string::npos;
           pos = text.find("Rejecting initial value", pos + 1)) {
        ++rejections;
      }
      std::printf("%s seed %2d status throw rejections %zu\n", label.c_str(),
                  seed, rejections);
      continue;
    }
    // missing entries against the first random draw
    double max_abs_missing = 0, max_rel_missing = 0;
    for (size_t n = 0; n < params.size(); ++n) {
      if (supplied_slots.count(slot_name(u_names[n])) == 0) {
        double d = std::fabs(params[n] - u_random[n]);
        max_abs_missing = std::max(max_abs_missing, d);
        if (u_random[n] != 0) {
          max_rel_missing = std::max(max_rel_missing, d / std::fabs(u_random[n]));
        }
      }
    }
    // supplied values after the round trip
    std::vector<int> params_i;
    std::vector<double> vars;
    stan::rng_t rng_wa = stan::services::util::create_rng(0, 1);
    model.write_array(rng_wa, params, params_i, vars, false, false);
    double max_rel_supplied = 0;
    // The flat entries of a slot come in var_context order for the
    // variables of this probe (1-d arrays of tuples with scalar slots).
    size_t unmatched = 0;
    for (const auto& v : supplied) {
      std::vector<size_t> idx;
      for (size_t n = 0; n < c_names.size(); ++n) {
        if (slot_name(c_names[n]) == v.name) {
          idx.push_back(n);
        }
      }
      if (idx.size() != v.values.size()) {
        ++unmatched;
        continue;
      }
      for (size_t k = 0; k < v.values.size(); ++k) {
        double d = std::fabs(vars[idx[k]] - v.values[k]);
        if (v.values[k] != 0) {
          d /= std::fabs(v.values[k]);
        }
        max_rel_supplied = std::max(max_rel_supplied, d);
      }
    }
    if (unmatched > 0) {
      std::printf("%s seed %2d WARNING %zu supplied variables not matched\n",
                  label.c_str(), seed, unmatched);
    }
    unsigned long long next = rng();  // the next draw: same stream?
    std::printf(
        "%s seed %2d status ok max_abs_missing %.3g max_rel_missing %.3g "
        "max_rel_supplied %.3g max_abs_unconstrained %.3g rng_next %llu\n",
        label.c_str(), seed, max_abs_missing, max_rel_missing,
        max_rel_supplied,
        std::fabs(*std::max_element(params.begin(), params.end(),
                                    [](double a, double b) {
                                      return std::fabs(a) < std::fabs(b);
                                    })),
        next);
    if (print_params) {
      std::printf("%s seed %2d params", label.c_str(), seed);
      for (double p : params) {
        std::printf(" %a", p);
      }
      std::printf("\n");
    }
  }
}

int main() {
  stan::io::empty_var_context empty;
  std::stringstream model_ss;
  stan_model model(empty, 12345, &model_ss);
  const int n_seeds = 20;

  std::vector<Var> full = {
      {"mu", {0.5}, {}},        {"sigma", {1.5}, {}},
      {"x", {0.2}, {}},         {"s", {0.2, 0.3, 0.5}, {3}},
      {"v", {0.1, -0.2}, {2}},  {"t.1", {0.3}, {}},
      {"t.2", {-0.4}, {}},      {"a", {1}, {}},
      {"b", {1.2}, {}},         {"c", {0.8}, {}},
      {"d", {1.1}, {}},         {"e", {0.9}, {}},
      {"lb", {0}, {}},          {"y", {2}, {}},
      {"at.1", {0.3, -0.2}, {2}}, {"at.2", {0.5, 0.6}, {2}},
      {"z", {0.1}, {}},         {"tv.1", {0.2}, {}},
      {"tv.2", {0.3, -0.3}, {2}}, {"w.1", {0.1}, {}},
      {"w.2.1", {0.2}, {}},     {"w.2.2", {0.3}, {}}};

  // runs that must not change (bitwise; compare_outputs.py)
  run_case(model, "none_r2", {}, 2, n_seeds, true);
  run_case(model, "none_r0", {}, 0, 1, true);
  run_case(model, "full_r2", full, 2, n_seeds, true);
  // partial inits
  run_case(model, "sigma_1e-20", {{"sigma", {1e-20}, {}}}, 2, n_seeds, false);
  run_case(model, "mu_5", {{"mu", {5}, {}}}, 2, n_seeds, false);
  run_case(model, "mu_5_r0", {{"mu", {5}, {}}}, 0, 1, false);
  run_case(model, "x_3", {{"x", {3}, {}}}, 2, n_seeds, false);
  run_case(model, "a_1_c_4", {{"a", {1}, {}}, {"c", {4}, {}}}, 2, n_seeds,
           false);
  run_case(model, "b_2_d_-3", {{"b", {2}, {}}, {"d", {-3}, {}}}, 2, n_seeds,
           false);
  run_case(model, "a_c_e", {{"a", {1}, {}}, {"c", {4}, {}}, {"e", {-2}, {}}},
           2, n_seeds, false);
  run_case(model, "z_2", {{"z", {2}, {}}}, 2, n_seeds, false);
  run_case(model, "t1_sigma", {{"sigma", {1e-20}, {}}, {"t.1", {0.7}, {}}}, 2,
           n_seeds, false);
  run_case(model, "w21_sigma", {{"sigma", {1e-20}, {}}, {"w.2.1", {0.4}, {}}},
           2, n_seeds, false);
  run_case(model, "lb_10", {{"lb", {10}, {}}}, 2, 3, false);
  return 0;
}
