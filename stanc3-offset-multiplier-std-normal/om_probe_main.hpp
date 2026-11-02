// Probe main for a model without data: reads unconstrained parameter values
// from stdin and prints, for (propto, jacobian) = (1,1), (1,0), (0,1), (0,0),
// one line "propto jacobian lp grad...".  Without input it prints the
// number of unconstrained parameters.  Link it instead of CmdStan's main.o:
// make CMDSTAN_MAIN_O=<this .o> <model>
#include <stan/model/model_base.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/math/rev.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

stan::model::model_base& new_model(stan::io::var_context& data_context,
                                   unsigned int seed, std::ostream* msg_stream);

int main() {
  stan::io::empty_var_context ctx;
  stan::model::model_base& model = new_model(ctx, 0, &std::cerr);
  std::vector<double> u;
  double v;
  while (std::cin >> v) {
    u.push_back(v);
  }
  if (u.empty()) {
    std::cout << model.num_params_r() << "\n";
    return 0;
  }
  if (u.size() != model.num_params_r()) {
    std::cerr << "expected " << model.num_params_r() << " values\n";
    return 1;
  }
  std::cout << std::setprecision(17);
  for (int k = 0; k < 4; ++k) {
    bool propto = k < 2;
    bool jacobian = k % 2 == 0;
    Eigen::Matrix<stan::math::var, -1, 1> x(u.size());
    for (size_t i = 0; i < u.size(); ++i) {
      x(i) = u[i];
    }
    stan::math::var lp;
    if (propto && jacobian) {
      lp = model.log_prob_propto_jacobian(x, &std::cerr);
    } else if (propto) {
      lp = model.log_prob_propto(x, &std::cerr);
    } else if (jacobian) {
      lp = model.log_prob_jacobian(x, &std::cerr);
    } else {
      lp = model.log_prob(x, &std::cerr);
    }
    lp.grad();
    std::cout << propto << " " << jacobian << " " << lp.val();
    for (int i = 0; i < x.size(); ++i) {
      std::cout << " " << x(i).adj();
    }
    std::cout << "\n";
    stan::math::recover_memory();
  }
  return 0;
}
