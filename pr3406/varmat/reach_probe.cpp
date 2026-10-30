// Instantiate the functions of stan-dev/math#3406 with var_value<VectorXd>
// and Matrix<var> arguments.  Compiled against trees whose two helpers carry
// a static_assert (reach.sh), it names the callers that reach the defects.
#include <stan/math/rev.hpp>

using stan::math::var;
using VV = stan::math::var_value<Eigen::VectorXd>;
using MV = Eigen::Matrix<var, -1, 1>;

template <typename Y, typename P>
void all(const Y& y, const P& mu, const P& s) {
  using namespace stan::math;
  normal_lcdf(y, mu, s);
  normal_lccdf(y, mu, s);
  normal_cdf(y, mu, s);
  lognormal_lcdf(y, mu, s);
  lognormal_lccdf(y, mu, s);
  lognormal_cdf(y, mu, s);
  exp_mod_normal_lpdf(y, mu, s, s);
  exp_mod_normal_lcdf(y, mu, s, s);
  exp_mod_normal_lccdf(y, mu, s, s);
  exp_mod_normal_cdf(y, mu, s, s);
  skew_normal_lpdf(y, mu, s, s);
}

template <typename Y>
void std_all(const Y& y) {
  using namespace stan::math;
  std_normal_lcdf(y);
  std_normal_lccdf(y);
  std_normal_cdf(y);
  Phi(y);
}

int main() {
  Eigen::VectorXd v = Eigen::VectorXd::Constant(3, 0.5);
  VV yv(v);
  MV ym = v.cast<var>();
  var mu = 0.1, s = 1.2;
  VV muv(v), sv(Eigen::VectorXd::Constant(3, 1.5));
  all(yv, mu, s);
  all(ym, mu, s);
  all(v, muv, sv);
  all(yv, muv, sv);
  std_all(yv);
  std_all(ym);
  return 0;
}
