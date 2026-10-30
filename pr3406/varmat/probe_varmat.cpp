// Probe the two helper defects named in the stan-dev/math#3406 description:
// as_array_or_scalar on a var_value that already holds an array, and
// ref_type_t / to_ref of an rvalue arena_matrix.  Build one case per -DCASE_x.
#include <stan/math/rev.hpp>
#include <cstdio>
#include <type_traits>

using stan::math::arena_matrix;
using stan::math::var;
using stan::math::var_value;

// sum(a) with the adjoint pushed back by a callback
template <typename T>
var sum_cb(const T& a) {
  return stan::math::make_callback_var(
      a.val().sum(), [a](auto& vi) mutable { a.adj().array() += vi.adj(); });
}

template <typename X>
void report(const char* name, const X& x, double sum_val) {
  std::printf("%s: sum value %g (expect 6), x.adj = %g %g %g (expect 1 1 1)\n",
              name, sum_val, x.adj()(0), x.adj()(1), x.adj()(2));
}

#ifdef CASE_AS_TWICE
// as_array_or_scalar applied twice to var_value<VectorXd>
int main() {
  Eigen::VectorXd v(3);
  v << 1, 2, 3;
  var_value<Eigen::VectorXd> x(v);
  auto a = stan::math::as_array_or_scalar(x);
  auto b = stan::math::as_array_or_scalar(a);
  std::printf("type of as_array_or_scalar(x): %s\n",
              stan::is_eigen_array<stan::value_type_t<decltype(a)>>::value
                  ? "var_value<array>" : "var_value<matrix>");
  var s = sum_cb(b);
  s.grad();
  report("as_array_or_scalar twice", x, s.val());
  return 0;
}
#endif

#ifdef CASE_AS_ARRAY
// as_array_or_scalar of a var_value<ArrayXd>
int main() {
  Eigen::ArrayXd v(3);
  v << 1, 2, 3;
  var_value<Eigen::ArrayXd> x(v);
  auto b = stan::math::as_array_or_scalar(x);
  var s = sum_cb(b);
  s.grad();
  report("as_array_or_scalar(var_value<ArrayXd>)", x, s.val());
  return 0;
}
#endif

#ifdef CASE_REF_TYPE
// Is ref_type_t of an rvalue arena_matrix a reference?
template <typename T>
void show(const char* name) {
  using R = stan::ref_type_t<T>;
  std::printf("%-40s reference: %s, const: %s\n", name,
              std::is_reference<R>::value ? "yes" : "no",
              std::is_const<std::remove_reference_t<R>>::value ? "yes" : "no");
}
int main() {
  show<arena_matrix<Eigen::VectorXd>&&>("ref_type_t<arena_matrix<VectorXd>&&>");
  show<arena_matrix<Eigen::VectorXd>>("ref_type_t<arena_matrix<VectorXd>>");
  show<arena_matrix<Eigen::VectorXd>&>("ref_type_t<arena_matrix<VectorXd>&>");
  show<arena_matrix<Eigen::Matrix<var, -1, 1>>&&>(
      "ref_type_t<arena_matrix<Matrix<var>>&&>");
  show<Eigen::VectorXd&&>("ref_type_t<VectorXd&&> (control)");
  return 0;
}
#endif

#ifdef CASE_DANGLE
// to_ref of a temporary arena_matrix, bound to const auto&, read later.
// On a tree where ref_type_t<arena_matrix&&> is a reference, the binding
// refers to the dead temporary: ASan reports stack-use-after-scope.
arena_matrix<Eigen::VectorXd> make_arena() {
  Eigen::VectorXd v(3);
  v << 1, 2, 3;
  return arena_matrix<Eigen::VectorXd>(v);
}

int main() {
  const auto& r = stan::math::to_ref(make_arena());
  const auto& r2 = stan::math::to_ref_if<true>(make_arena());
  double s = r.sum() + r2.sum();
  std::printf("to_ref of rvalue arena_matrix: sum %g (expect 12)\n", s);
  stan::math::recover_memory();
  return 0;
}
#endif
