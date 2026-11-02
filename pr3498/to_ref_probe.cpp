// Does develop's to_ref(as_value_column_array_or_scalar(x)) copy the values?
// to_ref returns ref_type_t<T&&>.  That type is a plain array when
// Eigen::Ref cannot match the expression, and the return statement then
// evaluates the expression into a new heap array.
// Build (M: stan-dev/math clone, T: its built lib/tbb):
//   g++ -std=c++17 -O2 -D_REENTRANT -DBOOST_DISABLE_ASSERTS -I$M \
//     -I$M/lib/eigen_5.0.1 -I$M/lib/boost_1.87.0 \
//     -I$M/lib/sundials_6.1.1/include -I$M/lib/tbb_2020.3/include \
//     to_ref_probe.cpp -o to_ref_probe -L$T -ltbb -Wl,-rpath,$T
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cstdio>
#include <cstdlib>
#include <type_traits>

static long big_allocs = 0;
static bool counting = false;
extern "C" void* __libc_malloc(size_t);
extern "C" void* malloc(size_t n) {
  if (counting && n >= 8000) {
    ++big_allocs;
  }
  return __libc_malloc(n);
}

template <typename T>
void report(const char* name, const T& x) {
  using T_ret = decltype(stan::math::to_ref(
      stan::math::as_value_column_array_or_scalar(x)));
  big_allocs = 0;
  counting = true;
  double s;
  {
    decltype(auto) v
        = stan::math::to_ref(stan::math::as_value_column_array_or_scalar(x));
    s = v.sum();
  }
  counting = false;
  std::printf("%-12s to_ref returns a plain type: %d   allocations >= 8000 "
              "bytes: %ld   (sum %g)\n",
              name, static_cast<int>(stan::is_plain_type<T_ret>::value),
              big_allocs, s);
}

int main() {
  const int N = 1000;
  Eigen::VectorXd xd = Eigen::VectorXd::Random(N);
  Eigen::Matrix<stan::math::var, -1, 1> xv = xd;
  report("Matrix<var>", xv);
  report("VectorXd", xd);
  return 0;
}
