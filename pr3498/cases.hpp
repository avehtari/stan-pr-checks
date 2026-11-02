#ifndef NORMAL_LPDF_SPEED_WORK_BRANCH_CASES_HPP
#define NORMAL_LPDF_SPEED_WORK_BRANCH_CASES_HPP

// Argument kinds and type cases for the develop / branch comparison of
// normal_lpdf.  stan::math::normal_lpdf is the branch version (from the
// tree on the include path); stan::math::normal_lpdf_dev is the develop
// version (dev_normal_lpdf.hpp, generated from origin/develop).

#include <stan/math/mix.hpp>
#include "dev_normal_lpdf.hpp"
#include <string>
#include <vector>

namespace nc {

using Eigen::VectorXd;
using stan::math::var;
using stan::math::var_value;

// DV data VectorXd, VV Matrix<AD, -1, 1>, VM var_value<VectorXd>,
// VS AD scalar, DS double, SDV std::vector<double>, SVV std::vector<AD>,
// RV Matrix<AD, 1, -1>, RM var_value<RowVectorXd>
enum kind { DV, VV, VM, VS, DS, SDV, SVV, RV, RM };

constexpr bool is_scalar_kind(int k) { return k == VS || k == DS; }
constexpr bool is_ad_kind(int k) { return k != DV && k != DS && k != SDV; }
constexpr bool is_varmat_kind(int k) { return k == VM || k == RM; }

struct type_case {
  const char* name;
  int y, mu, sigma;
};

// Cases A-H of the plan, then extra combinations.
constexpr type_case type_cases[] = {
    {"A", DV, VM, VS},  {"B", DV, VV, VS},  {"C", DV, VV, VV},
    {"D", VV, VS, VS},  {"E", DV, VV, DS},  {"F", DV, DV, VS},
    {"G", SDV, SVV, VS}, {"H", DV, DV, DS}, {"I", VV, VV, VV},
    {"J", DV, DV, VV},  {"K", VS, VS, VV},  {"L", VM, VM, VM},
    {"M", DV, VS, VS},  {"N", VV, DV, DV},  {"O", VS, VS, VS},
    {"P", VV, VV, VS},  {"Q", DV, VM, VM},  {"R", RV, RV, VS},
    {"S", DV, RM, VS},  {"T", DS, VV, DV},  {"U", SVV, SDV, SVV},
    {"V", DV, VS, VV},  {"W", VV, VS, DV},  {"X", VS, DV, DS},
    {"Y", DV, VV, DV},
};

// Make an argument of kind K from the values v (scalars use v[0]).  For
// forward mode, `tangent` sets d_ = 1 on every element.
template <int K, typename AD>
auto make(const VectorXd& v, bool tangent = false) {
  auto ad = [&](double x) {
    AD a(x);
    if constexpr (!std::is_same_v<AD, var>) {
      if (tangent) {
        a.d_ = 1.0;
      }
    }
    return a;
  };
  if constexpr (K == DV) {
    return VectorXd(v);
  } else if constexpr (K == VV) {
    Eigen::Matrix<AD, -1, 1> m(v.size());
    for (Eigen::Index i = 0; i < v.size(); ++i) {
      m[i] = ad(v[i]);
    }
    return m;
  } else if constexpr (K == RV) {
    Eigen::Matrix<AD, 1, -1> m(v.size());
    for (Eigen::Index i = 0; i < v.size(); ++i) {
      m[i] = ad(v[i]);
    }
    return m;
  } else if constexpr (K == VM) {
    return var_value<VectorXd>(v);
  } else if constexpr (K == RM) {
    return var_value<Eigen::RowVectorXd>(v.transpose());
  } else if constexpr (K == VS) {
    return ad(v.size() > 0 ? v[0] : 1.0);
  } else if constexpr (K == DS) {
    return v.size() > 0 ? v[0] : 1.0;
  } else if constexpr (K == SDV) {
    return std::vector<double>(v.data(), v.data() + v.size());
  } else {
    std::vector<AD> s;
    for (Eigen::Index i = 0; i < v.size(); ++i) {
      s.push_back(ad(v[i]));
    }
    return s;
  }
}

// Reverse-mode adjoints of an argument as a VectorXd (empty for data).
template <typename T>
VectorXd adjoints(const T& x) {
  if constexpr (std::is_same_v<T, var>) {
    return VectorXd::Constant(1, x.adj());
  } else if constexpr (stan::is_var_matrix<T>::value) {
    return VectorXd(Eigen::Map<const VectorXd>(x.adj().data(),
                                               x.adj().size()));
  } else if constexpr (stan::is_eigen<T>::value
                       && std::is_same_v<stan::value_type_t<T>, var>) {
    VectorXd a(x.size());
    for (Eigen::Index i = 0; i < x.size(); ++i) {
      a[i] = x.coeff(i).adj();
    }
    return a;
  } else if constexpr (stan::is_std_vector<T>::value
                       && std::is_same_v<stan::value_type_t<T>, var>) {
    VectorXd a(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      a[i] = x[i].adj();
    }
    return a;
  } else {
    return VectorXd(0);
  }
}

}  // namespace nc
#endif
