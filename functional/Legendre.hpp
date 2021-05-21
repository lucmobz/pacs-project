#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include "Traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief This namespace groups the functions used to compute the tensor product
  3D Legendre polynomial basis.

  The Legendre basis, derivatives and normalization constants can  be computed
  up to a composite polynomial degree equal to 10.
*/
namespace legendre {

using scalar_t = Traits::scalar_t;
using point_t = Traits::point_t;

/// \brief This namespace contains implemention details for the Legendre basis.
namespace detail {
// dimensionality of the problem
inline constexpr int DIM = Traits::DIM;
// maximal supported polynomial degree
inline constexpr int DEGREE = 10;
// number of degrees of freedom for the local polynomial space
inline constexpr int NDOF =
    utl::factorial(DIM + DEGREE) / utl::factorial(DEGREE) / utl::factorial(DIM);

/// \brief Evaluate 1D Legendre polynomials with the Horner's rule.
/// \param i Degree of the polynomial.
/// \param x 1D point to evaluate.
constexpr auto _phi(int i, scalar_t x) -> scalar_t {
  if (i == 0)
    return 1;
  else if (i == 1)
    return x;
  else if (i == 2)
    return 1.5 * x * x - 0.5;
  else if (i == 3)
    return (2.5 * x * x - 1.5) * x;
  else if (i == 4)
    return ((35 * x * x - 30) * x * x + 3) / 8.0;
  else if (i == 5)
    return (((63 * x * x - 70) * x * x + 15) * x) / 8.0;
  else if (i == 6)
    return (((231 * x * x - 315) * x * x + 105) * x * x - 5) / 16.0;
  else if (i == 7)
    return ((((429 * x * x - 693) * x * x + 315) * x * x - 35) * x) / 16.0;
  else if (i == 8)
    return ((((6435 * x * x - 12012) * x * x + 6930) * x * x - 1260) * x * x +
            35) /
           128.0;
  else if (i == 9)
    return (((((12155 * x * x - 25740) * x * x + 18018) * x * x - 4620) * x *
                 x +
             315) *
            x) /
           128.0;
  else if (i == 10)
    return (((((46189 * x * x - 109395) * x * x + 90090) * x * x - 30030) * x *
                 x +
             3465) *
                x * x -
            63) /
           256.0;
  else
    // never use exceptions they slow down the code
    return 0.0;
}

/// \brief Evaluate 1D Legendre polynomial derivatives with the Horner's rule.
/// \param i Degree of the polynomial.
/// \param x 1D point to evaluate.s
constexpr auto _dphi(int i, scalar_t x) -> scalar_t {
  if (i == 0)
    return 0.0;
  else if (i == 1)
    return 1.0;
  else if (i == 2)
    return 3.0 * x;
  else if (i == 3)
    return 7.5 * x * x - 1.5;
  else if (i == 4)
    return ((35 * 4 * x * x - 30 * 2) * x) / 8.0;
  else if (i == 5)
    return ((63 * 5 * x * x - 70 * 3) * x * x + 15) / 8.0;
  else if (i == 6)
    return (((231 * 6 * x * x - 315 * 4) * x * x + 105 * 2) * x) / 16.0;
  else if (i == 7)
    return (((429 * 7 * x * x - 693 * 5) * x * x + 315 * 3) * x * x - 35) /
           16.0;
  else if (i == 8)
    return ((((6435 * 8 * x * x - 12012 * 6) * x * x + 6930 * 4) * x * x -
             1260 * 2) *
            x) /
           128.0;
  else if (i == 9)
    return ((((12155 * 9 * x * x - 25740 * 7) * x * x + 18018 * 5) * x * x -
             4620 * 3) *
                x * x +
            315) /
           128.0;
  else if (i == 10)
    return (((((46189 * 10 * x * x - 109395) * 8 * x * x + 90090 * 6) * x * x -
              30030 * 4) *
                 x * x +
             3465 * 2) *
            x) /
           256.0;
  else
    // never use exceptions they slow down the code
    return 0.0;
}

/// \brief Compute the deegres of the polynomials in the tensor product Legendre
/// basis and arrange them hierarchically.
constexpr auto _compute_degrees() -> std::array<std::array<int, DIM>, NDOF> {
  std::array<std::array<int, DIM>, NDOF> ans;

  int index = 0;
  for (int i = 0; i <= DEGREE; ++i) {
    for (int j = 0; j <= DEGREE - i; ++j) {
      for (int k = 0; k <= DEGREE - i - j; ++k) {
        ans[index] = std::array<int, DIM>{i, j, k};
        ++index;
      }
    }
  }

  // compile time sort to order the basis functions in a hierarchical increasing
  // order, so that lower order polynomials appear first, this allow for easy
  // looping through a set of basis functions of maximal given polynomial
  // degree.
  std::sort(ans.begin(), ans.end(), [](const auto &lhs, const auto &rhs) {
    auto lhs_sum = lhs[0] + lhs[1] + lhs[2];
    auto rhs_sum = rhs[0] + rhs[1] + rhs[2];
    if (lhs_sum == rhs_sum) {
      return lhs > rhs;
    } else
      return lhs_sum < rhs_sum;
  });

  return ans;
}

/// Store at compile time the degrees of the polynomials appearing in the basis.
inline constexpr std::array<std::array<int, DIM>, NDOF> _degrees =
    _compute_degrees();

/// \brief Compute the normalization constants for the Legendre tensor product
/// basis.
constexpr auto _compute_normalization() -> std::array<scalar_t, NDOF> {
  std::array<scalar_t, NDOF> ans;
  for (int index = 0; index < NDOF; ++index) {
    auto [i, j, k] = _degrees[index];
    ans[index] = std::sqrt((2 * i + 1) * (2 * j + 1) * (2 * k + 1) / 8.0);
  }
  return ans;
}

/// Store at compile time all the normalization constants up to the maximal
/// degree.
inline constexpr std::array<scalar_t, NDOF> _normalization =
    _compute_normalization();
}  // namespace detail

using namespace detail;
/// \brief Returns the normalization constant for the particular Legendre 3D
/// basis function of given index.
inline constexpr auto normalization(int index) -> scalar_t {
  return _normalization[index];
}

/// \brief Returns the Legendre 3D basis function of given index evaluated at x.
inline auto phi(int index, const point_t &x) -> scalar_t {
  auto deg = _degrees[index];
  return _phi(deg[0], x[0]) * _phi(deg[1], x[1]) * _phi(deg[2], x[2]);
}
/// \brief Returns the Legendre 3D basis function x-derivative of given index
/// evaluated at x.
inline auto phix(int index, const point_t &x) -> scalar_t {
  auto deg = _degrees[index];
  return _dphi(deg[0], x[0]) * _phi(deg[1], x[1]) * _phi(deg[2], x[2]);
}
/// \brief Returns the Legendre 3D basis function y-derivative of given index
/// evaluated at x.
inline auto phiy(int index, const point_t &x) -> scalar_t {
  auto deg = _degrees[index];
  return _phi(deg[0], x[0]) * _dphi(deg[1], x[1]) * _phi(deg[2], x[2]);
}
/// \brief Returns the Legendre 3D basis function z-derivative of given index
/// evaluated at x.
inline auto phiz(int index, const point_t &x) -> scalar_t {
  auto deg = _degrees[index];
  return _phi(deg[0], x[0]) * _phi(deg[1], x[1]) * _dphi(deg[2], x[2]);
}
}  // namespace legendre
}  // namespace wave