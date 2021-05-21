#pragma once

#include <algorithm>
#include <array>
#include <boost/math/quadrature/gauss.hpp>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>
#include "Traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief Family of tensor product Gauss-Legedre quadrature rules on the
  reference simplex implemented as stateless singletons.

  A quadrature rule in the problem is any class with 4 static methods (two of
  which are constexpr):
  weights() get the weights.
  nodes() get the nodes.
  constexpr measure() get the measure of the reference simplex.
  constexpr size() get the number of quadrature nodes and weights.

  To use this class set the template parameters like this:
  DIM = the dimension of the quadrature points (1, 2 or 3)
  Ns = the number of quadrature nodes per dimension (in 3D set it to p + 2 per
  dimension to have exact integrals if using a method of order p).

  All information will be compute at compile time. A quadrature rule can then be
  passed to a LinearForm and BilinearForm int2d and int3d methods as a template
  parameter.

  The choice to have quadrature decided at compile time may be debated. It is
  less flexible but nodes and weights should be in the stack. (maybe implement
  dynamic quadrature in the future)
*/
// using a class instead of a namespace is because there are no templated
// namespaces
/// \tparam Dim The dimension of the quadrature nodes.
/// \tparam Ns The number of quadrature nodes per coordinate direction.
template <int DIM, int... Ns>
class QuadratureGaussLegendre {
 private:
  // total number of quadrature nodes
  static constexpr int SIZE = (Ns * ...);
  // wrap number of points per dimension in a tuple
  static constexpr auto Ns_tuple = std::make_tuple(Ns...);
  /// This class has no constructor, it's meant to be a singleton.
  QuadratureGaussLegendre() = delete;

 public:
  using scalar_t = Traits::scalar_t;
  /// A quadrature point is templated on the dimension, using type traits so
  /// that in 1D it is a scalar.
  using qpoint_t =
      std::conditional<DIM == 1, Traits::scalar_t, Traits::qpoint_t<DIM>>::type;
  using weights_t = std::array<scalar_t, SIZE>;
  using nodes_t = std::array<qpoint_t, SIZE>;

  /// \brief Get the weights array.
  static auto weights() -> const weights_t& { return _data.first; }
  /// \brief Get the node array.
  static auto nodes() -> const nodes_t& { return _data.second; }
  /// \brief Get the measure of the reference simplex.
  static constexpr auto measure() -> scalar_t { return _measure; }
  /// \brief Get the number of quadrature nodes.
  static constexpr auto size() -> int { return SIZE; }
  /// \brief Get the weights at given index
  static auto weight(int index) -> scalar_t { return _data.first[index]; }
  /// \brief Get the node at given index
  static auto node(int index) -> const qpoint_t& { return _data.second[index]; }
  /// \brief Get th degree of exactness (for now only 1D is correct)
  static constexpr auto exactness() -> int {
    if constexpr (DIM == 1) return 2 * std::get<0>(Ns_tuple) - 1;
    // not so trivial in this cases
    if constexpr (DIM == 2) {
      std::cerr << "unsupported\n";
      return -1;
    }
    if constexpr (DIM == 3) {
      std::cerr << "unsupported\n";
      return -1;
    }
  }

 private:
  /// \brief Computes the quadrature nodes and weights from the boost functions.
  // can't be constexpr because boost is not constexpr
  template <int N>
  static auto _compute_scalar_data()
      -> std::pair<std::array<scalar_t, N>, std::array<scalar_t, N>> {
    // compute scalar quadrature nodes and weights and reorder them, boosts only
    // gives half the nodes and weights as they are symmetric
    const auto& _w = boost::math::quadrature::gauss<scalar_t, N>::weights();
    const auto& _q = boost::math::quadrature::gauss<scalar_t, N>::abscissa();

    std::array<scalar_t, N> w{};
    std::array<scalar_t, N> q{};

    // fill the other half of the data
    std::copy(_w.rbegin(), _w.rend(), w.rbegin());
    std::copy(_w.rbegin(), _w.rend(), w.begin());
    std::copy(_q.rbegin(), _q.rend(), q.rbegin());
    std::transform(_q.rbegin(), _q.rend(), q.begin(),
                   [](auto x) { return -x; });

    return std::make_pair(w, q);
  }
  /// \brief Computes quadrature nodes and weights.
  static auto _compute_data() -> std::pair<weights_t, nodes_t> {
    // if false the code is ignored
    if constexpr (DIM == 1) return _compute_data_impl1();
    if constexpr (DIM == 2) return _compute_data_impl2();
    if constexpr (DIM == 3) return _compute_data_impl3();
  }
  /// \brief Computes quadrature nodes and weights in 1D.
  static auto _compute_data_impl1() -> std::pair<weights_t, nodes_t> {
    static_assert(std::tuple_size<decltype(Ns_tuple)>::value == 1);
    constexpr auto NX = std::get<0>(Ns_tuple);
    auto [wx, qx] = _compute_scalar_data<NX>();

    std::pair<weights_t, nodes_t> data;
    for (int i = 0; i < NX; ++i) {
      data.first[i] = wx[i];
      data.second[i] = qx[i];
    }
    return data;
  }
  /// \brief Computes tensor product quadrature nodes and weights in 2D.
  static auto _compute_data_impl2() -> std::pair<weights_t, nodes_t> {
    static_assert(std::tuple_size<decltype(Ns_tuple)>::value == 2);
    constexpr auto NX = std::get<0>(Ns_tuple);
    constexpr auto NY = std::get<1>(Ns_tuple);
    auto [wx, qx] = _compute_scalar_data<NX>();
    auto [wy, qy] = _compute_scalar_data<NY>();

    std::pair<weights_t, nodes_t> data;
    int index = 0;
    for (int i = 0; i < NX; ++i) {
      for (int j = 0; j < NY; ++j) {
        scalar_t x = (1 + qx[i]) / 2;
        scalar_t y = (1 - qx[i]) * (1 + qy[j]) / 4;
        scalar_t weight = (1 - qx[i]) / 8 * wx[i] * wy[j];

        data.first[index] = weight;
        data.second[index] = qpoint_t{x, y};

        ++index;
      }
    }
    return data;
  }
  /// \brief Computes tensor product quadrature nodes and weights in 3D.
  static auto _compute_data_impl3() -> std::pair<weights_t, nodes_t> {
    static_assert(std::tuple_size<decltype(Ns_tuple)>::value == 3);
    constexpr auto NX = std::get<0>(Ns_tuple);
    constexpr auto NY = std::get<1>(Ns_tuple);
    constexpr auto NZ = std::get<2>(Ns_tuple);
    auto [wx, qx] = _compute_scalar_data<NX>();
    auto [wy, qy] = _compute_scalar_data<NY>();
    auto [wz, qz] = _compute_scalar_data<NZ>();

    std::pair<weights_t, nodes_t> data;
    int index = 0;
    for (int i = 0; i < NX; ++i) {
      for (int j = 0; j < NY; ++j) {
        for (int k = 0; k < NZ; ++k) {
          scalar_t x = (1 + qx[i]) * (1 + qy[j]) * (1 + qz[k]) / 8.0;
          scalar_t y = (1 + qx[i]) * (1 + qy[j]) * (1 - qz[k]) / 8.0;
          scalar_t z = (1 + qx[i]) * (1 - qy[j]) / 4.0;
          scalar_t weight = (1 + qx[i]) * (1 + qx[i]) * (1 + qy[j]) / 64.0 *
                            wx[i] * wy[j] * wz[k];

          data.first[index] = weight;
          data.second[index] = qpoint_t{x, y, z};

          ++index;
        }
      }
    }
    return data;
  }

  // members
  /// Weights and nodes are stored as a pair of arrays.
  static inline std::pair<weights_t, nodes_t> _data = _compute_data();
  /// Measure of the reference simplex.
  static constexpr scalar_t _measure =
      1 / static_cast<scalar_t>(utl::factorial(DIM));
};
}  // namespace wave
