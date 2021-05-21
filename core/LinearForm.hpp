#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <iostream>
#include <set>
#include <type_traits>
#include <utility>
#include "Expression.hpp"
#include "FeElement.hpp"
#include "FeFace.hpp"
#include "FeSpace.hpp"
#include "Traits.hpp"
#include "expression_type_traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief This class represents a LinearForm as in the variational formulation.

  A LinearForm is built on top of a FeSpace, it can be given an expression as
  imput to the int2d and int3d methods as a mathematical integral expression
  written on paper (these two methods takes as template parameter the quadrature
  rule as in QuadratureGaussLegendre or other rules). Once the expression is
  given, the assembly happens and all the vector coefficients are saved. The
  method vector will return the dense underlying vector to be used for further
  computations.
*/
class LinearForm {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;
  using dense_t = Traits::dense_t;

  /// \brief Construct a LinearForm from an FeSpace
  explicit LinearForm(const FeSpace &fespace);

  /// \brief Get the number of degrees of freedom (the length of the vector)
  int size() const;
  /// \brief Get the form associated vector.
  const dense_t &vector() const;
  /// \brief See BilinearForm int3d method.
  /// The only difference is the absence of the symmetry flag.
  template <typename Quadrature, typename Op, typename... Args>
  void int3d(const Expression<Op, Args...> &expr);
  /// \brief See BilinearForm int2d method.
  /// The only difference is the absence of the symmetry flag.
  template <typename Quadrature, typename Op, typename... Args>
  void int2d(const Expression<Op, Args...> &expr,
             const std::set<int> boundaries);

 private:
  /// The FeSpace on which the form is built.
  const FeSpace &_fespace;
  /// The vector associated with the form.
  dense_t _vector;
};

//------------------------------------------------------------------------------
inline LinearForm::LinearForm(const FeSpace &fespace) : _fespace(fespace) {
  _vector.resize(fespace.ndof());
  _vector.setZero();
};

inline int LinearForm::size() const { return _fespace.ndof(); }

inline const LinearForm::dense_t &LinearForm::vector() const { return _vector; }

//------------------------------------------------------------------------------
// See BilinearForm int3d method for info
template <typename Quadrature, typename Op, typename... Args>
void LinearForm::int3d(const Expression<Op, Args...> &expr) {
#ifdef VERBOSE
  utl::Timer timer("linear int3d");
#endif

  // prefetch mesh
  const auto &mesh = _fespace.mesh();

  // prefetch quadrature information
  constexpr auto QSIZE = Quadrature::size();
  constexpr auto reference_measure = Quadrature::measure();
  const auto weights = Quadrature::weights();
  const auto nodes = Quadrature::nodes();
  // helper to treat vectors as matrices (like sub2ind)
  auto at = [&QSIZE](auto q, auto i) { return q + i * QSIZE; };

  // loop over the elements
  for (ElementIterator e(mesh); !e.end(); ++e) {
    const auto &fe = _fespace(*e);

    for (const auto &[jacobian, translation, simplex_measure] :
         fe.triangulation()) {
      // precompute ratio
      scalar_t measure_ratio = simplex_measure / reference_measure;
      // precompute transformed nodes
      std::array<point_t, QSIZE> x;
      std::array<point_t, QSIZE> xf;

      int index = 0;
      for (const auto &node : nodes) {
        xf[index] = jacobian * node + translation;
        x[index] = fe.map(xf[index]);
        ++index;
      }

      // precompute basis functions on nodes (seg-fault access never
      // happens and this avoids realloc)
      std::vector<scalar_t> phi;
      std::vector<point_t> dphi;

      if constexpr (has_Trial_or_Test_v<decltype(expr)>) {
        phi = fe.eval_basis(x);
      }
      if constexpr (has_GradTrial_or_GradTest_v<decltype(expr)>) {
        dphi = fe.eval_basis_derivatives(x);
      }

      // compute values
      for (int i = 0; i < fe.ndof(); ++i) {
        scalar_t value = 0.0;

        for (int q = 0; q < QSIZE; ++q) {
          value += expr.eval(phi[at(q, i)], dphi[at(q, i)], xf[q]) * weights[q];
        }
        _vector[fe.dof(i)] += value * measure_ratio;
      }
    }
  }
}
//------------------------------------------------------------------------------
// See BilinearForm int2d method for info
template <typename Quadrature, typename Op, typename... Args>
void LinearForm::int2d(const Expression<Op, Args...> &expr,
                       const std::set<int> boundaries) {
#ifdef VERBOSE
  utl::Timer timer("linear int2d");
#endif

  // prefetch mesh
  const auto &mesh = _fespace.mesh();

  // prefetch quadrature information
  constexpr auto QSIZE = Quadrature::size();
  constexpr auto reference_measure = Quadrature::measure();
  const auto weights = Quadrature::weights();
  const auto nodes = Quadrature::nodes();
  // helper to treat vectors as matrices (like sub2ind)
  auto at = [&QSIZE](auto q, auto i) { return q + i * QSIZE; };

  // loop over the faces (signed)
  for (FaceIterator f(mesh); !f.end(); ++f) {
    const auto &fef = _fespace(*f);

    if (fef.is_exterior() && boundaries.contains(fef.boundary())) {
      ElementIterator e(*f);
      const auto &fe = _fespace(*e);

      for (const auto &[jacobian, translation, simplex_measure] :
           fef.triangulation()) {
        // precompute ratio
        scalar_t measure_ratio = simplex_measure / reference_measure;

        // precompute transformed nodes
        std::array<point_t, QSIZE> x;
        std::array<point_t, QSIZE> xf;

        int index = 0;
        for (const auto &node : nodes) {
          xf[index] = jacobian * node + translation;
          x[index] = fe.map(xf[index]);
          ++index;
        }

        // precompute basis functions on nodes (seg-fault access never
        // happens and this avoids realloc)
        std::vector<scalar_t> phi;
        std::vector<point_t> dphi;

        if constexpr (has_Trial_or_Test_v<decltype(expr)>) {
          phi = fe.eval_basis(x);
        }
        if constexpr (has_GradTrial_or_GradTest_v<decltype(expr)>) {
          dphi = fe.eval_basis_derivatives(x);
        }

        for (int i = 0; i < fe.ndof(); ++i) {
          scalar_t value = 0.0;

          for (int q = 0; q < QSIZE; ++q) {
            value += expr.eval(fef, FeFace::Side::POSITIVE, phi[at(q, i)],
                               dphi[at(q, i)], xf[q]) *
                     weights[q];
          }
          _vector[fe.dof(i)] += value * measure_ratio;
        }
      }
    }
  }
}

}  // namespace wave