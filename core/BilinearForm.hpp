#pragma once

#include <Eigen/SparseCore>
#include <algorithm>
#include <array>
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
#include "mesh_include.hpp"

namespace wave {
/**
  \brief This class represents a BilinearForm as in the variational formulation.

  A BilinearForm is built on top of a FeSpace, it can be given an expression as
  imput to the int2d and int3d methods as a mathematical integral expression
  written on paper (these two methods takes as template parameter the quadrature
  rule as in QuadratureGaussLegendre or other rules). Once the expression is
  given, the assembly happens and all the matrix coefficients are saved in
  triplets. The method finalize will build the sparse matrix. If the symmetry
  flag is set to true for every time the int3d or int2d methods are called the
  resulting matrix will be upper triangular, as only half the values will be
  computed. The method matrix will return the sparse underlying matrix to be
  used for further computations.
*/
class BilinearForm {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;
  using triplet_t = Traits::triplet_t;
  using sparse_t = Traits::sparse_t;

  /// \brief Construct a bilinear form associated to a FeSpace
  explicit BilinearForm(const FeSpace &fespace);

  /// \brief Get the number of degrees of freedom (number of rows) in the form.
  int size() const;
  /// \brief Get the underlying sparse matrix associated to the form.
  const sparse_t &matrix() const;
  /// \brief Tests if the form has been built with the symmetry flags.
  bool is_symmetric() const;
  /// \brief Turns the triplets stored into a sparse matrix.
  void finalize();

  /**
    \brief Assembly the triplets related to volume integral terms in the
    formulation.

    \tparam Quadrature The 3D quadrature rule to be used to compute integral (as
    in QuadratureGaussLegendre or QuadratureKeastTetrahedron45). \tparam Op This
    parameter is deduced from the expression argument. \tparam Args These
    parameters are deduced from the expression argument. \param expr An
    Expression built through the expression system (first construct the
    terminals: Trial, Test, Normal, Penalty, Field, then combine them with
    operations such as: grad(), dot(,), avg(), jump(), *, +, -). It must be
    constructed in the argument of the method to work. \param is_symmetric Flag
    to set if the bilinear form is symmetric, this will only compute the upper
    triangular part of the form matrix.
  */
  template <typename Quadrature, typename Op, typename... Args>
  void int3d(const Expression<Op, Args...> &expr, bool is_symmetric = false);
  /**
    \brief Assembly the triplets related to surface integral terms in the
    formulation.

    \tparam Quadrature The 2D quadrature rule to be used to compute integral (as
    in QuadratureGaussLegendre or QuadratureDunavantTriangle25). \tparam Op This
    parameter is deduced from the expression argument. \tparam Args These
    parameters are deduced from the expression argument. \param expr An
    Expression built through the expression system (first construct the
    terminals: Trial, Test, Normal, Penalty, Field, then combine them with
    operations such as: grad(), dot(,), avg(), jump(), *, +, -). It must be
    constructed in the argument of the method to work. \param boundaries These
    are the labels of the element boundaries that are used to control what faces
    are counted in the algorithm. If {..., 0, ...} is included the interior
    faces will be considered, if {..., 1, 2, 3, ...} are included then boundary
    number 1, 2 and 3 exterior faces will be also considered. \param
    is_symmetric Flag to set if the bilinear form is symmetric, this will only
    compute the upper triangular part of the form matrix.
  */
  template <typename Quadrature, typename Op, typename... Args>
  void int2d(const Expression<Op, Args...> &expr,
             const std::set<int> &boundaries, bool is_symmetric = false);

 private:
  /// The FeSpace on which the BilinearForm is built.
  const FeSpace &_fespace;
  /// The associated sparse matrix to the form.
  sparse_t _matrix;
  /// The vector of triplets used to build the matrix.
  std::vector<std::pair<bool, std::vector<triplet_t> > > _triplets;
};

//------------------------------------------------------------------------------
inline BilinearForm::BilinearForm(const FeSpace &fespace) : _fespace(fespace) {
  _matrix.resize(fespace.ndof(), fespace.ndof());
};

inline int BilinearForm::size() const { return _fespace.ndof(); }

inline const BilinearForm::sparse_t &BilinearForm::matrix() const {
  return _matrix;
}

inline bool BilinearForm::is_symmetric() const {
  // the form is symmetric if all the integrations were set to be symmetric
  for (const auto &ts : _triplets) {
    if (ts.first == false) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// if symmetric only the upper triangular part is computed
template <typename Quadrature, typename Op, typename... Args>
void BilinearForm::int3d(const Expression<Op, Args...> &expr,
                         bool is_symmetric) {
#ifdef VERBOSE
  utl::Timer timer("bilinear int3d");
#endif

  // prefetch mesh
  const auto &mesh = _fespace.mesh();

  // prefetch quadrature information (because of templates nodes and weights
  // are in the stack, as they are known at compile time)
  constexpr auto QSIZE = Quadrature::size();
  constexpr auto reference_measure = Quadrature::measure();
  const auto weights = Quadrature::weights();
  const auto nodes = Quadrature::nodes();
  // helper to treat vectors as matrices (like sub2ind)
  auto at = [&QSIZE](auto q, auto i) { return q + i * QSIZE; };

  // estimate number of triplets
  _triplets.emplace_back(is_symmetric, std::vector<triplet_t>{});
  int ntriplets = 0;
  auto compute_ntriplets =
      (is_symmetric ? [](int n) { return (n * (n + 1)) / 2; }
                    : [](int n) { return n * n; });

  for (ElementIterator e(mesh); !e.end(); ++e) {
    const auto &fe = _fespace(*e);
    ntriplets += fe.triangulation().size() * compute_ntriplets(fe.ndof());
  }
  _triplets.back().second.reserve(ntriplets);

  // loop over the elements
  for (ElementIterator e(mesh); !e.end(); ++e) {
    // get the associated finite element
    const auto &fe = _fespace(*e);

    // loop over the triangulation of an element
    for (const auto &[jacobian, translation, simplex_measure] :
         fe.triangulation()) {
      // precompute measure ratio
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

      // precompute basis functions on nodes (no initializaion here to avoid
      // wasting time, seg-fault never happens as if one of these is not used in
      // the expression it will not be accessed)
      std::vector<scalar_t> phi;
      std::vector<point_t> dphi;

      // parse the expression at compile time and compute only what is needed
      if constexpr (has_Trial_or_Test_v<decltype(expr)>) {
        phi = fe.eval_basis(x);
      }
      if constexpr (has_GradTrial_or_GradTest_v<decltype(expr)>) {
        dphi = fe.eval_basis_derivatives(x);
      }

      // compute matrix entries
      // loop over dof
      for (int i = 0; i < fe.ndof(); ++i) {
        for (int j = (is_symmetric ? i : 0); j < fe.ndof(); ++j) {
          scalar_t value = 0.0;

          // loop over quadrature nodes in the tetrahedron
          for (int q = 0; q < QSIZE; ++q) {
            value += expr.eval(phi[at(q, j)], phi[at(q, i)], dphi[at(q, j)],
                               dphi[at(q, i)], xf[q]) *
                     weights[q];
          }

          _triplets.back().second.emplace_back(fe.dof(i), fe.dof(j),
                                               value * measure_ratio);
        }
      }
    }
  }
  _triplets.back().second.shrink_to_fit();
}

//------------------------------------------------------------------------------
// if symmetric only the upper triangular part is computed
template <typename Quadrature, typename Op, typename... Args>
void BilinearForm::int2d(const Expression<Op, Args...> &expr,
                         const std::set<int> &boundaries, bool is_symmetric) {
#ifdef VERBOSE
  utl::Timer timer("bilinear int2d");
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

  // estimate number of triplets
  _triplets.emplace_back(is_symmetric, std::vector<triplet_t>{});
  int ntriplets = 0;
  auto compute_ntriplets =
      (is_symmetric ? [](int n) { return (n * (n + 1)) / 2; }
                    : [](int n) { return n * n; });
  for (FaceIterator f(mesh); !f.end(); ++f) {
    ElementIterator e(*f);
    const auto &fef = _fespace(*f);
    const auto &fe = _fespace(*e);
    ntriplets += fef.triangulation().size() * 2 * compute_ntriplets(fe.ndof());
  }
  _triplets.back().second.reserve(ntriplets);

  // loop over faces
  for (FaceIterator f(mesh); !f.end(); ++f) {
    // get the associated finite element face and finite element that shares the
    // face
    const auto &fef = _fespace(*f);
    ElementIterator e(*f);
    const auto &fe = _fespace(*e);

    // here only faces that appear in the set of boundaries will be
    // considered, 0 will mark interior faces.
    if (boundaries.contains(fef.boundary())) {
      if (fef.is_interior()) {
        // get the second element that shares the face (if interior)
        ElementIterator en(*f, 1);
        const auto &fen = _fespace(*en);

        // loop over triangulation
        for (const auto &[jacobian, translation, simplex_measure] :
             fef.triangulation()) {
          // precompute ratio
          scalar_t measure_ratio = simplex_measure / reference_measure;
          // precompute transformed nodes
          std::array<point_t, QSIZE> x;
          std::array<point_t, QSIZE> xf;
          std::array<point_t, QSIZE> xn;

          int index = 0;
          for (const auto &node : nodes) {
            xf[index] = jacobian * node + translation;
            x[index] = fe.map(xf[index]);
            xn[index] = fen.map(xf[index]);
            ++index;
          }

          // precompute basis functions on nodes (no initializaion here to avoid
          // wasting time, seg-fault never happens as if one of these is not
          // used in the expression it will not be accessed)
          std::vector<scalar_t> phi;
          std::vector<point_t> dphi;
          std::vector<scalar_t> phin;
          std::vector<point_t> dphin;

          // parse the expression at compile time and compute only what's needed
          if constexpr (has_Trial_or_Test_v<decltype(expr)>) {
            phi = fe.eval_basis(x);
            phin = fen.eval_basis(xn);
          }
          if constexpr (has_GradTrial_or_GradTest_v<decltype(expr)>) {
            dphi = fe.eval_basis_derivatives(x);
            dphin = fen.eval_basis_derivatives(xn);
          }

          // compute same element terms for current element
          for (int i = 0; i < fe.ndof(); ++i) {
            for (int j = (is_symmetric ? i : 0); j < fe.ndof(); ++j) {
              scalar_t value = 0.0;

              // loop over quadrature nodes
              for (int q = 0; q < QSIZE; ++q) {
                value += expr.eval(fef, FeFace::Side::POSITIVE,
                                   FeFace::Side::POSITIVE, phi[at(q, j)],
                                   phi[at(q, i)], dphi[at(q, j)],
                                   dphi[at(q, i)], xf[q]) *
                         weights[q];
              }
              _triplets.back().second.emplace_back(
                  fe.dof(i), fe.dof(j),
                  value * simplex_measure / reference_measure);
            }
          }
          // compute same element terms for neighboring element
          for (int i = 0; i < fen.ndof(); ++i) {
            for (int j = (is_symmetric ? i : 0); j < fen.ndof(); ++j) {
              scalar_t value = 0.0;

              for (int q = 0; q < QSIZE; ++q) {
                value += expr.eval(fef, FeFace::Side::NEGATIVE,
                                   FeFace::Side::NEGATIVE, phin[at(q, j)],
                                   phin[at(q, i)], dphin[at(q, j)],
                                   dphin[at(q, i)], xf[q]) *
                         weights[q];
              }
              _triplets.back().second.emplace_back(fen.dof(i), fen.dof(j),
                                                   value * measure_ratio);
            }
          }

          // compute neighbor terms
          // if the computation is symmetric only half the entries need
          // to be saved, in particular if the degree of freedom of the
          // current positive element are bigger (just by comparing the
          // smallest in this case as they are incremental) then those
          // of the negative neighboring element, the block is going to
          // be in the upper triangular part of the full matrix, and so
          // it is computed, otherwise the neighboring one is computed
          bool is_current_block_upper_triangular = fe.dof()[0] >= fen.dof()[0];
          bool do_compute_current_block =
              !is_symmetric || is_current_block_upper_triangular;
          bool do_compute_neighboring_block =
              !is_symmetric || !is_current_block_upper_triangular;

          if (do_compute_current_block) {
            for (int i = 0; i < fen.ndof(); ++i) {
              for (int j = 0; j < fe.ndof(); ++j) {
                scalar_t value = 0.0;
                for (int q = 0; q < QSIZE; ++q) {
                  value += expr.eval(fef, FeFace::Side::POSITIVE,
                                     FeFace::Side::NEGATIVE, phi[at(q, j)],
                                     phin[at(q, i)], dphi[at(q, j)],
                                     dphin[at(q, i)], xf[q]) *
                           weights[q];
                }
                _triplets.back().second.emplace_back(fen.dof(i), fe.dof(j),
                                                     value * measure_ratio);
              }
            }
          }
          if (do_compute_neighboring_block) {
            for (int i = 0; i < fe.ndof(); ++i) {
              for (int j = 0; j < fen.ndof(); ++j) {
                scalar_t value = 0.0;
                for (int q = 0; q < QSIZE; ++q) {
                  value += expr.eval(fef, FeFace::Side::NEGATIVE,
                                     FeFace::Side::POSITIVE, phin[at(q, j)],
                                     phi[at(q, i)], dphin[at(q, j)],
                                     dphi[at(q, i)], xf[q]) *
                           weights[q];
                }
                _triplets.back().second.emplace_back(fe.dof(i), fen.dof(j),
                                                     value * measure_ratio);
              }
            }
          }
        }  // end loop triangulation
      }    // end if face interior
      else {
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

          // precompute basis functions on nodes (no initializaion here to avoid
          // wasting time, seg-fault never happens as if one of these is not
          // used in the expression it will not be accessed)
          std::vector<scalar_t> phi;
          std::vector<point_t> dphi;

          if constexpr (has_Trial_or_Test_v<decltype(expr)>) {
            phi = fe.eval_basis(x);
          }
          if constexpr (has_GradTrial_or_GradTest_v<decltype(expr)>) {
            dphi = fe.eval_basis_derivatives(x);
          }

          // compute same element terms
          for (int i = 0; i < fe.ndof(); ++i) {
            for (int j = (is_symmetric ? i : 0); j < fe.ndof(); ++j) {
              scalar_t value = 0.0;
              for (int q = 0; q < QSIZE; ++q) {
                value += expr.eval(fef, FeFace::Side::POSITIVE,
                                   FeFace::Side::POSITIVE, phi[at(q, j)],
                                   phi[at(q, i)], dphi[at(q, j)],
                                   dphi[at(q, i)], xf[q]) *
                         weights[q];
              }
              _triplets.back().second.emplace_back(fe.dof(i), fe.dof(j),
                                                   value * measure_ratio);
            }
          }
        }  // end loop triangulation
      }    // end else face exterior
    }      // end if face is interior or boundary
  }        // end loop over faces
  _triplets.back().second.shrink_to_fit();
}

}  // namespace wave
