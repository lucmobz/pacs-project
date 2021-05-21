#pragma once

#include <cmath>
#include <utility>
#include <vector>
#include "FeElement.hpp"
#include "FeFace.hpp"
#include "Traits.hpp"
#include "mesh_include.hpp"

namespace wave {
/**
  \brief A FeSpace represents the finite element space used in the discontinuous
  Galerking method.

  This class stores the finite elements (FeElement's) and finite element faces
  (FeFace's) used in the algorithm that assembles the linear system. It is
  constructed from a mesh and a polynomial degree, or a mesh and a vector of
  polynomial degrees as the degree can vary from element to element.
*/
class FeSpace {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;
  using dense_t = Traits::dense_t;

  /// \brief Construct the FeSpace from a Mesh and a polynomial degree.
  FeSpace(const Mesh &mesh, int degree);
  /// \brief Construct the FeSpace from a Mesh and a vector of polynomial
  /// degrees. Depending on the index of each finite element (FeElement) the
  /// vector sets the degree of such element. This way different degrees can be
  /// used locally.
  FeSpace(const Mesh &mesh, const std::vector<int> &degrees);
  /// \brief Get the Mesh on which the FeSpace is built.
  auto mesh() const -> const Mesh &;
  /// \brief Get the number of degrees of freedom.
  auto ndof() const -> int;
  /// \brief Get the underlying FeElement associated with a mesh Element.
  auto operator()(const Element &element) const -> const FeElement &;
  /// \brief Get the underlying FeFace associated with a mesh Face.
  auto operator()(const Face &face) const -> const FeFace &;
  /**
    \brief From the vector of degrees of freedom of a finite element function,
    get the nodal values on the vertices of the mesh. \param solution The vector
    of degrees of freedom of any finite element function. \return The vector
    containing the function evaluated at all the nodes in the mesh (in order to
    represent discontinuities, each node is included as many times as many
    elements share it).
  */
  auto nodal(const dense_t &solution) const -> dense_t;
  /**
    \brief Compute the L^2 and H^1_0 errors.
    \param dof The vector of degrees of freedom of the solution.
    \param solution Any callable object representing the solution function (a
    Field can be used by passing a lambda). \param gradient Any callable object
    representing the gradient of the solution function (it must return a
    point_t, a Field can be used by passing an arbitrary amount of lambdas)
    \return The error is returned as a pair where: error.first is the L^2 error
    and error.second is the H^1_0 error.
  */
  template <typename Quadrature, typename Solution, typename Gradient>
  auto error(const dense_t &dof, const Solution &solution,
             const Gradient &gradient) const -> std::pair<scalar_t, scalar_t>;

 private:
  /// The Mesh on which the FeSpace is built.
  const Mesh &_mesh;
  /// The number of degrees of freedom
  int _ndof{0};
  /// The vector containing the FeElement's
  std::vector<FeElement> _elements;
  /// The vector containing the FeFace's
  std::vector<FeFace> _faces;
};

//------------------------------------------------------------------------------
inline auto FeSpace::mesh() const -> const Mesh & { return _mesh; }

inline auto FeSpace::ndof() const -> int { return _ndof; }

inline auto FeSpace::operator()(const Element &element) const
    -> const FeElement & {
  return _elements[element.index()];
}

inline auto FeSpace::operator()(const Face &face) const -> const FeFace & {
  return _faces[face.index()];
}

template <typename Quadrature, typename Solution, typename Gradient>
auto FeSpace::error(const dense_t &dof, const Solution &solution,
                    const Gradient &gradient) const
    -> std::pair<scalar_t, scalar_t> {
#ifdef VERBOSE
  utl::Timer timer("error");
#endif

  scalar_t L2{0.0};
  scalar_t H10{0.0};

  // prefetch mesh
  const auto &mesh = this->mesh();

  // prefetch quadrature information
  constexpr int QSIZE = Quadrature::size();
  constexpr scalar_t reference_measure = Quadrature::measure();
  const auto weights = Quadrature::weights();
  const auto nodes = Quadrature::nodes();
  // helper to treat vectors as matrices (like sub2ind)
  auto at = [&QSIZE](auto q, auto i) { return q + i * QSIZE; };

  // loop over the elements
  for (ElementIterator e(mesh); !e.end(); ++e) {
    const auto &fe = (*this)(*e);

    // loop over the triangulation of an element
    for (const auto &[jacobian, translation, simplex_measure] :
         fe.triangulation()) {
      // precompute ratio
      scalar_t measure_ratio = simplex_measure / reference_measure;
      // precompute transformed nodes
      std::array<point_t, QSIZE> x;
      std::array<point_t, QSIZE> xf;
      std::array<scalar_t, QSIZE> solution_values;
      std::array<point_t, QSIZE> gradient_values;

      // take the nodes back to the reference space and compute the solution and
      // gradient at those points
      int index = 0;
      for (const auto &node : nodes) {
        xf[index] = jacobian * node + translation;
        x[index] = fe.map(xf[index]);
        solution_values[index] = solution(xf[index]);
        gradient_values[index] = gradient(xf[index]);
        ++index;
      }

      // precompute basis functions on nodes
      std::vector<scalar_t> phi = fe.eval_basis(x);
      std::vector<point_t> dphi = fe.eval_basis_derivatives(x);

      // loop over quadrature nodes
      for (int q = 0; q < QSIZE; ++q) {
        point_t du;
        du.setZero();
        scalar_t u{0.0};

        // reconstruct the approximation and its gradient at the nodes
        for (int i = 0; i < fe.ndof(); ++i) {
          du += dphi[at(q, i)] * dof[fe.dof(i)];
          u += phi[at(q, i)] * dof[fe.dof(i)];
        }

        // error computation
        L2 += (solution_values[q] - u) * (solution_values[q] - u) * weights[q] *
              measure_ratio;
        H10 += (gradient_values[q] - du).dot(gradient_values[q] - du) *
               weights[q] * measure_ratio;
      }
    }
  }
  return std::make_pair(std::sqrt(L2), std::sqrt(H10));
}

}  // namespace wave