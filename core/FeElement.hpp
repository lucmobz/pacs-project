#pragma once

#include <numeric>
#include <utility>
#include "Traits.hpp"
#include "functional_include.hpp"
#include "mesh_include.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief A FeElement is a finite element in a finite element space (FeSpace)

  This class is built on top of an underlying mesh Element (which is nothing but
  an index with methods to compute associated information) and stores the
  algorithm specific information associated with such element.
*/
class FeElement {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;

  /// \brief Construct a finite element from an Element and a polynomial degree
  FeElement(const Element &element, int degree);

  /// \brief Get the polynomial degree.
  auto degree() const -> int;
  /// \brief Get the number of degrees of freedom on the element.
  auto ndof() const -> int;
  /// \brief Get the indices of the degrees of freedom of the element.
  auto dof() const -> const std::vector<int> &;
  /// \brief Get the index of the i-th degree of freedom of the element.
  auto dof(int i) const -> int;
  /// \brief Compute the location in the reference cartesian bounding box from a
  /// physical point. \param x The physical point. \return The reference point
  /// in the box (-1,1)^3.
  auto map(const point_t &x) const -> point_t;
  /// \brief Get the jacobian of the transfomration that maps a reference point
  /// into a physical point.
  auto jacobian() const -> const point_t &;
  /// \brief Get the diameter of the underlying mesh Element.
  auto diameter() const -> scalar_t;
  /// \brief Set the diameter of the underlying mesh Element (normally by using
  /// the Element method).
  void diameter(scalar_t diameter);
  /// \brief Get the Triangulation of the Element used in quadrature rules.
  auto triangulation() const -> const Triangulation<Element> &;
  /// \brief Compute the Triangulation of the Element.
  void triangulate();

  /// \brief Compute the basis functions on a set of points.
  /// \tparam Points The container type for the set of points.
  /// \param points The container for the set of points.
  /// \return Returns the values of the basis functions in a vector of scalars.
  template <typename Points>
  auto eval_basis(const Points &points) const -> std::vector<scalar_t>;
  /// \brief Compute the basis function derivatives on a set of points.
  /// \tparam Points The container type for the set of points.
  /// \param points The container for the set of points.
  /// \return Returns the gradients of the basis functions computed at the
  /// points as a vector of points.
  template <typename Points>
  auto eval_basis_derivatives(const Points &points) const
      -> std::vector<point_t>;
  /// \brief Compute the i-th basis function on a point given by x.
  auto eval_basis(int i, const point_t &x) const -> scalar_t;
  /// \brief Compute the i-th basis function gradient on a point given by x.
  /// \return Returns the gradient as a point.
  auto eval_basis_derivatives(int i, const point_t &x) const -> point_t;

 private:
  /// The underlying mesh Element.
  Element _element;
  /// The indices of the degrees of freedom.
  std::vector<int> _dof;
  /// The polynomial degree
  int _degree{0};
  /// The reference to physical box map, represented as the jacobian and
  /// midpoint.
  std::pair<point_t, point_t> _map{point_t::Zero(), point_t::Zero()};
  /// The Triangulation of the underlying Element.
  Triangulation<Element> _triangulation;
  /// The diameter of the element.
  scalar_t _diameter{0.0};

  /// \brief Compute the number of degrees of freedom.
  static constexpr auto _compute_ndof(int dim, int degree) -> int;
};

//------------------------------------------------------------------------------
inline FeElement::FeElement(const Element &element, int degree)
    : _element(element), _degree(degree) {
  int ndof = _compute_ndof(element.dim(), degree);
  _dof.resize(ndof);
  std::iota(_dof.begin(), _dof.end(), element.index() * ndof);

  const auto &bounds = element.bounds();
  _map = std::make_pair((bounds.second - bounds.first) / 2,
                        (bounds.first + bounds.second) / 2);
}

inline auto FeElement::degree() const -> int { return _degree; }

inline auto FeElement::ndof() const -> int { return _dof.size(); }

inline auto FeElement::dof() const -> const std::vector<int> & { return _dof; }

inline auto FeElement::dof(int i) const -> int { return _dof[i]; }

inline auto FeElement::diameter() const -> scalar_t { return _diameter; }

inline void FeElement::diameter(scalar_t diameter) { _diameter = diameter; }

inline auto FeElement::triangulation() const -> const Triangulation<Element> & {
  assert(!_triangulation.empty());
  return _triangulation;
}

inline void FeElement::triangulate() {
  _triangulation = std::move(Triangulation<Element>(_element));
}

inline auto FeElement::map(const point_t &x) const -> point_t {
  return (x - _map.second).array() / _map.first.array();
}

inline auto FeElement::jacobian() const -> const point_t & {
  return _map.first;
}

constexpr auto FeElement::_compute_ndof(int dim, int degree) -> int {
  return utl::factorial(degree + dim) / utl::factorial(degree) /
         utl::factorial(dim);
}

inline auto FeElement::eval_basis(int i, const point_t &x) const -> scalar_t {
  return legendre::normalization(i) * legendre::phi(i, x);
}

inline auto FeElement::eval_basis_derivatives(int i, const point_t &x) const
    -> point_t {
  return point_t{
      legendre::normalization(i) * legendre::phix(i, x) / jacobian()[0],
      legendre::normalization(i) * legendre::phiy(i, x) / jacobian()[1],
      legendre::normalization(i) * legendre::phiz(i, x) / jacobian()[2]};
}

template <typename Points>
inline auto FeElement::eval_basis(const Points &points) const
    -> std::vector<scalar_t> {
  int npoints = points.size();
  std::vector<scalar_t> phi;
  phi.reserve(ndof() * npoints);

  for (int i = 0; i < ndof(); ++i) {
    for (int q = 0; q < npoints; ++q) {
      phi.emplace_back(eval_basis(i, points[q]));
    }
  }
  return phi;
}

template <typename Points>
inline auto FeElement::eval_basis_derivatives(const Points &points) const
    -> std::vector<point_t> {
  int npoints = points.size();
  std::vector<point_t> dphi;
  dphi.reserve(ndof() * npoints);

  for (int i = 0; i < ndof(); ++i) {
    for (int q = 0; q < npoints; ++q) {
      dphi.emplace_back(eval_basis_derivatives(i, points[q]));
    }
  }
  return dphi;
}

}  // namespace wave