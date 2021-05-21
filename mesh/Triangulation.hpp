#pragma once

#include <tuple>
#include <vector>
#include "Element.hpp"
#include "Face.hpp"
#include "MeshEntity.hpp"
#include "Traits.hpp"
#include "mesh_type_traits.hpp"

namespace wave {
/**
  \brief A Triangulation is a partition of an entity (Element or Face) into
  simplices.

  It consists of a vector of tuples, each tuple has a matrix that is the
  jacobian of the affine transformation that maps the reference simplex into the
  physical one, a vector that represents the translation in the affine
  transofrmation, and the measure of the physical tetrahedron. It can be
  iterated using a range for loop with tuple unpacking.
*/
// Class template most methods are defined in body, except explicitely defined
// specializations
template <typename Entity>
class Triangulation {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;
  using entity_t = Entity;
  using matrix_t = Traits::matrix_t<entity_dim<Entity>>;
  /// Constant to exclude tetrahedral sleeves or too small simplices.
  static constexpr scalar_t EPSCGAL = Traits::EPSCGAL;

  // copy operations have been deleted to be sure to avoid copies
  Triangulation() = default;
  /// \brief Constructs a triangulation for the given entity (Element or Face).
  Triangulation(const entity_t& entity);
  Triangulation(const Triangulation&) = delete;
  Triangulation(Triangulation&&) = default;
  ~Triangulation() = default;
  Triangulation& operator=(const Triangulation&) = delete;
  Triangulation& operator=(Triangulation&&) = default;

  /// \brief Returns the tuple that represents the affine transformation for the
  /// simplex at given index.
  auto map(int index) const -> const std::tuple<matrix_t, point_t, scalar_t>& {
    return _maps[index];
  }
  /// \brief Returns the jacobian of the affine transformation for the simplex
  /// at given index.
  auto jacobian(int index) const -> const point_t& {
    return std::get<0>(_maps[index]);
  }
  /// \brief Returns the translation vector in the affine transfomration for the
  /// simplex at given index.
  auto translation(int index) const -> const point_t& {
    return std::get<1>(_maps[index]);
  }
  /// \brief Returns the measure of the simplex at given index.
  auto measure(int index) const -> const scalar_t& {
    return std::get<2>(_maps[index]);
  }
  /// \brief Returns the number of simplices in the triangulation
  auto size() const -> int { return _maps.size(); }
  /// \brief Tests if the triangulation is empty.
  auto empty() const -> bool { return _maps.empty(); }
  /// \brief Returns an iterator to the first simplex affine map.
  decltype(auto) begin() const { return _maps.begin(); }
  /// \brief Returns an iterator past the last simplex affine map.
  decltype(auto) end() const { return _maps.end(); }

 private:
  /// Contains the affine transformation in a vector of tuples, the first is the
  /// jacobian, then the translation vector then the measure of the simplex.
  std::vector<std::tuple<matrix_t, point_t, scalar_t>> _maps{};
};

template <>
Triangulation<Face>::Triangulation(const Face& face);

template <>
Triangulation<Element>::Triangulation(const Element& element);

}  // namespace wave