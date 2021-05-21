#pragma once

#include <type_traits>
#include <utility>
#include <vector>
#include "Mesh.hpp"
#include "MeshEntityIterator.hpp"
#include "mesh_type_traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief A MeshFunction applies a lambda to a set of entities and stores the
  resulting values.

  Instead of looping over the mesh entities it's possible to create a mesh
  function that applies a lambda to the entities and stores the values, all the
  relevant info is deduced so that no template parameters are needed. For
  example if the lambda takes an Element as input, the lambda will be applied to
  all elements in the Mesh. It also works locally on the set of sub-entities of
  an entity.
*/
// the class is a template to allow deduction of types using type traits (class
// is an experiment, maybe remove it later)
template <typename Function>
class MeshFunction {
 public:
  // use type trait to get the type of the first argument of the lambda to
  // deduce the entity type
  using entity_t = std::remove_cvref_t<utl::first_argument<Function>>;
  // deduce the return type from the lambda in order to get the type to store in
  // the vector
  using value_t = decltype(std::declval<Function>()(std::declval<entity_t>()));
  using entity_iterator_t = MeshEntityIterator<entity_t>;
  // use type traits to figure out the entity dimension
  static constexpr int DIM = entity_dim<entity_t>;

 public:
  /**
    \brief Constructs a MeshFunction on a Mesh by taking a lambda

    From the lambda argument type the entity type is deduce, the lambda then
    applies it self to all the entities in the mesh and the values are stored in
    a vector, that can be manipulated later. \param mesh The Mesh. \param
    function The lambda must have an Element, Face or Node as argument type.
  */
  MeshFunction(const Mesh& mesh, const Function& function)
      : _mesh{mesh}, _size{mesh.topology(DIM, 0).size()} {
    // compute all the values by iterating over the mesh
    _values.reserve(_size);
    for (entity_iterator_t it(_mesh); !it.end(); ++it) {
      _values.emplace_back(function(*it));
    }
  }
  /**
    \brief Constructs a MeshFunction on a MeshEntity (Node, Face, Element) by
    taking a lambda.

    From the lambda argument type the type of the subentity is deduced, this
    allows to apply the mesh function only to the sub-entities of the given
    entity (for example if entity is an Element and the lambda takes a Face as
    argument, the lambda is applied to all the faces of the element) \param
    entity The entity whose sub-entities will be passed to the lambda. \param
    function The lambda that will be called on the sub-entities.
  */
  MeshFunction(const MeshEntity& entity, const Function& function)
      : _mesh(entity.mesh()),
        _size(entity.mesh().topology(entity.dim(), DIM).size(entity.index())) {
    // compute all the values by iterating over the mesh
    _values.reserve(_size);
    for (entity_iterator_t it(entity); !it.end(); ++it) {
      _values.emplace_back(function(*it));
    }
  }
  /// \brief Get the mesh.
  auto mesh() const -> const Mesh& { return _mesh; }
  /// \brief Get the topological dimension of the entities
  auto dim() const -> int { return DIM; }
  /// \brief Get the number of entities.
  auto size() const -> int { return _size; }
  /// \brief Get the vector of values.
  auto values() const -> const std::vector<value_t>& { return _values; }
  /// \brief Get the specific value for the specific entity
  auto operator()(const entity_t& entity) const -> const value_t& {
    return _values[entity.index()];
  }
  /// \brief Get an iterator at the beginning of the set of values
  decltype(auto) begin() const { return _values.begin(); }
  /// \brief Get and iterator at the end of the set of values
  decltype(auto) end() const { return _values.end(); }

 private:
  /// The Mesh
  const Mesh& _mesh;
  /// The number of entities
  int _size{0};
  /// The vector of values after having called the function
  std::vector<value_t> _values{};
};

}  // namespace wave