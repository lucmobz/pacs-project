#pragma once

#include <span>
#include "Mesh.hpp"
#include "MeshEntity.hpp"

namespace wave {

/**
  \brief A MeshEntityIterator provides a way to iterate over MeshEntity's
  belonging to a Mesh.

  Such iterator can be used both to iterate over Element's, Face's and Node's in
  a Mesh, and to iterate over the adjacent entities of a given entity. The class
  is a template that is specialized for actual concrete entities (NodeIterator,
  FaceIterator and ElementIterator). The syntax is the following: for
  (ElementIterator e(mesh); !e.end(); ++e) { e->foo(); } This iterates over all
  the elements in the mesh. Another example is iterating over all the neighbors
  of an element: ElmentIterator e(mesh, 0) // get an iterator to the first
  element of the mesh for (ElementIterator en(*e); !en.end(); ++en) { en->foo();
  } The syntax for looping over the nodes of face is: FaceIterator f(mesh) //
  get an iterator to the first face of the mesh for (NodeIterator n(*f);
  !n.end(); ++n) { n->foo(); }
*/
// Class template all definitions are in body.
template <typename Entity>
class MeshEntityIterator {
 public:
  using entity_t = Entity;

  /**
    \brief Constructs an iterator over a mesh starting at the given index.
    \param mesh The Mesh to iterate on.
    \param index The index of the entity to start iterating.
  */
  explicit MeshEntityIterator(const Mesh& mesh, int index = 0)
      : _index{index}, _entity{mesh, index} {
    _index_end = mesh.size(_entity.dim());
  };
  /**
    \brief Constructs an iterator over the local sub-entities of an entity
    \param entity The entity on which to iterate.
    \param index The starting local index of the sub-entity iterator.
  */
  explicit MeshEntityIterator(const MeshEntity& entity, int index = 0)
      : _index{index}, _entity{entity.mesh(), index} {
    const MeshConnectivity& connectivity =
        entity.mesh().topology(entity.dim(), _entity.dim());
    _index_end = connectivity.size(entity.index());
    _index_map = connectivity(entity.index());
  };
  /// \brief Tests if the iterator has reached the end (or is invalid with a
  /// negative index)
  // it works both for iterating forward and backward
  auto end() const -> bool { return _index >= _index_end || _index < 0; }
  /// \brief Dereferencing returns the unerlying MeshEntity.
  auto operator*() -> entity_t& { return *(this->operator->()); }
  /// \brief Allows access to methods of a MeshEntity.
  // the arrow operator overload reapplies itself so it needs to return a
  // pointer
  auto operator-> () -> entity_t* {
    // test if the iterator is local to an entity or global to the mesh, sets
    // the underlying entity index accordingly and returns the entity
    _entity._index = _index_map.empty() ? _index : _index_map[_index];
    return &_entity;
  }
  /// \brief Prefix increment operator.
  auto operator++() -> MeshEntityIterator& {
    ++_index;
    return *this;
  }
  /// \brief Postfix increment operator.
  auto operator++(int) -> MeshEntityIterator {
    MeshEntityIterator other = *this;
    ++*this;
    return other;
  }
  /// \brief Advance the iterator by n positions.
  auto operator+=(int n) -> MeshEntityIterator& {
    _index += n;
    return *this;
  }
  /// \brief Prefix decrement operator
  auto operator--() -> MeshEntityIterator& {
    --_index;
    return *this;
  }
  /// \brief Postfix decrement operator
  auto operator--(int) -> MeshEntityIterator {
    MeshEntityIterator other = *this;
    --*this;
    return other;
  }
  /// \brief Moves iterator back by n positions (if the index becomes negative
  /// negative it will become the end iterator)
  auto operator-=(int n) -> MeshEntityIterator& {
    _index -= n;
    return *this;
  }
  /// \brief plus operator to return an iterator advanced by n positions
  friend auto operator+(const MeshEntityIterator& lhs, int n)
      -> MeshEntityIterator {
    MeshEntityIterator ans = lhs;
    ans += n;
    return ans;
  }
  /// \brief minus operator to return an iterator falling back by n positions
  /// (if index becomes negative it will become the end iterator)
  friend auto operator-(const MeshEntityIterator& lhs, int n)
      -> MeshEntityIterator {
    MeshEntityIterator ans = lhs;
    ans -= n;
    return ans;
  }

 private:
  /// Current local index of the iterator
  int _index{0};
  /// Final index in the set of entities on which one is iterating
  int _index_end{0};
  /// Underlying MeshEntity (entities are implementd as fancy indices they don't
  /// store data)
  entity_t _entity;
  /// A view on the local connectivity of the entity on which one is iterating
  /// (empty span if one is iterating on all the entities in a mesh of a
  /// specific dimension)
  std::span<const int> _index_map{};
};

}  // namespace wave