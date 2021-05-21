#pragma once

#include "Mesh.hpp"
#include "MeshEntity.hpp"
#include "MeshEntityIterator.hpp"
#include "Traits.hpp"

namespace wave {
/**
  \brief A Node is a MeshEntity of topological dimension 0.
*/
class Node : public MeshEntity {
 public:
  using scalar_t = Traits::scalar_t;
  /// A map to a point (as in Eigen), accessing coordinates will return a view
  /// of the point coordinates in the mesh and will not copy it.
  using mpoint_t = Traits::mpoint_t;
  /**
    \brief Constructs an Node.
    \param mesh The Mesh to which the Node belongs.
    \param index The global index of the Node in the Mesh.
  */
  Node(const Mesh& mesh, int index);
  /// \brief Get the coordinates of the Node (as a view on the total coordinates
  /// vector stored in the Mesh).
  auto point() const -> mpoint_t;
  /// \brief Get coordinate i of the Node.
  auto point(int i) const -> scalar_t;
  /// \brief Get coordinate i of the Node.
  auto operator[](int i) const -> scalar_t;
};

/// Specialization of a MeshEntityIterator on nodes.
using NodeIterator = MeshEntityIterator<Node>;

//------------------------------------------------------------------------------
inline Node::Node(const Mesh& mesh, int index) : MeshEntity{mesh, 0, index} {}

inline auto Node::point() const -> mpoint_t { return _mesh.geometry(_index); }

inline auto Node::point(int i) const -> scalar_t {
  // nodes are stored in a vector inside the mesh
  return _mesh.geometry(_index, i);
}

inline auto Node::operator[](int i) const -> scalar_t {
  return _mesh.geometry(_index, i);
}

}  // namespace wave