#pragma once

#include <iostream>
#include "Mesh.hpp"
#include "MeshEntity.hpp"
#include "MeshEntityIterator.hpp"

namespace wave {
/*!
  \brief An Element is a MeshEntity of topological dimension 3 (in 3d).
*/
class Element : public MeshEntity {
 public:
  /**
    \brief Constructs an Element.
    \param mesh The Mesh to which the Element belongs.
    \param index The global index of the Element in the Mesh.
  */
  Element(const Mesh& mesh, int index);

  /// \brief Get the number of Face's.
  auto nfaces() const -> int;
  /// \brief Get the number of neighboring Element's.
  auto nneighbors() const -> int;
};

/// Specialization of a MeshEntityIterator on elements.
using ElementIterator = MeshEntityIterator<Element>;

//------------------------------------------------------------------------------
inline Element::Element(const Mesh& mesh, int index)
    : MeshEntity{mesh, mesh.dim(), index} {}

inline auto Element::nfaces() const -> int {
  // the number of faces is extracted from the 3d to 2d connectivity.
  return _mesh.topology(_dim, _dim - 1).size(_index);
}

inline auto Element::nneighbors() const -> int { return nfaces(); }

}  // namespace wave