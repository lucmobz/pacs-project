#pragma once

#include "Element.hpp"
#include "Mesh.hpp"
#include "MeshEntity.hpp"
#include "MeshEntityIterator.hpp"
#include "Traits.hpp"

namespace wave {
/*!
  \brief A Face is a MeshEntity of topological dimension 2 (in 3d).
*/
class Face : public MeshEntity {
 public:
  using point_t = Traits::point_t;
  /**
    \brief Constructs a Face.
    \param mesh The Mesh to which the Face belongs.
    \param index The global index of the Face in the Mesh.
  */
  Face(const Mesh& mesh, int index);
  /// \brief Get the normal (orientation depends on the order of the Node's).
  auto normal() const -> point_t;
  /**
    \brief Get the normal oriented a certain way.
    \param orientation The normal will be oriented on the side of the specified
    point (if not coplanar).
  */
  auto normal(const point_t& orientation) const -> point_t;
  /// \brief Tests if the Face is interior.
  auto is_interior() const -> bool;
  /// \brief Tests if the Face is exterior.
  auto is_exterior() const -> bool;
  /// \brief Get the boundary label (from 1 on) to which the Face belongs, 0 if
  /// interior.
  auto boundary() const -> int;
};

/// Specialization of a MeshEntityIterator on faces.
using FaceIterator = MeshEntityIterator<Face>;

//------------------------------------------------------------------------------
inline Face::Face(const Mesh& mesh, int index)
    : MeshEntity{mesh, mesh.dim() - 1, index} {}

inline auto Face::is_interior() const -> bool {
  // get the second element associated to the face and test if it is valid (if
  // the index is positive, if negative -index is the label of the boundary)
  ElementIterator en(*this, 1);
  return en->is_valid();
}

inline auto Face::is_exterior() const -> bool { return !is_interior(); }

inline auto Face::boundary() const -> int {
  if (is_interior())
    // convention to give boundary 0 to interior faces
    return 0;
  else {
    // if the face is exterior the second element associated to the face will be
    // an invalid element (with a negative index), -index will be the label of
    // the boundary. This is dependent on the mesh format. The decision was to
    // keep invalid elements so that the number of faces and neighbors of an
    // elments stay the same (some neighbors will therefor be invalid).
    ElementIterator en(*this, 1);
    return -en->index();
  }
}

}  // namespace wave