#pragma once

#include "Mesh.hpp"
#include "Traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief A MeshEntity is an abstract entity in a Mesh.

  It is defined by a reference to a Mesh, a topological dimension and a global
  index into the set of all entities of given topological dimension (a
  MeshConnectivity object stored inside the Mesh). A MeshEntity doesn't contain
  any other data and provides methods to access entity specific information
  stored inside the Mesh. It can be thought as a view of a Mesh or as a fancy
  index.
*/
class MeshEntity {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;
  /**
    \brief Constructs a MeshEntity belonging to a Mesh.
    \param mesh The Mesh to which the MeshEntity belongs.
    \param dim The topological dimension.
    \param index The global index.
  */
  MeshEntity(const Mesh& mesh, int dim, int index);
  /// \brief Pure virtual destructor to make the class abstract.
  virtual ~MeshEntity() = 0;
  /// \brief Get the Mesh to which the MeshEntity belongs.
  auto mesh() const -> const Mesh&;
  /// \brief Get the topological dimension.
  auto dim() const -> int;
  /// \brief Get the number of Node's connected to the MeshEntity.
  auto nnodes() const -> int;
  /// \brief Get the global index.
  auto index() const -> int;
  /**
    \brief Get the local index of the MeshEntity into the connectivity of
    another specified MeshEntity. \param entity A MeshEntity assumed to be
    adjacent to the current one. \return Returns the local index or -1, if the
    current MeshEntity is not found in the connectivity of the specified
    MeshEntity.
  */
  auto index(const MeshEntity& entity) const -> int;
  /*
    \brief Tests if the current MeshEntity has a positive index.

    Negative index entities appear if an element has a boundary face. The
    neighbor through that Face is an entity with negative index. This way the
    number of Face's coincides with the number of neighbors for all Element's
    and the negative index can be used to label boundary Face's.
  */
  auto is_valid() const -> bool;
  /// \brief Get the diameter of the element.
  auto diameter() const -> scalar_t;
  /// \brief Get the barycent of the element.
  auto barycenter() const -> point_t;
  /**
    \brief Get the bounding box bounds.
    \return Returns a pair of points, the minimal and maximal box coordinates.
  */
  auto bounds() const -> std::pair<point_t, point_t>;

 protected:
  /// The mesh to which the entity belongs
  const Mesh& _mesh;
  /// The topological dimension of the entity.
  int _dim{0};
  /// The global index of the entity.
  int _index{0};

  // A mesh iterator allows iterating on the mesh and requires access to the
  // internals.
  template <typename Entity>
  friend class MeshEntityIterator;
};

//------------------------------------------------------------------------------
inline MeshEntity::MeshEntity(const Mesh& mesh, int dim, int index)
    : _mesh{mesh}, _dim{dim}, _index{index} {}

inline auto MeshEntity::mesh() const -> const Mesh& { return _mesh; }

inline auto MeshEntity::dim() const -> int { return _dim; }

inline auto MeshEntity::index() const -> int { return _index; }

inline auto MeshEntity::is_valid() const -> bool { return _index >= 0; }

inline auto MeshEntity::nnodes() const -> int {
  // if the entity is valid and is not a node return the size of the entity to
  // node connectivity of the specific entity. Otherwise return 1 if it's a
  // node. Else return 0. Exceptions have been avoided to not slow down the code
  // (which slows down by 25% in some points).
  if (this->is_valid() && _dim > 0)
    return _mesh.topology(_dim, 0).size(_index);
  else if (_dim == 0)
    return 1;
  else
    return 0;
}

}  // namespace wave