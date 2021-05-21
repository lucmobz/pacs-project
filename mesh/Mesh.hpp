#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "MeshConnectivity.hpp"
#include "Traits.hpp"

namespace wave {

/**
  \brief A Mesh stores geometrical and topological informations sequentially.

  Geometrical information (the coordinates of the Node's) is stored in a vector
  contiguously and can be accessed both globally from the Mesh and through
  Node's objects (which themselves don't store anything a part from their
  index). Topological informations, i.e. the connectivity relations among
  MeshEntity's are stored in MeshConnectivity's objects (which store them
  contiguously) and can be accessed by specifying the topological dimensions
  associated to the connectivity relations (3d to 3d for the element to element
  connectivty, 2d to 0d for the face to node, 3d to 2d for the element to face
  etc...) Iteration over the Mesh is done using MeshEntityIterator's, which also
  allow to iterate of a subset of MeshEntity's adjacent to a MeshEntity. A Mesh
  can be used to output itself in VTK format, and the same method can be used to
  output associated data as well (such as a solution computed at the nodes).
*/
class Mesh {
 public:
  using scalar_t = Traits::scalar_t;
  /// A map to a point (as an Eigen map), this way the node coordinates stored
  /// in the vector are not copied
  using mpoint_t = Traits::mpoint_t;
  static constexpr int DIM = Traits::DIM;
  /// The VTK default cell type id for a polyhedron needed in the VTK format.
  static constexpr int VTK_CELL_TYPE = 42;
  /**
    \brief Constructs a mesh from a file and a dimension (can be omitted).

    The Mesh file must have a compatible format (such as .voro), the format is
    extracted from the path and using a MeshLoader the mesh is constructed
    through a factory method of the MeshLoader (for example it can call the
    concrete MeshLoaderVoro for the voro format). This way only the MeshLoader
    has to be modified if new formats are to be supported.

    \param file The string containing the file path or the file name.
    \param dim The dimension of the Mesh.
  */
  explicit Mesh(const std::string& file, int dim = DIM);

  /// \brief Get the dimension.
  auto dim() const -> int;
  /// \brief Get the number of MeshEntity's of given topological dimension.
  auto size(int dim) const -> int;
  /// \brief Get the number of Node's.
  auto nnodes() const -> int;
  /// \brief Get the number of Element's.
  auto nelements() const -> int;
  /// \brief Get the number of Face's.
  auto nfaces() const -> int;
  /// \brief Get the vector of all node coordinates.
  auto geometry() const -> const std::vector<scalar_t>&;
  /// \brief Get the point (as a view, Eigen map) of given index.
  auto geometry(int index) const -> mpoint_t;
  /// \brief Get the coordinate i of point of given index.
  auto geometry(int index, int i) const -> scalar_t;
  /**
    \brief Get the connectivity relation between the two topological dimensions.

    3d to 0d is the element to node connectivity.
    2d to 0d is the face to node.
    3d to 3d is the element to element (neighbors).
    3d to 2d is the element to face (the face indices for a given element).
    2d to 3d is the face to element (the element indices that share a face).

    \param d0 First topological dimension.
    \param d1 Second topological dimension
    \return A MeshConnectivity object representing the relation
  */
  auto topology(int d0, int d1) const -> const MeshConnectivity&;
  /**
    \brief Output the Mesh in VTK format (use ParaView to visualize).

    \param file The name of the output file without the format extension (the
    extesion will be appended) \param cell_type VTK cell type, leave as default
    for general polyhedron mesh.
  */
  void vtk(const std::string& file, int cell_type = VTK_CELL_TYPE) const;
  /**
    \brief Output the Mesh in VTK format as well as the nodal data specified

    \tparam Data The container type for the nodal data (a vector for example)
    \param file The file name to outuput without extension.
    \param data The container for nodal data.
    \param cell_type VTK cell type, leave as default for general polyhedron
    mesh.
  */
  template <typename Data>
  void vtk(const std::string& file, const Data& data,
           int cell_type = VTK_CELL_TYPE) const;

 private:
  /// Mesh dimension
  int _dim{0};
  /// Number of Node's
  int _nnodes{0};
  /// Number of Element's
  int _nelements{0};
  /// Number of Face's
  int _nfaces{0};

  /// Stores contiguously the coordinates of the nodes (x0, y0, z0, x1, y1, z1,
  /// ...).
  std::vector<scalar_t> _nodes{};
  /// 3d to 0d connectivity, element to node.
  MeshConnectivity _elements{DIM, 0};
  /// 3d to 3d connectivity, element to element (neighbors).
  MeshConnectivity _neighbors{DIM, DIM};
  /// 2d to 0d connectivity, face to node.
  MeshConnectivity _faces{DIM - 1, 0};
  /// 3d to 0 connectivity, element to face.
  MeshConnectivity _element_faces{DIM, DIM - 1};
  /// 2d to 3d connectivity, face to element.
  MeshConnectivity _face_elements{DIM - 1, DIM};

  /// \brief Get the index in the node vector of the coordinate i of node at the
  /// specified index.
  static constexpr auto _pos(int index, int i) -> int;

  /// A MeshLoader needs to access the Mesh internals to fill the Mesh.
  friend class MeshLoader;
};

//------------------------------------------------------------------------------
inline auto Mesh::dim() const -> int { return _dim; }

inline auto Mesh::nnodes() const -> int { return _nnodes; }

inline auto Mesh::nelements() const -> int { return _nelements; }

inline auto Mesh::nfaces() const -> int { return _nfaces; }

inline auto Mesh::size(int dim) const -> int {
  // return the number of topological entities of given dimension
  if (dim == 0)
    return _nnodes;
  else if (dim == _dim)
    return _nelements;
  else if (dim == _dim - 1)
    return _nfaces;
  else
    // return 0 as default case, exceptions slow down the code
    return 0;
}

inline auto Mesh::geometry() const -> const std::vector<scalar_t>& {
  return _nodes;
}

inline auto Mesh::geometry(int index) const -> mpoint_t {
  // return a view of a point (the actual data stays in the vector and is not
  // copied)
  return mpoint_t{_nodes.data() + _pos(index, 0)};
}

inline auto Mesh::geometry(int index, int i) const -> scalar_t {
  return _nodes[_pos(index, i)];
}

inline auto Mesh::topology(int d0, int d1) const -> const MeshConnectivity& {
  if (d0 == _dim && d1 == 0)
    return _elements;
  else if (d0 == _dim && d1 == _dim)
    return _neighbors;
  else if (d0 == _dim - 1 && d1 == 0)
    return _faces;
  else if (d0 == _dim && d1 == _dim - 1)
    return _element_faces;
  else if (d0 == _dim - 1 && d1 == _dim)
    return _face_elements;
  else
    // return the element connectivity as default, because exceptions slow down
    // the code
    return _elements;
}

template <typename Data>
void Mesh::vtk(const std::string& file, const Data& data, int cell_type) const {
  // fill the file with the mesh information and then append the nodal data
  // information
  this->vtk(file, cell_type);
  std::ofstream output(file, std::ios_base::app);

  output << "POINT_DATA " << data.size() << "\n"
         << "SCALARS solution double\n"
         << "LOOKUP_TABLE default\n";

  for (int index = 0; index < data.size(); ++index) {
    output << data[index] << "\n";
  }

  output.close();
}

constexpr auto Mesh::_pos(int index, int i) -> int { return index * DIM + i; }

}  // namespace wave