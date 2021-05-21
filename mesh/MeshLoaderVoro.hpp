#pragma once

#include <string>
#include "Mesh.hpp"
#include "MeshLoader.hpp"

namespace wave {
/**
  \brief Concrete MeshLoader that loads a file into a Mesh (file written in the
  VoroCrust format .voro)

  The VoroCrust format .voro consists in a mesh file where cells are listed with
  nodes, faces and neighbors in a sequence. Therefore node coordiantes are
  inserted as many time as they appear (otherwise comparing floating point
  number would be too costly). While connectivities are manipulated to have in
  the mesh: the set of all elements, the set of all faces (counted only once),
  the set of all neighbors (where negative indices in the connectivity signify
  boundary neighbors, this way the number of faces and number of neighbors is
  the same, and the negative index can be used as a label for the boundary
  faces). Also the set of faces for each element and the set of elements that
  share a face are represented.
*/
class MeshLoaderVoro : public MeshLoader {
 public:
  /// \brief Loads a file (given as file path) in the VoroCrust .voro format
  /// into a Mesh.
  virtual void load(const std::string& file, Mesh& mesh) const override;
};

}  // namespace wave
