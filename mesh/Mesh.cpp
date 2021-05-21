#include "Mesh.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include "Element.hpp"
#include "MeshFunction.hpp"
#include "MeshLoader.hpp"
#include "utilities_include.hpp"

namespace wave {

Mesh::Mesh(const std::string& file, int dim)
    : _dim{dim},
      _elements{dim, 0},
      _neighbors{dim, dim},
      _faces{dim - 1, 0},
      _element_faces{dim, dim - 1},
      _face_elements{dim - 1, dim} {
  // figure out the format from the file path (tried regex and filesystem, but
  // is super buggy)
  std::string format = file.substr(file.find_last_of(".") + 1);
  // use a factory method to build the proper concrete loader and load the mesh
  MeshLoader::make(format)->load(file, *this);
}

void Mesh::vtk(const std::string& file, int cell_type) const {
  std::ofstream output(file);

  output << "# vtk DataFile Version 3.0\n"
         << "output\n"
         << "ASCII\n"
         << "DATASET UNSTRUCTURED_GRID\n";

  output << "POINTS " << _nnodes << " double\n";

  // insert node coordinates
  for (int index = 0; index < _nnodes; ++index) {
    for (int i = 0; i < _dim; ++i) output << geometry(index, i) << " ";

    // padding with extra zeros according to mesh dimension
    if (_dim == 1)
      output << 0.0 << " " << 0.0 << "\n";
    else if (_dim == 2)
      output << 0.0 << "\n";
    else
      output << "\n";
  }

  // Here the total number of entries needs to be computed for the VTK format,
  // this means that the total number of numbers appearing after CELLS needs to
  // be known according to the specific VTK_CELL_TYPE that is used (42 is for
  // polyehdrons in this case). Each row will have a number telling how many
  // numbers are in that row, then the number of faces, then the number of face
  // nodes, per each face, then the nodes indices. To compute this, a
  // composition of MeshFunction objects can be used, avoiding to write a loop.
  // Each MeshFunction applies a function to the entities in a mesh (the type of
  // entity is derived from the lambda), and stores the values in a vector that
  // can be manipulated later, MeshFunctions can also be called on entities, in
  // that case only the adjacent entities will be considered. Here a
  // MeshFunction is built on the elements, and for each element another mesh
  // function is built on the faces.
  MeshFunction nelement_faces(*this, [](const Element& element) {
    MeshFunction nface_nodes(element,
                             [](const Face& face) { return face.nnodes(); });
    return 2 + element.nfaces() +
           std::accumulate(nface_nodes.begin(), nface_nodes.end(), 0);
  });

  // the total entry can be obtained by summing all the row entries
  int nentries =
      std::accumulate(nelement_faces.begin(), nelement_faces.end(), 0);

  const auto& element_faces = topology(_dim, _dim - 1);
  const auto& faces = topology(_dim - 1, 0);

  // output cell information following the VTK format for polyhedrons
  output << "CELLS " << _nelements << " " << nentries << "\n";

  // for each cell output the number of entry on the line, the number of faces,
  // the number of nodes of each face, followed by the node indices
  for (int index = 0; index < _nelements; ++index) {
    int nfaces = element_faces.size(index);
    int nface_nodes = 0;

    for (const auto& f : element_faces(index)) nface_nodes += faces.size(f);
    output << (1 + nfaces + nface_nodes) << " " << nfaces << " ";

    for (const auto& f : element_faces(index)) {
      output << faces.size(f) << " ";
      for (const auto& n : faces(f)) output << n << " ";
    }
    output << "\n";
  }

  // output cell types as per VTK format
  output << "CELL_TYPES " << _nelements << "\n";
  for (int index = 0; index < _nelements; ++index) output << cell_type << "\n";

  output.close();
}

}  // namespace wave