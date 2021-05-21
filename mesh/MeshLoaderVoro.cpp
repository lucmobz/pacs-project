#include "MeshLoaderVoro.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "Mesh.hpp"
#include "MeshConnectivity.hpp"

namespace wave {
// this function just loads the mesh file into the mesh data structure,
// manipulating the VoroCrust .voro format.
void MeshLoaderVoro::load(const std::string& file, Mesh& mesh) const {
  std::ifstream input(file);

  if (input.is_open()) {
    const int dim = mesh.dim();
    int node_counter = 0;
    int element_counter = 0;
    int face_counter = 0;

    while (skip(input, "VoronoiCell")) {
      // add nodes to the mesh (nodes are duplicated in VoroCrust)
      int nnodes;
      input >> nnodes;
      skip(input, 1);

      ref_geometry(mesh).reserve(ref_geometry(mesh).size() + nnodes);
      std::copy_n(std::istream_iterator<double>(input), nnodes * dim,
                  std::back_inserter(ref_geometry(mesh)));

      // create all connectivities
      // create element to node connectivity
      std::vector<int> element(nnodes);
      std::iota(element.begin(), element.end(), node_counter);

      // for each face, create face to node connectivity
      int nfaces;
      input >> nfaces;
      skip(input, 1);

      std::vector<std::vector<int>> faces;
      faces.reserve(nfaces);

      for (int i = 0; i < nfaces; ++i) {
        std::string line;
        std::getline(input, line, '#');
        skip(input, 1);

        std::istringstream sline(line);
        std::istream_iterator<int> eos;

        std::vector<int> face(std::istream_iterator<int>(sline), eos);
        std::for_each(face.begin(), face.end(),
                      [&element](auto& i) { i = element[i]; });

        faces.emplace_back(std::move(face));
      }

      // create element to face connectivity, fill later
      std::vector<int> element_faces(nfaces);

      // create element to element connectivity
      std::vector<int> neighbors;
      std::copy_n(std::istream_iterator<int>(input), nfaces,
                  std::back_inserter(neighbors));

      // add all the connectivities to the mesh
      ref_topology(mesh, dim, 0).push(element);

      for (int i = 0; i < nfaces; ++i) {
        // Take care of duplicated faces, if the neighbor of the current element
        // through face i has an index that is lower then the one of the current
        // element, it means that face i already exists (with a different
        // orientation), because it was already added when such neighbor was the
        // current element. Therefore the face to node connectivity and the face
        // to element connectivity is not added. And the current, element to
        // face connectivity has to be modified accordingly. Otherwise the face
        // is added and the face counter is increased. The case of a negative
        // index (which means a boundary face) is also taken care of.
        if (neighbors[i] < element_counter && neighbors[i] >= 0) {
          const auto& connectivity = mesh.topology(dim, dim)(neighbors[i]);
          auto it = std::find(connectivity.begin(), connectivity.end(),
                              element_counter);
          int pos = std::distance(connectivity.begin(), it);

          element_faces[i] = mesh.topology(dim, dim - 1)(neighbors[i])[pos];
        } else {
          ref_topology(mesh, dim - 1, 0).push(faces[i]);
          // create face to element connectivity
          ref_topology(mesh, dim - 1, dim)
              .push(std::vector<int>{element_counter, neighbors[i]});
          element_faces[i] = face_counter;

          ++face_counter;
        }
      }

      ref_topology(mesh, dim, dim - 1).push(element_faces);
      ref_topology(mesh, dim, dim).push(neighbors);

      // increment counters
      node_counter += nnodes;
      ++element_counter;
    }
  } else
    std::cerr << "error in opening the file\n";

  input.close();

  finalize(mesh);
};

}  // namespace wave