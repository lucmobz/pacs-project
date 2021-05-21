#include "MeshLoader.hpp"
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include "Mesh.hpp"
#include "MeshLoaderVoro.hpp"

namespace wave {

// factory method, it returns a unique pointer to a concrete loader, or the
// nullptr if the format is not supported
auto MeshLoader::make(const std::string& format)
    -> std::unique_ptr<MeshLoader> {
  if (format == "voro") {
    return std::make_unique<MeshLoaderVoro>();
    // to add new formats modify this if statement
  } else {
    std::cerr << "format error\n";
    return nullptr;
  }
}
// skips the number of lines in the stream, setting it at the beginning of the
// new line past those.
auto MeshLoader::skip(std::ifstream& input, int nlines) -> std::ifstream& {
  for (int i = 0; i < nlines; ++i)
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  return input;
}
// finds the target line as a substrnig of a given line and sets the stream at
// the beginning of the new line
auto MeshLoader::skip(std::ifstream& input, const std::string& target_line)
    -> std::ifstream& {
  std::string line;
  while (input) {
    std::getline(input, line);
    if (line.find(target_line) != std::string::npos) {
      break;
    }
  }
  return input;
}

void MeshLoader::finalize(Mesh& mesh) {
  mesh._nnodes = mesh._nodes.size() / mesh._dim;
  mesh._nelements = mesh._elements.size();
  mesh._nfaces = mesh._faces.size();
}

}  // namespace wave