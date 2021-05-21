#pragma once

#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include "Mesh.hpp"
#include "MeshConnectivity.hpp"

namespace wave {
/**
  \brief A MeshLoader is an abstract loader for loading mesh files.

  A pure abstract class that represents an object with access to the internals
  of a Mesh. using the make method a concrete loader can be instantiated (such
  as MeshLoader Voro), which can then load a Mesh according to the specific
  format. Addition of a new concrete loader for a different mesh format just
  requires to add the concrete loader to the make method and implement it.
*/
class MeshLoader {
 public:
  /// \brief virtual Destructor to have inheritance.
  virtual ~MeshLoader() = default;
  /**
    \brief Creates a specific concrete loader according to the format specified.
    \param format The string containing the format.
    \return Returns a unique pointer to the concrete loader instance.
  */
  static auto make(const std::string& format) -> std::unique_ptr<MeshLoader>;
  /**
    \brief Pure virtual load method.
    \param file The string containing the file path.
    \param mesh The Mesh to load the file into.
  */
  virtual void load(const std::string& file, Mesh& mesh) const = 0;

 protected:
  /// \brief Skips the given number of lines from the input stream.
  static auto skip(std::ifstream& input, int nlines) -> std::ifstream&;
  /// \brief Sets the input stream past the target line
  static auto skip(std::ifstream& input, const std::string& target_line)
      -> std::ifstream&;
  /// \brief Access the internal geometrical information in a Mesh (the
  /// coordinate vector).
  static auto ref_geometry(Mesh& mesh) -> std::vector<double>&;
  /**
    \brief Access the internal topological connectivity information in a Mesh.

    For example to get access to the elements in a mesh set d0 = 3, d1 = 0.
    \param mesh The Mesh.
    \param d0 The first topological dimension.
    \param d1 The second topological dimension.
    \return Returns a reference to a MeshConnectivity inside the Mesh.
  */
  static auto ref_topology(Mesh& mesh, int d0, int d1) -> MeshConnectivity&;
  /// \brief Finalizes the mesh by computing the sizes and last parameters.
  static void finalize(Mesh& mesh);
};

//------------------------------------------------------------------------------
inline auto MeshLoader::ref_geometry(Mesh& mesh) -> std::vector<double>& {
  return const_cast<std::vector<double>&>(mesh.geometry());
}

inline auto MeshLoader::ref_topology(Mesh& mesh, int d0, int d1)
    -> MeshConnectivity& {
  return const_cast<MeshConnectivity&>(mesh.topology(d0, d1));
}

}  // namespace wave