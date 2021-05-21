#pragma once

#include <cassert>
#include <utility>
#include "Traits.hpp"
#include "mesh_include.hpp"
#include "utilities_include.hpp"

namespace wave {
/**
  \brief A FeFace represents the face between two finite elements in a FeSpace.

  This class is built on top of a mesh Face (which is basically an index with
  methods to compute associated information), and stores all the algorithm
  specific information associated with the Face.
*/
class FeFace {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;

  /// \brief Nested class that represents the side of a face (used in expressions to
  /// determine the orientation of the normal)
  // an enum class slows down the code considerably, but a multiplication by +1
  // or -1 is faster
  struct Side {
    static constexpr int POSITIVE = 1;
    static constexpr int NEGATIVE = -1;
  };
  /// \brief Construct a FeFace from a mesh Face.
  FeFace(const Face &face);
  /// \brief Tests if the underlying Face is interior.
  bool is_interior() const;
  /// \brief Tests if the underlying Face is exterior.
  bool is_exterior() const;
  /// \brief Get the boundary tag of the FeFace (0 interior, greater than 0 one
  /// of the boundaries)
  int boundary() const;
  /// \brief Get the DG penalty parameter associated with the FeFace.
  scalar_t penalty() const;
  /// \brief Set the DG penalty parameter.
  void penalty(scalar_t penalty);
  /// \brief Get the normal to the Face.
  const point_t &normal() const;
  /// \brief Set the normal to the Face (normally by using the Face specific
  /// method that also allows to set the orientation)
  void normal(const point_t &normal);
  /// \brief Get the Triangulation of the Face used in quadrature rules.
  const Triangulation<Face> &triangulation() const;
  /// \brief Compute the Triangulation of the Face.
  void triangulate();

 private:
  /// Underlying mesh Face.
  Face _face;
  /// Boundary of the face, 0 is interior, 1 2 3... are exterior
  int _boundary{0};
  /// Value of the DG penalty.
  scalar_t _penalty{0.0};
  /// Value of the normal vector.
  point_t _normal{0.0, 0.0, 0.0};
  /// Triangulation of the underlying Face.
  Triangulation<Face> _triangulation;
};

//------------------------------------------------------------------------------
inline FeFace::FeFace(const Face &face) : _face(face) {
  if (face.is_exterior()) _boundary = face.boundary();
}

inline bool FeFace::is_interior() const { return _boundary == 0; }

inline bool FeFace::is_exterior() const { return !is_interior(); }

inline int FeFace::boundary() const { return _boundary; }

inline void FeFace::penalty(scalar_t penalty) { _penalty = penalty; }

inline FeFace::scalar_t FeFace::penalty() const { return _penalty; }

inline const FeFace::point_t &FeFace::normal() const { return _normal; }

inline void FeFace::normal(const point_t &normal) { _normal = normal; }

inline void FeFace::triangulate() {
  _triangulation = std::move(Triangulation<Face>(_face));
}

inline const Triangulation<Face> &FeFace::triangulation() const {
  assert(!_triangulation.empty());
  return _triangulation;
}

}  // namespace wave