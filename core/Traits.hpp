#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include <limits>
#include <type_traits>

namespace wave {
/// \brief This namespace collects the main types and constants used in the
/// problem.
namespace Traits {
/// \brief Scalar type used in the problem can be set here
using scalar_t = double;
/// \brief Tolerance for excluding tetrahedral sleeves from triangulations.
/// Tetrahedrons with a volume smaller than the polyhedron volume times the
/// tolerance are excluded.
inline constexpr scalar_t EPSCGAL = 1e-15;
/// \brief Dimensionality of the problem
inline constexpr int DIM = 3;
/// \brief Point type used in the problem (as many Eigen methods are used,
/// better to leave at default)
using point_t = Eigen::Vector3d;
/// \brief Proxy for a point, used to avoid copying coordinates when not needed.
using mpoint_t = Eigen::Map<const point_t>;
/// \brief Quadrature point type, it is templated on the number of coordinates
/// as both 2d and 3d points are needed. \tparam NROWS number of coordinates of
/// the point
template <int NROWS>
using qpoint_t = Eigen::Matrix<scalar_t, NROWS, 1>;
/// \brief Small dense matrix type to represent affine transformations from
/// tetrahedrons or triangles.
/// \tparam NCOLS number of columns in the matrix
template <int NCOLS>
using matrix_t = Eigen::Matrix<double, DIM, NCOLS, Eigen::RowMajor>;
/// \brief Triplet type to fill sparse matrices.
using triplet_t = Eigen::Triplet<scalar_t>;
/// \brief Sparse matrix type.
using sparse_t = Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>;
/// \brief Dense vector type used for right hand sides in linear systems.
using dense_t = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
/// \brief Format used to print Eigen matrices, useful for debugging.
inline const Eigen::IOFormat IOFormat(Eigen::StreamPrecision, 0, " ", "\n", "[",
                                      "]");
}  // namespace Traits
}  // namespace wave