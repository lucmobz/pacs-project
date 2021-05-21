#include <Eigen/Core>
#include <tuple>

namespace wave {
namespace {
/**
  This function is the implementation of compute_map. It uses variadic templates
  and the index sequence, int... variadic trick to iterate at compile time.
*/
template <typename Scalar, int NROWS, int... OPTIONS,
          template <typename, int, int...> typename... Points,
          std::size_t... Is>
constexpr auto compute_map_impl(
    const std::tuple<Points<Scalar, NROWS, OPTIONS...>...>& points,
    std::index_sequence<Is...>)
    -> Eigen::Matrix<Scalar, NROWS, sizeof...(Points) - 1, Eigen::RowMajor> {
  // the implementation expects points in a tuple and an index sequence from
  // which Is is deduced as a pack of integers. construct the jacobian matrix by
  // deducing the scalar type the number of rows and the number of columns from
  // the set of points
  Eigen::Matrix<Scalar, NROWS, sizeof...(Points) - 1> map{};
  // iterate using a fold expression and the index sequence, int... variadic
  // trick. Each column will be the difference of two consecutive points.
  ((map.col(Is) = (std::get<Is + 1>(points) - std::get<0>(points))), ...);
  // now the matrix will be like [x1 - x0 | x2 - x0 | x3 - x0]
  return map;
}
}  // namespace
/** This function computes the jacobian of the affine transformation that maps
  the reference simplex into the physical simplex, an aribitraty number of
  points of arbitrary dimension can be passed. If 4 3d points are passed the
  map refers to a tetrahedron in 3d, if 2 3d points are passed the map takes a
  2d triangle into a 3d triangle (one that lives in 3d), if 3 2d points are
  passed, the map refers to a triangle in 2d, and so on.
*/
template <typename Scalar, int NROWS, int... OPTIONS,
          template <typename, int, int...> typename... Points>
constexpr auto compute_map(const Points<Scalar, NROWS, OPTIONS...>&... points)
    -> Eigen::Matrix<Scalar, NROWS, sizeof...(Points) - 1, Eigen::RowMajor> {
  // This functions takes a variadic number of points (mainly Eigen points),
  // each point is a template, the scalar type is deduced and also the number of
  // rows (options are variadic to allow points with more template parameters),
  // this way a matrix is outputted with the same number of rows; in order to
  // set the rows of the matrix as the difference of two consecutive points
  // (this is how the affine map is built mathematically) an index sequence is
  // passed to iterate at compile time through a fold expression. Therefore the
  // implementation is called by packing the points into a tuple and passing an
  // index sequence of the proper size, from which the integer parameter pack
  // can be deduced.
  return compute_map_impl(std::make_tuple(points...),
                          std::make_index_sequence<(sizeof...(Points)) - 1>{});
}

}  // namespace wave