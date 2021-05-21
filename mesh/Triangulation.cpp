#ifndef CGAL_NO_GMP
#define CGAL_NO_GMP 1
#endif
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <Eigen/Core>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>
#include "MeshEntity.hpp"
#include "Node.hpp"
#include "Traits.hpp"
#include "Triangulation.hpp"
#include "geometry.hpp"

namespace wave {

/**
  \brief Specialization of the Triangulation constructor for a Face.

  A triangulation is built by creating (using CGAL) a constrained Delaunay
  triangulation of a planar polygon embedded in 3d. The jacobians are assembled
  using the compute_map function.
*/
template <>
Triangulation<Face>::Triangulation(const Face& face) {
  // cgal aliases
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Traits = CGAL::Triangulation_2_projection_traits_3<Kernel>;
  using Delaunay = CGAL::Constrained_Delaunay_triangulation_2<Traits>;
  using Point = Kernel::Point_3;
  using Vector = Kernel::Vector_3;

  // compute points, normal and element barycenter
  int npoints = face.nnodes();

  std::vector<Point> points;
  points.reserve(npoints);

  for (NodeIterator n(face); !n.end(); ++n) {
    const auto& x = *n;
    points.emplace_back(x[0], x[1], x[2]);
  }

  ElementIterator e(face);
  auto normal = face.normal(e->barycenter());

  // compute constraints
  std::vector<std::pair<int, int>> constraints;
  constraints.reserve(npoints);

  for (int i = 0; i < npoints - 1; ++i) constraints.emplace_back(i, i + 1);
  constraints.emplace_back(npoints - 1, 0);

  // compute cgal constrained triangulation, it needs a normal vector
  Vector plane_orthogonal_vector(normal[0], normal[1], normal[2]);
  Traits traits(plane_orthogonal_vector);

  Delaunay dt(traits);
  dt.insert_constraints(points.begin(), points.end(), constraints.begin(),
                        constraints.end());

  // compute data
  for (Delaunay::Face_handle f : dt.finite_face_handles()) {
    const auto& x0 = f->vertex(0)->point();
    const auto& x1 = f->vertex(1)->point();
    const auto& x2 = f->vertex(2)->point();

    _maps.emplace_back(
        compute_map(point_t(x0[0], x0[1], x0[2]), point_t(x1[0], x1[1], x1[2]),
                    point_t(x2[0], x2[1], x2[2])),
        point_t(x0[0], x0[1], x0[2]), std::sqrt(dt.triangle(f).squared_area()));
  }
}

/**
  \brief Specialization of a Triangulation of an Element (which is a
  polyhedron).

  A triangulation is built by using CGAL Delaunay triangulation methods. Sleeves
  simplices are removed if their volume is smaller than a tolerance times the
  element volume.  The jacobians are assembled using the compute_map function.
*/
template <>
Triangulation<Element>::Triangulation(const Element& element) {
  // cgal aliases
  using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Delaunay = CGAL::Delaunay_triangulation_3<Kernel>;
  using Point = Kernel::Point_3;

  // compute points
  int npoints = element.nnodes();

  std::vector<Point> points;
  points.reserve(npoints);

  for (NodeIterator n(element); !n.end(); ++n) {
    const auto& x = *n;
    points.emplace_back(x[0], x[1], x[2]);
  }

  // compute cgal triangulation
  Delaunay dt(points.begin(), points.end());

  // compute data
  std::vector<scalar_t> measures;
  measures.reserve(dt.number_of_finite_cells());

  // compute total volume and the minimal allowed volume
  for (Delaunay::Cell_handle c : dt.finite_cell_handles())
    measures.emplace_back(dt.tetrahedron(c).volume());
  scalar_t min_measure = std::accumulate(measures.begin(), measures.end(),
                                         static_cast<scalar_t>(0)) *
                         EPSCGAL;

  // build the triangulation by computing the maps
  int i = 0;
  for (Delaunay::Cell_handle c : dt.finite_cell_handles()) {
    if (measures[i] > min_measure) {
      const auto& x0 = c->vertex(0)->point();
      const auto& x1 = c->vertex(1)->point();
      const auto& x2 = c->vertex(2)->point();
      const auto& x3 = c->vertex(3)->point();

      _maps.emplace_back(compute_map(point_t(x0[0], x0[1], x0[2]),
                                     point_t(x1[0], x1[1], x1[2]),
                                     point_t(x2[0], x2[1], x2[2]),
                                     point_t(x3[0], x3[1], x3[2])),
                         point_t(x0[0], x0[1], x0[2]), measures[i]);
    }
    ++i;
  }
}

}  // namespace wave