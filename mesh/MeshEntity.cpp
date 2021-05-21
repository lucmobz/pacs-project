#include "MeshEntity.hpp"
#include <algorithm>
#include <iterator>
#include <span>
#include <utility>
#include "Mesh.hpp"
#include "Node.hpp"
#include "Traits.hpp"

namespace wave {

// definition of destructor is needed as it's pure virtual
MeshEntity::~MeshEntity() {}

auto MeshEntity::index(const MeshEntity& entity) const -> int {
  // get the connectivity of the specific entity as a view
  std::span<const int> entities{_mesh.topology(_dim, entity._dim)(_index)};
  // find the index of the current entity inside the connectivity of the
  // super-entity.
  auto it{std::find(entities.begin(), entities.end(), entity._index)};
  if (it != end(entities))
    // return the local index of the current entity, if found.
    return std::distance(entities.begin(), it);
  else
    // return an invalid index if the entity is not found.
    return -1;
}

auto MeshEntity::diameter() const -> scalar_t {
  if (_dim == 0)
    return 0;
  else {
    // loop over all the nodes while keeping one fixed, as the diameter for a
    // polyhedron must be one of the main diagonals. Then switch to another
    // node, and only look on the next nodes in order to recompute already
    // computed diagonals.
    scalar_t diameter = 0;

    for (NodeIterator n0(*this); !n0.end(); ++n0) {
      // by summing one to the iterator the next one is obtained, (or the end
      // iterator if no more are there).
      for (NodeIterator n1 = n0 + 1; !n1.end(); ++n1) {
        const auto& x0 = n0->point();
        const auto& x1 = n1->point();

        scalar_t current_diameter = (x1 - x0).norm();
        if (current_diameter > diameter) diameter = current_diameter;
      }
    }
    return diameter;
  }
}

auto MeshEntity::barycenter() const -> point_t {
  if (_dim == 0)
    return point_t{_mesh.geometry(_index).data()};
  else {
    // compute the mean coordinates of all the nodes
    point_t barycenter{0.0, 0.0, 0.0};

    for (NodeIterator n(*this); !n.end(); ++n) {
      const auto& x = n->point();
      barycenter += x;
    }
    return barycenter / nnodes();
  }
}

auto MeshEntity::bounds() const -> std::pair<point_t, point_t> {
  NodeIterator n(*this);

  point_t xmin = n->point();
  point_t xmax = xmin;

  // compute the minimal and maximal coordinates of all the nodes
  for (n += 1; !n.end(); ++n) {
    const auto& x = n->point();

    for (int i = 0; i < mesh().dim(); ++i) {
      if (x[i] < xmin[i])
        xmin[i] = x[i];
      else if (x[i] > xmax[i])
        xmax[i] = x[i];
    }
  }
  return std::make_pair(xmin, xmax);
}

}  // namespace wave