#include "Face.hpp"
#include <Eigen/Core>
#include "Node.hpp"

namespace wave {

auto Face::normal() const -> point_t {
  // get the first two nodes in the face
  NodeIterator n0(*this, 0);
  NodeIterator n1(*this, 1);

  // get the coordinates and barycenter of the faces
  const auto& x0 = n0->point();
  const auto& x1 = n1->point();
  const auto& barycenter = this->barycenter();

  // compute normal
  return (x0 - barycenter).cross(x1 - barycenter).normalized();
}

auto Face::normal(const point_t& orientation) const -> point_t {
  const auto& normal = this->normal();
  // get a node that is on the plane of the face to decide if the normal is
  // outward or inward, of course the orientation point doesn't have to be
  // coplanar, it's best to give the barycenter of a convex element.
  NodeIterator n(*this, 0);

  // flip normal to match the orientation point direction
  if (normal.dot(orientation - n->point()) < 0)
    return normal;
  else
    return -normal;
}

}  // namespace wave