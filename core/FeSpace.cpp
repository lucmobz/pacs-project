#include "FeSpace.hpp"
#include <cmath>
#include "FeElement.hpp"
#include "FeFace.hpp"
#include "mesh_include.hpp"

namespace wave {

FeSpace::FeSpace(const Mesh &mesh, int degree)
    : FeSpace(mesh, std::vector<int>(mesh.nelements(), degree)) {}

FeSpace::FeSpace(const Mesh &mesh, const std::vector<int> &degrees)
    : _mesh(mesh) {
  _elements.reserve(mesh.nelements());
  _faces.reserve(mesh.nfaces());

  // create elements
  for (ElementIterator e(mesh); !e.end(); ++e) {
    // create a finite element from an element and a degree
    _elements.emplace_back(*e, degrees[e->index()]);
    // dof indices are just a sequence in this algorithm as no continuity is
    // required
    _ndof += _elements.back().ndof();
    // compute algorithm specific info for the element
    _elements.back().triangulate();
    _elements.back().diameter(e->diameter());
  }

  // create faces
  for (FaceIterator f(mesh); !f.end(); ++f) {
    // get the first element
    ElementIterator e(*f);
    const auto &fe = (*this)(*e);

    // compute penalty parameter
    scalar_t penalty = fe.degree() * fe.degree() / fe.diameter();

    if (f->is_interior()) {
      // get the second element
      ElementIterator en(*e, e->index(*f));
      const auto &fen = (*this)(*en);
      penalty = std::max(penalty, fen.degree() * fen.degree() / fen.diameter());
    }

    // compute algortihm specific information for faces
    _faces.emplace_back(*f);
    _faces.back().penalty(penalty);
    _faces.back().triangulate();
    // each normal will be exterior with respect to the first element adjacent
    // to the face
    _faces.back().normal(f->normal(e->barycenter()));
  }
}

FeSpace::dense_t FeSpace::nodal(const dense_t &dof_values) const {
  dense_t nodal_values;
  // Here nodes are counted as they appear in each cell (so they appear many
  // times), using a MeshFunction is helpful and also covers the case where the
  // underlying mesh doesn't have its nodes duplicated. So meshes with multiple
  // definition of nodes are treated as meshes with nodes appearing only once,
  // in particular mesh.nnodes() in that case will return the number of nodes,
  // and here we need a different thing, as a node is counted as many times as
  // it appears in the cells.
  MeshFunction nnodes(mesh(),
                      [](const Element &element) { return element.nnodes(); });
  nodal_values.resize(std::accumulate(nnodes.begin(), nnodes.end(), 0));

  // loop over the elements
  int index = 0;
  for (ElementIterator e(mesh()); !e.end(); ++e) {
    // get the finite element
    const auto &fe = (*this)(*e);

    // loop over the nodes of the element
    for (NodeIterator n(*e); !n.end(); ++n) {
      const auto &x = n->point();
      scalar_t value = 0.0;

      // fromt he dof of the solution compute the value of the solution at the
      // nodes
      for (int i = 0; i < fe.ndof(); ++i) {
        value += fe.eval_basis(i, fe.map(x)) * dof_values[fe.dof(i)];
      }
      nodal_values[index++] = value;
    }
  }

  return nodal_values;
}

}  // namespace wave
