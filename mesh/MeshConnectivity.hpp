#pragma once

#include <span>
#include <vector>

namespace wave {
/**
  \brief A MeshConnectivity represents the connectivity relations among
  MeshEntity's.

  Given a set of mesh entities of topological dimension d0 (for example,
  elements with d0 = 3), and another set of topological dimensions d1 (faces
  with d1 = 2, or nodes with d1 = 0), this class stores, for each entity the
  indices of the sub-entities that are connected to such entity.
*/
class MeshConnectivity {
 public:
  /**
    \brief Constructor that takes two topological dimensions.
    \param d0 First topological dimension.
    \param d1 Second topological dimension.
  */
  MeshConnectivity(int d0, int d1);
  /**
    \brief Push an entity connectivity into the global connectivity.
    \param connectivity The indices of the sub-entities of second topological
    dimension connected to the parent entity.
  */
  void push(const std::vector<int>& connectivity);
  /// \brief Test if the MeshConnectivity doesn't contain any entity.
  auto empty() const -> bool;
  /// \brief Get the total number of entities of first topological dimension.
  auto size() const -> int;
  /// \brief Get the total number of sub-entities of an entity.
  /// \param index The index of the entity of first topological dimension.
  auto size(int index) const -> int;
  /**
    \brief Get the indices of the sub-entities of a given entity

    The call operator returns a span, which is a view of a sequential container
    in order to avoid copies (it is a wrapper made by a pointer and a size).
    \param index The index of the entity of first topological dimension.
    \return The indices of the sub-entities connected to the specified entity.
  */
  auto operator()(int index) const -> std::span<const int>;

 private:
  /// The first topological dimension.
  int _d0{0};
  /// The second topological dimension.
  int _d1{0};
  /// The global connectivity vector, connectivities are contiguous and use an
  /// offset vector to discriminate among entities.
  std::vector<int> _connectivity{};
  /// The offset vector into the global connectivity. It starts at 0 and ends
  /// with the size of the last entity.
  std::vector<int> _offsets{};
};

//------------------------------------------------------------------------------
// _offsets is init to 0 as the first entity always starts at zero, its size
// will be the number of entities + 1.
inline MeshConnectivity::MeshConnectivity(int d0, int d1)
    : _d0{d0}, _d1{d1}, _offsets{0} {}

inline void MeshConnectivity::push(const std::vector<int>& connectivity) {
  _connectivity.insert(_connectivity.end(), connectivity.begin(),
                       connectivity.end());
  _offsets.emplace_back(_offsets.back() + connectivity.size());
}

inline auto MeshConnectivity::size() const -> int {
  return _offsets.size() - 1;
}

inline auto MeshConnectivity::size(int index) const -> int {
  return _offsets[index + 1] - _offsets[index];
}

inline auto MeshConnectivity::empty() const -> bool {
  return _connectivity.empty();
}
// A std::span is returned to avoid copies, it is a view of a sequential
// container.
inline auto MeshConnectivity::operator()(int index) const
    -> std::span<const int> {
  return {_connectivity.begin() + _offsets[index],
          _connectivity.begin() + _offsets[index] + size(index)};
}

}  // namespace wave