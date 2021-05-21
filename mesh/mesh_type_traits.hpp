#pragma once

#include <type_traits>
#include "Element.hpp"
#include "Face.hpp"
#include "Node.hpp"
#include "Traits.hpp"

namespace wave {

/// type trait to deduce the dimension from entity type
template <typename T>
inline constexpr int entity_dim;

template <>
inline constexpr int entity_dim<Node> = 0;

template <>
inline constexpr int entity_dim<Face> = Traits::DIM - 1;

template <>
inline constexpr int entity_dim<Element> = Traits::DIM;

}  // namespace wave