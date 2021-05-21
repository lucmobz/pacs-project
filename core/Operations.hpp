#pragma once

#include <Eigen/Core>
#include <tuple>
#include <type_traits>
#include "Expression.hpp"
#include "expression_type_traits.hpp"

namespace wave {

// sum
template <typename Lhs, typename Rhs>
constexpr bool is_sum_valid =
    is_Expression_or_Terminal_v<Lhs> || is_Expression_or_Terminal_v<Rhs>;

/// \brief Dummy operator that builds a sum expression
template <typename Lhs, typename Rhs,
          typename = std::enable_if_t<is_sum_valid<Lhs, Rhs> > >
inline auto operator+(const Lhs &lhs, const Rhs &rhs) {
  // return an Expression by passing the lambda concrete operator and the two
  // arguments
  return Expression{[](auto l, auto r) { return l + r; }, lhs, rhs};
}

// product
template <typename Lhs, typename Rhs>
constexpr bool is_product_valid =
    is_Expression_or_Terminal_v<Lhs> || is_Expression_or_Terminal_v<Rhs>;

/// \brief Dummy operator that builds a product expression (or matrix vector
/// product)
template <typename Lhs, typename Rhs,
          typename = std::enable_if_t<is_product_valid<Lhs, Rhs> > >
inline auto operator*(const Lhs &lhs, const Rhs &rhs) {
  // The trailing return type is necessary to distinguish the case where an
  // Eigen matrix is given in an expression (for example a diffusion matrix), or
  // any other product operation with Eigen types because Eigen uses its own
  // expression template system The return type of the lambda will be deduced
  // from the left and right operand (this is tested only through one level of
  // recursion) so it shouldn't be an Eigen expression.
  return Expression{
      [](auto l, auto r)
          -> std::conditional<
              is_Eigen_v<decltype(l)> || is_Eigen_v<decltype(r)>,
              typename std::conditional<is_Eigen_v<decltype(l)>, decltype(l),
                                        decltype(r)>::type,
              decltype(l * r)>::type{return l * r;
}
, lhs, rhs
};  // namespace wave
}

// these overloads are to handle Eigen expression templates (old version)
/*
template <typename Lhs,
          typename = std::enable_if_t<is_product_valid<Lhs, GradTrial>>>
inline auto operator*(const Lhs& lhs, const GradTrial& rhs) {
  return Expression{[](auto l, auto r) -> decltype(r) { return l * r; }, lhs,
rhs};
}

template <typename Rhs,
          typename = std::enable_if_t<is_product_valid<GradTrial, Rhs>>>
inline auto operator*(const GradTrial& lhs, const Rhs& rhs) {
  return Expression{[](auto l, auto r) -> decltype(l) { return l * r; }, lhs,
rhs};
}
*/

// difference
template <typename Lhs, typename Rhs>
constexpr bool is_difference_valid =
    is_Expression_or_Terminal_v<Lhs> || is_Expression_or_Terminal_v<Rhs>;

/// \brief Dummy operator to build a difference expression
template <typename Lhs, typename Rhs,
          typename = std::enable_if_t<is_difference_valid<Lhs, Rhs> > >
inline auto operator-(const Lhs &lhs, const Rhs &rhs) {
  return Expression{[](auto l, auto r) { return l - r; }, lhs, rhs};
}

// negation
template <typename Rhs>
constexpr bool is_negation_valid = is_Expression_or_Terminal_v<Rhs>;

/// \brief Dummy operator to build a negation unary expression
template <typename Rhs, typename = std::enable_if_t<is_negation_valid<Rhs> > >
inline auto operator-(const Rhs &rhs) {
  return Expression{[](auto r) { return -r; }, rhs};
}

// dot
template <typename Lhs, typename Rhs>
constexpr bool is_inner_product_valid =
    is_Expression_or_Terminal_v<Lhs> || is_Expression_or_Terminal_v<Rhs>;

/// \brief Dummy operator to build a dot product expression
template <typename Lhs, typename Rhs,
          typename = std::enable_if_t<is_inner_product_valid<Lhs, Rhs> > >
inline auto dot(const Lhs &lhs, const Rhs &rhs) {
  // here Eigen types don't create problems as the dot product is evaluated
  // immediately
  return Expression{[](const auto &l, const auto &r) { return l.dot(r); }, lhs,
                    rhs};
}

// gradient
template <typename Rhs>
constexpr bool is_gradient_valid = is_Trial_v<Rhs> || is_Test_v<Rhs>;

/// \brief Dummy class used in the grad dummy operator
struct Grad {};

/// \brief Dummy operator to build a gradient expression from a Trial or Test
/// object.
template <typename Rhs, typename = std::enable_if_t<is_gradient_valid<Rhs> > >
inline auto grad(const Rhs &rhs) {
  return Expression{Grad(), rhs};
}

// average
template <typename Rhs>
constexpr bool is_average_valid = is_GradTrial_v<Rhs> || is_GradTest_v<Rhs>;

/// \brief Dummy class used in the avg dummy operator
struct Avg {};
/// \brief Dummy operator to build an average of a gradient of a Trial or Test
/// expression
template <typename Rhs, typename = std::enable_if_t<is_average_valid<Rhs> > >
inline auto avg(const Rhs &rhs) {
  return Expression{Avg(), rhs};
}

// jump
template <typename Rhs>
constexpr bool is_jump_valid = is_Trial_v<Rhs> || is_Test_v<Rhs>;

/// \brief Dummy class used in the jump dummy operator
struct Jump {};
/// \brief Dummy operator to build a jump of a Trial or Test expression
template <typename Rhs, typename = std::enable_if_t<is_jump_valid<Rhs> > >
inline auto jump(const Rhs &rhs) {
  return Expression{Jump(), rhs};
}
}