#pragma once

#include <Eigen/Core>
#include <type_traits>

namespace wave {
// This file contains type traits for the expression system, it is used to check
// if a type is in the expression system, and also to parse an expression at
// compile time to figure out if basis functions or basis gradients have to be
// computed (for example if the expression is int3d<>(u*v) then the gradients
// will not be computed by the int3d<> method)
// ----------------------------------------------------------------
// forward decl are useful in type traits
struct TerminalBase;
struct ExpressionBase;

template <typename Op, typename... Args>
class Expression;

template <typename... Functions>
class Field;

class Penalty;
class Normal;

class Trial;
class Test;

struct Grad;
struct Avg;
struct Jump;

using GradTrial = Expression<Grad, Trial>;
using GradTest = Expression<Grad, Test>;

template <>
class Expression<Grad, Trial>;
template <>
class Expression<Grad, Test>;
template <>
class Expression<Avg, GradTrial>;
template <>
class Expression<Avg, GradTest>;
template <>
class Expression<Jump, Trial>;
template <>
class Expression<Jump, Test>;

using AvgGradTrial = Expression<Avg, GradTrial>;
using AvgGradTest = Expression<Avg, GradTest>;
using JumpTrial = Expression<Jump, Trial>;
using JumpTest = Expression<Jump, Test>;

// -----------------------------------------------------------------------------
// type traits to test if a type is in the expression system
template <typename T>
inline constexpr bool is_Expression_v =
    std::is_base_of_v<ExpressionBase, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_Terminal_v =
    std::is_base_of_v<TerminalBase, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_Expression_or_Terminal_v =
    is_Expression_v<std::remove_cvref_t<T> > ||
    is_Terminal_v<std::remove_cvref_t<T> >;

// -----------------------------------------------------------------------------
// type traits to identify types in the expression system
template <typename T>
inline constexpr bool is_Eigen_v =
    std::is_base_of_v<Eigen::MatrixBase<std::remove_cvref_t<T> >,
                      std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_Trial_v =
    std::is_same_v<Trial, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_Test_v = std::is_same_v<Test, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_GradTrial_v =
    std::is_same_v<GradTrial, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_GradTest_v =
    std::is_same_v<GradTest, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_AvgGradTrial_v =
    std::is_same_v<AvgGradTrial, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_AvgGradTest_v =
    std::is_same_v<AvgGradTest, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_JumpTrial_v =
    std::is_same_v<JumpTrial, std::remove_cvref_t<T> >;

template <typename T>
inline constexpr bool is_JumpTest_v =
    std::is_same_v<JumpTest, std::remove_cvref_t<T> >;

//------------------------------------------------------------------------------
// type traits to parse the expression types
// this contains the impl of the has_Trial_or_Test_v type trait (see below)
namespace {
template <typename T>
inline constexpr bool has_Trial_or_Test_v_impl = false;

template <>
inline constexpr bool has_Trial_or_Test_v_impl<Trial> = true;

template <>
inline constexpr bool has_Trial_or_Test_v_impl<Test> = true;

// use fold expressions to parse the expression
template <typename Op, typename... Args>
inline constexpr bool has_Trial_or_Test_v_impl<Expression<Op, Args...> > =
    (has_Trial_or_Test_v_impl<Args> || ...) &&
    !is_GradTrial_v<Expression<Op, Args...> > &&
    !is_GradTest_v<Expression<Op, Args...> >;
}  // namespace

// type trait to check if an expression contains a Trial or a Test
template <typename T>
inline constexpr bool has_Trial_or_Test_v =
    has_Trial_or_Test_v_impl<std::remove_cvref_t<T> >;

// this contains the impl of the has_GradTrial_or_GradTest_v type trait (see
// below)
namespace {
template <typename T>
inline constexpr bool has_GradTrial_or_GradTest_v_impl = false;

template <>
inline constexpr bool has_GradTrial_or_GradTest_v_impl<GradTrial> = true;

template <>
inline constexpr bool has_GradTrial_or_GradTest_v_impl<GradTest> = true;

// use fold expressions to parse the expression
template <typename Op, typename... Args>
inline constexpr bool
    has_GradTrial_or_GradTest_v_impl<Expression<Op, Args...> > =
        (has_GradTrial_or_GradTest_v_impl<Args> || ...);
}  // namespace

// type trait to check if an expression contains a GradTrial or a GradTest
template <typename T>
inline constexpr bool has_GradTrial_or_GradTest_v =
    has_GradTrial_or_GradTest_v_impl<std::remove_cvref_t<T> >;

}  // namespace wave