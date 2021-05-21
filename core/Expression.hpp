#pragma once

#include <tuple>
#include <type_traits>
#include "Traits.hpp"
#include "expression_type_traits.hpp"

namespace wave {

/// \brief Dummy base class to be used in type traits, to check if something
/// derives from it means it is an Expression.
struct ExpressionBase {};

/**
  \brief Core class of the Expression template system.

  Represents a generic mathematical or integral expression (unary, binary,
  etc...), it is templated on an operator and a set of arguments. Expression's
  are built by suitable operator overloading. They store the concrete operation
  to do and the arguments, and thanks to inlining in the actual code there will
  be the full Expression tree unfolded.

  \tparam Op The operator object, can be any lambda or functor.
  \tparam Args The arguments of the expression.
*/
template <typename Op, typename... Args>
// All expressions derive from a base class (used in type traits)
class Expression : public ExpressionBase {
 public:
  /// \brief Constructs an Expression from a lambda or functor and a set of
  /// arguments.
  explicit Expression(const Op &op, const Args &... args)
      : _op(op), _args(args...) {}

  /**
    \brief Evaluate an Expression on a set of parameters (Analog of evaluating a
    vector expression at a particular subscript)

    This is a generic method to handle all possible contexts in which an
    Expression is evaluated. For the main problem in the library, there is a
    need to evaluate Expression's with different sets of parameters (for example
    boundary integrals Expression's have different parameters than volume
    integrals Expression's).

    \tparam Params Parameter pack that defines the generic list of arguments
    where an Expression can be evaluated.
  */
  template <typename... Params>
  auto eval(const Params &... params) const {
    // Apply the operator to the set of arguments in the expression, but first
    // subscript them (i.e. check if any of them is a scalar, if not reapply
    // eval() to the arguments in order to move through the Expression tree).
    auto eval_at = [this, &params...](const auto &... args) {
      return _op(_subscript(args, params...)...);
    };
    // apply function that acts on several arguments by first unpacking a tuple
    // into the set of arguments
    return std::apply(eval_at, _args);
  }

 private:
  /// \brief Store a reference to the operator.
  const Op &_op;
  /// \brief Store arguments in a tuple of constant references.
  std::tuple<const Args &...> _args;

  /**
    \brief Discriminate between Expression types and scalars or Eigen types.

    Check if arg is an Expression, if yes, evaluate it, if not, return it as is.
  */
  template <typename Arg, typename... Params>
  static constexpr auto _subscript(const Arg &arg, const Params &... params) {
    if constexpr (is_Expression_or_Terminal_v<Arg>)
      // if the argument is an expression recurse into an evaluation.
      return arg.eval(params...);
    else
      // if the argument is a scalar or Eigen type just return it
      return arg;
  }
};

}  // namespace wave