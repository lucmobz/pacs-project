#pragma once

#include <tuple>
#include <type_traits>
#include "Expression.hpp"
#include "FeFace.hpp"
#include "Traits.hpp"
#include "functional_include.hpp"

namespace wave {

/// \brief Dummy class to check if a type is the terminal leaf of an Expression
/// tree.
struct TerminalBase {};

// field
/**
  \brief Main function type in the Expression system.

  A Field represents a function or vector of functions that acts on a point in
  3d. If it's constructed by passing a single lambda it's a scalar function, if
  by passing 3 lambdas it is a vector function and will return a vector. When an
  Expression involves a function then a Field must be used. Just create a Field
  by passing some lambdas and use it in an Expression.

  \tparam Functions The functions passed are deduced in the constructor.
*/
template <typename... Functions>
class Field : public TerminalBase {
 private:
  static constexpr int NFUNCTIONS = sizeof...(Functions);

 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;
  template <int DIM>
  using qpoint_t = Traits::qpoint_t<DIM>;
  using return_t =
      std::conditional<NFUNCTIONS == 1, scalar_t, qpoint_t<NFUNCTIONS> >::type;

  /// \brief Constructs a generic function by taking a set of lambdas (the
  /// components of the function).
  Field(const Functions &... functions)
      : _functions(std::make_tuple(functions...)) {}

  /// \brief Evaluate the scalar or vector function at a specific point.
  return_t operator()(const point_t &x) const {
    auto eval_at = [&x](const auto &... functions) {
      return return_t{functions(x)...};
    };
    return std::apply(eval_at, _functions);
  }
  /// \brief Evaluate the Field as a terminal leaf in an expression tree.
  template <typename... Params>
  return_t eval(const Params &... params) const {
    // any expression evaluation context involving a Field will have as the last
    // parameter the point where to evaluate the field
    return (*this)(std::get<sizeof...(Params) - 1>(std::tie(params...)));
  }

 private:
  /// The underlying functions used by the Field are stored in a tuple.
  std::tuple<Functions...> _functions;
};

// penalty
/// \brief A Penalty represents the DG method penalty function in the Expression
/// system. Just define a Penalty variable with a value and use it in an
/// Expression.
class Penalty : public TerminalBase {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;

  /// \brief Construct a Penalty with a given value, the final penalty function
  /// will be the classic DG method penalty: max(c * p^2 / h) between the two
  /// faces.
  Penalty(scalar_t penalty) : _penalty(penalty) {}
  /// \brief Evaluate the penalty on a finite element face.
  scalar_t operator()(const FeFace &fef) const {
    return _penalty * fef.penalty();
  }
  /// \brief Evaluate the Penalty as a terminal leaf in an expression tree.
  template <typename... Params>
  scalar_t eval(const Params &... params) const {
    // any expression evaluating context involving a Penalty will have as first
    // parameter the finite element face from which to extract the values of p
    // and h.
    const auto &arg0 = std::get<0>(std::tie(params...));
    return (*this)(arg0);
  }

 private:
  scalar_t _penalty{0.0};
};

// normal
/// \brief A Normal represents the outward normal vector to FeElement's and
/// FeFace's in the Expression system. Just declare a Normal variable and use it
/// in an Expression.
class Normal : public TerminalBase {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;

  /// \brief Evalate the normal on a FeFace and turn it to the given side
  /// \param fef The FeFace to extract the normal from
  /// \param side Either +1 or -1 to set the side of the normal.
  point_t operator()(const FeFace &fef, const int &side) const {
    return fef.normal() * side;
  }
  /// \brief Evaluate the Normal as a terminal leaf in an expression tree.
  template <typename... Params>
  point_t eval(const Params &... params) const {
    // first argument in the expression context is the finite element face
    const auto &arg0 = std::get<0>(std::tie(params...));
    // second argument in the expression context is the side of the face (+-1)
    const auto &arg1 = std::get<1>(std::tie(params...));
    return (*this)(arg0, arg1);
  }
};

// trial
/// \brief It represents a Trial function in an Expression. Just declare a Trial
/// variable and use it in an Expression.
class Trial : public TerminalBase {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;

  /**
    \brief Evaluate the Trial function in a bilinear volume integration
    expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf. \param u Value of the basis function
    related to the Trial function. \param v Value of the basis function related
    to the Test function. \param du Value of the gradient of the basis function
    related to the Trial function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  scalar_t eval(const scalar_t &u, const scalar_t &v, const point_t &du,
                const point_t &dv, const point_t &x) const {
    return u;
  }
};

// test
/// \brief It represents a Test function in an Expression. Just declare a Test
/// variable and use it in an Expression.
class Test : public TerminalBase {
 public:
  using scalar_t = Traits::scalar_t;
  using point_t = Traits::point_t;

  /**
    \brief Evaluate the Test function in a bilinear volume integration
    expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf. \param u Value of the basis function
    related to the Trial function. \param v Value of the basis function related
    to the Test function. \param du Value of the gradient of the basis function
    related to the Trial function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  scalar_t eval(const scalar_t &u, const scalar_t &v, const point_t &du,
                const point_t &dv, const point_t &x) const {
    return v;
  }
  /**
    \brief Evaluate the Test function in a linear volume integration expression
    context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf. \param u Value of the basis function
    related to the Trial function. \param v Value of the basis function related
    to the Test function. \param du Value of the gradient of the basis function
    related to the Trial function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  scalar_t eval(const scalar_t &v, const point_t &dv, const point_t &x) const {
    return v;
  }
  /// \brief Evaluate the Trial function in a linear surface integration
  /// expression context.
  scalar_t eval(const FeFace &fef, const int &v_side, const scalar_t &v,
                const point_t &dv, const point_t &x) const {
    return v;
  }
};

// specialized terminals
// Dummy classes to represents a composite operations (forward decl)
struct Grad;
struct Avg;
struct Jump;

using GradTrial = Expression<Grad, Trial>;
using GradTest = Expression<Grad, Test>;

// gradient of trial
/**
  \brief This Expression represents the gradient of a Trial function.

  It is used when in an Expression the term grad(u) appears (where u is a Trial
  function). The specilization is needed to turn a particular Expression into a
  terminal leaf of the tree and perform the evaluation in a specific way
  breaking the recursion.
*/
template <>
class Expression<Grad, Trial> : public TerminalBase {
 public:
  using point_t = Traits::point_t;
  using scalar_t = Traits::scalar_t;

  /// Constructor for this specialized Expression.
  Expression(const Grad &, const Trial &) {}
  /**
    \brief Evaluate the gradient of a Trial function in a bilinear volume
    integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf. \param u Value of the basis function
    related to the Trial function. \param v Value of the basis function related
    to the Test function. \param du Value of the gradient of the basis function
    related to the Trial function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  point_t eval(const scalar_t &u, const scalar_t &v, const point_t &du,
               const point_t &dv, const point_t &x) const {
    return du;
  }
};

// gradient of test
/**
  \brief This Expression represents the gradient of a Test function.

  It is used when in an expression the term grad(v) appears (where v is a Test
  function). The specilization is needed to turn a particular Expression into a
  terminal leaf of the tree and perform the evaluation in a specific way
  breaking the recursion.
*/
template <>
class Expression<Grad, Test> : public TerminalBase {
 public:
  using point_t = Traits::point_t;
  using scalar_t = Traits::scalar_t;

  /// Constructor for this specialized Expression.
  Expression(const Grad &, const Test &) {}
  /**
    \brief Evaluate the gradient of Test function in a bilinear volume
    integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf. \param u Value of the basis function
    related to the Trial function. \param v Value of the basis function related
    to the Test function. \param du Value of the gradient of the basis function
    related to the Trial function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  point_t eval(const scalar_t &u, const scalar_t &v, const point_t &du,
               const point_t &dv, const point_t &x) const {
    return dv;
  }

  /**
    \brief Evaluate the gradient of Test function in a linear volume integration
    expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf. \param v Value of the basis function
    related to the Test function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  point_t eval(const scalar_t &v, const point_t &dv, const point_t &x) const {
    return dv;
  }

  /**
    \brief Evaluate the gradient of Test function in a linear surface
    integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf.
    \param fef FeFace used in the expression context (for example to fetch the
    penalty or normal). \param v_side The side of the FeFace to be considered
    (either 1 or -1, it will multiply the normal to get the correct result)
    \param v Value of the basis function related
    to the Test function. \param dv Value of the gradient of the basis
    function related to the Test function. \param x Physical point.
  */
  point_t eval(const FeFace &fef, const int &v_side, const scalar_t &v,
               const point_t &dv, const point_t &x) const {
    return dv;
  }
};

// average of trial gradient
/**
  \brief This Expression represents the average of the gradient of a Trial
  function.

  It is used when in an expression the term avg(grad(u)) appears (where u is a
  Trial function). The specilization is needed to turn a particular Expression
  into a terminal leaf of the tree and perform the evaluation in a specific way
  breaking the recursion.
*/
template <>
class Expression<Avg, GradTrial> : public TerminalBase {
 public:
  using point_t = Traits::point_t;
  using scalar_t = Traits::scalar_t;
  /// Constructor for this specialized Expression.
  Expression(const Avg &, const GradTrial &) {}

  /**
    \brief Evaluate the average of the gradient of a Trial function in a
    bilinear surface integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf.
    \param fef FeFace used in the expression context (for example to fetch the
    penalty or normal). \param u_side The side of the FeFace to be considered
    for the Trial function (either 1 or -1, it will multiply the normal to get
    the correct result) \param v_side The side of the FeFace to be considered
    for the Test function (either 1 or -1, it will multiply the normal to get
    the correct result) \param u Value of the basis function related to the
    Trial function. \param v Value of the basis function related to the Test
    function. \param du Value of the gradient of the basis related to the Trial
    function. \param dv Value of the gradient of the basis related to the Test
    function. \param x Physical point.
  */
  point_t eval(const FeFace &fef, const int &u_side, const int &v_side,
               const scalar_t &u, const scalar_t &v, const point_t &du,
               const point_t &dv, const point_t &x) const {
    scalar_t coefficient = (fef.is_interior() ? 0.5 : 1.0);
    return coefficient * du;
  }
};

// average of test gradient
/**
  \brief This Expression represents the average of the gradient of a Test
  function.

  It is used when in an expression the term avg(grad(v)) appears (where v is a
  Test function). The specilization is needed to turn a particular Expression
  into a terminal leaf of the tree and perform the evaluation in a specific way
  breaking the recursion.
*/
template <>
class Expression<Avg, GradTest> : public TerminalBase {
 public:
  using point_t = Traits::point_t;
  using scalar_t = Traits::scalar_t;
  /// Constructor for this specialized Expression.
  Expression(const Avg &, const GradTest &) {}

  /**
    \brief Evaluate the average of the gradient of a Test function in a bilinear
    surface integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf.
    \param fef FeFace used in the expression context (for example to fetch the
    penalty or normal). \param u_side The side of the FeFace to be considered
    for the Trial function (either 1 or -1, it will multiply the normal to get
    the correct result) \param v_side The side of the FeFace to be considered
    for the Test function (either 1 or -1, it will multiply the normal to get
    the correct result) \param u Value of the basis function related to the
    Trial function. \param v Value of the basis function related to the Test
    function. \param du Value of the gradient of the basis related to the Trial
    function. \param dv Value of the gradient of the basis related to the Test
    function. \param x Physical point.
  */
  point_t eval(const FeFace &fef, const int &u_side, const int &v_side,
               const scalar_t &u, const scalar_t &v, const point_t &du,
               const point_t &dv, const point_t &x) const {
    scalar_t coefficient = (fef.is_interior() ? 0.5 : 1.0);
    return coefficient * dv;
  }
};

// jump of trial
/**
  \brief This Expression represents the jump of a Trial function.

  It is used when in an expression the term jump(u) appears (where u is a Trial
  function). The specilization is needed to turn a particular Expression into a
  terminal leaf of the tree and perform the evaluation in a specific way
  breaking the recursion.
*/
template <>
class Expression<Jump, Trial> : public TerminalBase {
 public:
  using point_t = Traits::point_t;
  using scalar_t = Traits::scalar_t;
  /// Constructor for this specialized Expression.
  Expression(const Jump &, const Trial &) {}
  /**
    \brief Evaluate the average of the jump of a Trial function in a bilinear
    surface integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf.
    \param fef FeFace used in the expression context (for example to fetch the
    penalty or normal). \param u_side The side of the FeFace to be considered
    for the Trial function (either 1 or -1, it will multiply the normal to get
    the correct result) \param v_side The side of the FeFace to be considered
    for the Test function (either 1 or -1, it will multiply the normal to get
    the correct result) \param u Value of the basis function related to the
    Trial function. \param v Value of the basis function related to the Test
    function. \param du Value of the gradient of the basis related to the Trial
    function. \param dv Value of the gradient of the basis related to the Test
    function. \param x Physical point.
  */
  point_t eval(const FeFace &fef, const int &u_side, const int &v_side,
               const scalar_t &u, const scalar_t &v, const point_t &du,
               const point_t &dv, const point_t &x) const {
    return u * fef.normal() * u_side;
  }
};

// jump of test
/**
  \brief This Expression represents the jump of a Test function.

  It is used when in an expression the term jump(v) appears (where v is a Test
  function). The specilization is needed to turn a particular Expression into a
  terminal leaf of the tree and perform the evaluation in a specific way
  breaking the recursion.
*/
template <>
class Expression<Jump, Test> : public TerminalBase {
 public:
  using point_t = Traits::point_t;
  using scalar_t = Traits::scalar_t;
  /// Constructor for this specialized Expression.
  Expression(const Jump &, const Test &) {}
  /**
    \brief Evaluate the average of the jump of a Test function in a bilinear
    surface integration expression context.

    All the parameters must be pre-evaluated and passed to the expression
    evaluation context, that, according to the terminal leaf type, will extract
    the used information for that leaf.
    \param fef FeFace used in the expression context (for example to fetch the
    penalty or normal). \param u_side The side of the FeFace to be considered
    for the Trial function (either 1 or -1, it will multiply the normal to get
    the correct result) \param v_side The side of the FeFace to be considered
    for the Test function (either 1 or -1, it will multiply the normal to get
    the correct result) \param u Value of the basis function related to the
    Trial function. \param v Value of the basis function related to the Test
    function. \param du Value of the gradient of the basis related to the Trial
    function. \param dv Value of the gradient of the basis related to the Test
    function. \param x Physical point.
  */
  point_t eval(const FeFace &fef, const int &u_side, const int &v_side,
               const scalar_t &u, const scalar_t &v, const point_t &du,
               const point_t &dv, const point_t &x) const {
    return v * fef.normal() * v_side;
  }
};

}  // namespace wave