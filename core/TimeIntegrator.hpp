#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <functional>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>
#include "Traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/// \brief Helper class that contains the time integration options (it's an
/// aggregate).
struct TimeIntegratorOptions {
  using scalar_t = Traits::scalar_t;
  scalar_t time = 0.0;
  scalar_t timestep = 0.0;
  int ntimes = 0;
};

/**
  \brief A TimeIntegrator performs the second order explicit finite difference
  leap-frog scheme.

  This class is used by passing a solver (either an Eigen sparse solver or the
  BlockDiagonalSolver)as template parameter. With the compute method it will
  compute the factorization, while the solve method will solve the sequence of
  systems that give the solution time stepping, and it will return the solution
  at the final step. It takes a set of TimeIntegratorOptions to set the
  parameters. \tparam Solver The solver to be used to solve the sequence of
  linear systems.
*/
template <typename Solver>
class TimeIntegrator {
 public:
  using scalar_t = Traits::scalar_t;
  using dense_t = Traits::dense_t;
  using sparse_t = Traits::sparse_t;
  using solver_t = Solver;
  /// \brief Default constructor
  TimeIntegrator() = default;
  /// \brief Construct the time solver by passing the options
  TimeIntegrator(const TimeIntegratorOptions &options) : _options(options) {}
  /// \brief Get the options (and also enable to set them).
  TimeIntegratorOptions &options() { return _options; }
  /// \brief Get the options.
  const TimeIntegratorOptions &options() const { return _options; }
  /// \brief Get the the underlying solver (useful to initialize or configure
  /// the solver).
  solver_t &solver() { return _solver; }
  /// \brief Get the solver.
  const solver_t &solver() const { return _solver; }
  /// \brief Compute the factorization of the sparse matrix according to the
  /// solver.
  auto &compute(const sparse_t &matrix) {
#ifdef VERBOSE
    utl::Timer timer("factorization");
#endif

    _solver.compute(matrix);
    return *this;
  }

  /**
    \brief Apply the leap-frog scheme and get the solution at the final time

    The leap-frog scheme uses the mass matrix as problem matrix to be inverted,
    while the right hand side is built iteratively, starting from the initial
    condition vectors using the stiffness matrix and the set of boundary
    condition functions (here the space-time dependence is by product of
    functions). The number of right hand sides varies from at least a source
    term, and a variadic number representing other separate sources or Dirichlet
    and Neumann boundary terms. \tparam Function The source function type
    (normally a Field or lambda). \tparam Functions The boundary function type
    (normally Field's or lambdas). \param mass The mass matrix (sparse). \param
    stiffness The stiffness matrix (sparse). \param solution0 The L^2 projection
    of the initial condition on the function (this vector is obtained by solving
    M*x=U0 where U0 is the vector associated to the linear form u0 * v, u0 the
    initial condition function and v a test function. \param solution1 The L^2
    projection of the initial condition on the function time derivative (this
    vector is obtained by solving M*x=U1 where U1 is the vector associated to
    the linear form u1 * v, u1 the initial condition function derivative and v a
    test function. \param function The first right hand side vector (for example
    the source term, or the zero vector). Tie together in a tuple the vector of
    the space linear form and the lambda giving the time dependence. \param
    functions The other right hand side vectors (for example the boundary
    terms). Tie together in a tuple the vectors of the space linear forms and
    the lambdas giving the time dependence. \return Return the solution degrees
    of freedom in a vector (use nodal method of FeSpace to compute the solution
    at the mesh nodes, and use the vtk method of the Mesh to output it in VTK
    format and view with ParaView).
  */
  template <typename Function, typename... Functions>
  dense_t solve(const sparse_t &mass, const sparse_t &stiffness,
                const dense_t &solution0, const dense_t &solution1,
                const Function &function,
                const Functions &... functions) const {
#ifdef VERBOSE
    utl::Timer timer("solver");
#endif

    // init solution
    int size = mass.rows();
    dense_t approximation(size);
    approximation.setZero();

    // project initial condition vectors
    // this will also hold the most recent timestep
    dense_t approximation0 = _solver.solve(solution0);
    // this will also hold the oldest timestep
    dense_t approximation1 = _solver.solve(solution1);

    // init parameters
    scalar_t t = _options.time;
    scalar_t dt = _options.timestep;
    int niterations = _options.ntimes;

    // first step
    dense_t rhs =
        0.5 * dt * dt *
            (std::get<0>(function) * std::get<1>(function)(t) +
             ((std::get<0>(functions) * std::get<1>(functions)(t)) + ...) -
             stiffness * approximation0) +
        mass * (approximation0 + dt * approximation1);

    approximation = _solver.solve(rhs);
    int iteration = 1;

    // subsequent steps
    while (iteration < niterations) {
      // update
      approximation1 = approximation0;
      approximation0 = approximation;
      t += dt;

      rhs = dt * dt *
                (std::get<0>(function) * std::get<1>(function)(t) +
                 ((std::get<0>(functions) * std::get<1>(functions)(t)) + ...) -
                 stiffness * approximation0) +
            mass * (2 * approximation0 - approximation1);

      approximation = _solver.solve(rhs);
      ++iteration;
    }
    return approximation;
  }

 private:
  /// Time parameters options (initial time, timestep, number of steps)
  TimeIntegratorOptions _options{};
  /// The underlying solver
  solver_t _solver{};
};

}  // namespace wave
