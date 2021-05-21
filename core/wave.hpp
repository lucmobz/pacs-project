#pragma once

#include "Traits.hpp"
#include "core_include.hpp"
#include "functional_include.hpp"
#include "mesh_include.hpp"
#include "utilities_include.hpp"

/// \brief This is the main namespace used in the library (everything except for
/// utilities, which are in the utl namespace, is contained here)
namespace wave {}

/**
  \mainpage Wave Equation Library

  Here is described how the library functions (for general information or installation see [this link](README.md)). The tests are also commented.

  \section Workflow

  Here a brief overview on how to use the library is given.

  First of all a meshfile is needed, at the moment the library only supports meshes built with the VoroCrust software (not open-source, but a few meshes are included in \c ./tests/data, they are the unit cube with boundaries labeled starting from 1). 

  The main library header is \c wave.hpp. To use the \ref utl::Timer the utilities header is \c utilities_include.hpp. The timer will print the time elapsed when it exits a scope or by using the method \ref utl::Timer::stop.

  The first thing to do is to import the namespaces:
  \code
    using namespace wave // import the general namespace
    using namespace wave::Traits // import the types namespace
  \endcode

  A mesh is created by simply passing the path to the meshfile (a \ref wave::MeshLoader abstract class will extract the format and delegate the building of the mesh to a concrete loader, in the case of VoroCrust meshes it is \ref wave::MeshLoaderVoro).

  \code
    Mesh mesh("./data/voro1000.voro");
  \endcode

  A finite element space (\ref wave::FeSpace) related to the DG method is created by passing a mesh and a polynomial degree (or an \ref std::vector of polynomial degrees with the size of the elements and varying degrees):

  \code
    FeSpace fespace(mesh, 2);
  \endcode

  At this point bilinear and linear forms can be created on the finite element space:

  \code
    BilinearForm M(fespace); // Mass
    BilinearForm K(fespace); // Stiffness 
    LinearForm F(fespace); // Source term
    LinearForm G(fespace); // Dirichlet boundary datum or Neumann
  \endcode

  Now the expressions must be given. First of all create the terms in the expressions. The basic terms are:

  \code
    Trial u; // term representing a trial function
    Test v; // term representing a test function
    Normal n; // term representing the outward normal from the faces
    Penalty sigma(10.0); // term representing the penalty mesh function defined on faces (with underlying constant equal to 10)
  \endcode

  Constant scalars, vectors (or matrices in the diffusion part) are given directly as arithmetic types or **Eigen** types. Non constant coefficients, source or boundary functions are given as \ref wave::Field's. A \ref wave::Field is a term representing a function in an expression. To make such variable simply pass a few lambdas representing the component of the function.

  \code
    Field f1([](const point_t& x){ return x[0] + x[1] + x[2]; }); // scalar function
    Field f2([](const point_t& x){ return x[0]; },
             [](const point_t& x){ return x[1]; },
             [](const point_t& x){ return x[2]; }); // vector function (returns a point_t)
  \endcode

  An expression is given directly to the integrating methods of the forms, which also take as template parameter a quadrature rule (for example \ref wave::QuadratureGaussLegendre), here a few examples:

  \code
    using Q3 = QuadratureGaussLegendre<3, 5, 5, 5>; // gauss quadrature in 3D with 5 nodes per dimension
    using Q2 = QuadratureDunavantTriangle25; // some other quadrature rule in 2D
    
    Eigen::Vector3d beta{1, 1, 1};
    K.int3d<Q3>(dot(grad(u), grad(v)) + dot(beta, v) * u, true); // volume expression (true means the form is symmetric)
    K.int2d<Q2>(-dot(avg(grad(u)), jump(v))
                -dot(avg(grad(v)), jump(u)),
                + sigma * dot(jump(u), jump(v)), // surface expression
                {0, 1, 2, 3, 4, 5, 6}, // label of the faces boundaries where to integrate (0 is for interior faces)
                true); // symmetry flag

    F.int3d<Q3>(1 * v);
    F.int3d<Q3>(f * v);
    G.int2d<Q2>(g * (-dot(grad(v), n) + sigma * v), {1, 2, 3, 4, 5, 6});
  \endcode

  The boundary labels are given in a \ref std::set, 0 must be included for \ref BilinearForm::int2d for interior faces, the other numbers are the labels of the boundaries as given in the meshfile.

  To get the the **Eigen** types corresponding to the discretized linear system there are the \ref wave::BilinearForm::matrix and \ref wave::LinearForm::vector methods:

  \code
    K.finalize(); // transform the triplets into a sparse matrix
    K.matrix(); // Eigen sparse matrix
    F.vector(); // Eigen dense vector
  \endcode

  At this point a manual loop can be created with the **Eigen** types to implement a solution of the wave equation (if the problem is time dependent, if not one can solve the problem with an **Eigen** solver). A built-in wave equation solver that uses the second order finite difference leap-frog scheme is included, the syntax is the following as the linear system sequence is block diagonal:

  \code
    TimeIntegratorOptions options{0.0, 0.0005, 2000}; // options aggregate class, t0, dt, num of timesteps

    using solver_t = BlockDiagonalSolver<Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>>; // set the solver to be used in the time stepping
      
    TimeIntegrator<solver_t> integrator(options); // create time solver
  \endcode

  Here the custom \ref wave::BlockDiagonalSolver is used. It is templated on the **Eigen** dense solver that is used to solve the block systems. The integrator can take instead also a generic **Eigen** sparse solver as template parameter.

  To solve, the initial conditions are needed. They are given by the L^2 projection of the initial condition functions onto the discrete space, which can be obtained by solving the system M*x=U, where U is the linear form associated to the function and M is the mass matrix, here is an example:
  
  \code
    Field u_t0([](const point_t &x) { return 1 * x[0]; }); // some random initial condition
    LinearForm U(fespace); // get the associated linear form
    U.int3d<Q3>(u_t0 * v); 
    U.vector() // this is the initial condition vector before projection

    BilinearForm M(fespace);
    M.int3d<Q3>(u * v);

    // solve the system ...
  \endcode
  
  The \ref wave::TimeIntegrator will need the initial condition vector to solve (it will compute the projection itself). To solve, use the appropriate method after having initialized or configured the underlying solver.

  \code
    integrator.solver().init(mesh.nelements(), block_size); // init or config the solver using the solver() reference method

    // solve the problem with the leap-frog finite difference scheme
    integrator.compute(M.matrix().selfadjointView<Eigen::Upper>());

    Eigen::VectorXd uh = integrator.solve(M.matrix().selfadjointView<Eigen::Upper>(),
                                          K.matrix().selfadjointView<Eigen::Upper>(), U0.vector(),
                                          U1.vector(), std::tie(F.vector(), f_t),
                                          std::tie(G.vector(), g_t));
  \endcode
  
  Some comments. If the symmetry flag is used in all the integrals the matrices will be upper triangular (only half of them will be computed) and will need to be passed as self-adjoint views, if an integral is not symmetric for a given form, the matrix will be double sided. U0.vector() and U1.vector() are the **Eigen** types associated with the intial conditions (if one has the zero initial condition, the zero vector can be passed instead). Source and boundary functions are supported if their space-time dependence is by product (like x*y*z*sin(t) is supported). In the linear forms only the space part must be given. In the solving, the space vector and the lambda consisting in the time term must be tied together. This is done for as many source or boundary terms are present (the function is variadic).

  To finish the solution can be exported in the VTK format (viewable with ParaView) and the error computed if the exact solution is known:

  \code
    Field u_exact([](const point_t& x){ return std::exp(x[0]); }); // some random exact solution

    Field gradu_exact([](const point_t& x){ return std::exp(x[0]); },
                      [](const point_t& x){ return 0.0; }, 
                      [](const point_t& x){ return 0.0; }); // its gradient

    mesh.vtk("output.vtk", fespace.nodal(uh)); // compute solution at mesh nodes and export it
    auto error = fespace.error<Q3>(uh, u_exact, gradu_exact); // return the error as a std::pair of scalars (L^2, H^1_0)

    std::cout << error.first << "\n" << error.second << "\n";
  \endcode 


*/