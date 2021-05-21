#include <iostream>
#include "wave.hpp"

int main(int argc, const char *argv[]) {
  // main namespace
  using namespace wave;
  // use types without scope operator
  using namespace wave::Traits;
  // start the timer
  utl::Timer timer("total");

  // polynomial degree used
  constexpr int p = 2;
  // read the mesh from a file
  std::string meshfile = argv[1];
  Mesh Th(meshfile);
  // construct finite element space from mesh and degree
  FeSpace Vh(Th, p);

  // flag to set the symmetry of the bilinear for (true to compute only upper
  // diagonal part)
  bool sym = false;
  // expression terms
  Trial u;
  Test v;
  Normal n;
  Penalty sigma(10.0);

  // space dependent part of the source and dirichlet data
  Field f([](const point_t &x) {
    return -2 * (1 + x[0]) * x[1] * x[2] * std::exp(x[0]);
  });
  Field g([](const point_t &x) { return x[0] * x[1] * x[2] * std::exp(x[0]); });
  // time dependent part of the source and dirichlet data
  auto f_t = [](scalar_t t) { return std::sin(t); };
  auto g_t = [](scalar_t t) { return std::sin(t); };

  // create the forms
  BilinearForm M(Vh);
  BilinearForm K(Vh);
  LinearForm F(Vh);
  LinearForm G(Vh);

  // set the quadrature rules
  using Q3 = QuadratureKeastTetrahedron45;
  using Q2 = QuadratureDunavantTriangle25;

  // assign the mathematical expressions
  M.int3d<Q3>(u * v, sym);
  K.int3d<Q3>(dot(grad(u), grad(v)), sym);
  K.int2d<Q2>(-dot(avg(grad(u)), jump(v)) - dot(avg(grad(v)), jump(u)) +
                  sigma * dot(jump(u), jump(v)),
              {0, 1, 2, 3, 4, 5, 6}, sym);
  F.int3d<Q3>(f * v);
  G.int2d<Q2>(g * (-dot(grad(v), n) + sigma * v), {1, 2, 3, 4, 5, 6});

  // build matrices from triplets
  K.finalize();
  M.finalize();

  // initial condition on the function
  dense_t U0(Vh.ndof());
  U0.setZero();

  // initial condition on the function derivative (the L2 projection is obtained
  // solving M*x=U1.vector())
  Field uet_t0(
      [](const point_t &x) { return x[0] * x[1] * x[2] * std::exp(x[0]); });
  LinearForm U1(Vh);
  U1.int3d<Q3>(uet_t0 * v);

  // Time integration options (starting time, timestep, number of steps)
  TimeIntegratorOptions options{0.0, 0.0005, 2000};

  // set the solver to be used in the time stepping
  using solver_t =
      BlockDiagonalSolver<Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>>;

  // build time integrator
  TimeIntegrator<solver_t> integrator(options);

  // initialize the solver by passing the number of blocks and the dimension of
  // the square block
  int nlocal_dof = (p + 3) * (p + 2) * (p + 1) / 6;
  integrator.solver().init(Th.nelements(), nlocal_dof);

  // solver the problem with the leap-frog finite difference scheme
  integrator.compute(M.matrix().selfadjointView<Eigen::Upper>());
  dense_t uh = integrator.solve(M.matrix().selfadjointView<Eigen::Upper>(),
                                K.matrix().selfadjointView<Eigen::Upper>(), U0,
                                U1.vector(), std::tie(F.vector(), f_t),
                                std::tie(G.vector(), g_t));

  // compute the errors using the solutions at the final time
  scalar_t t1 = 1.0;

  Field ue_t1([&](const point_t &x) {
    return std::sin(t1) * x[0] * x[1] * x[2] * std::exp(x[0]);
  });
  Field due_t1(
      [&](const point_t &x) {
        return std::sin(t1) * (1 + x[0]) * x[1] * x[2] * std::exp(x[0]);
        ;
      },
      [&](const point_t &x) {
        return std::sin(t1) * x[0] * x[2] * std::exp(x[0]);
        ;
      },
      [&](const point_t &x) {
        return std::sin(t1) * x[0] * x[1] * std::exp(x[0]);
        ;
      });

  auto error = Vh.error<Q3>(uh, ue_t1, due_t1);

  // compute mesh size
  MeshFunction diameters(Th, [](Element e) { return e.diameter(); });
  scalar_t h = *std::max_element(diameters.begin(), diameters.end());

  // print info
  timer.stop();
  std::cout << "meshfile = " << meshfile << "\n";
  std::cout << "polynomial degree = " << p << "\n";
  std::cout << "number of dof = " << Vh.ndof() << "\n";
  std::cout << "mesh diameter = " << h << "\n";
  std::cout << "L2 error = " << error.first << "\n";
  std::cout << "H10 error = " << error.second << "\n";
  std::cout << std::endl;

  // output the solution in the VTK format after projecting it at the nodes
  Th.vtk("./tests/output.vtk", Vh.nodal(uh));

  return 0;
}
