#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "GetPot"
#include "wave.hpp"

int main(int argc, char **argv) {
  using namespace wave;
  // use library types without scope operator
  using namespace wave::Traits;
  // set quadrature rules to have exact integrals
  using Q3 = QuadratureGaussLegendre<3, 10, 10, 10>;
  using Q2 = QuadratureGaussLegendre<2, 10, 10>;
  // set the solver for the sequence of linear systems
  using solver_t =
      BlockDiagonalSolver<Eigen::LLT<Eigen::MatrixXd, Eigen::Upper>>;

  // read input data
  GetPot cl(argc, argv);
  GetPot ifl(cl.follow("./pconv.pot", 2, "-f", "--file"));

  std::string meshfile = ifl("data", "./data/voro27.voro");

  std::vector<int> degrees;
  for (int i = 0; i < ifl.vector_variable_size("p"); ++i)
    degrees.emplace_back(ifl("p", 0, i));

  TimeIntegratorOptions options(ifl("t0", 0.0), ifl("dt", 0.0005),
                                ifl("nt", 2000));

  // loop over polynomial degrees
  for (const auto &degree : degrees) {
    utl::Timer timer("total");
    // create problem
    int p = degree;
    Mesh Th(meshfile);
    FeSpace Vh(Th, p);

    // compute mesh size
    MeshFunction diameters(Th, [](Element e) { return e.diameter(); });
    scalar_t h = *std::max_element(diameters.begin(), diameters.end());

    // create forms
    BilinearForm M(Vh), K(Vh);
    LinearForm F(Vh), G(Vh);

    // expression terms
    bool sym = true;
    Trial u;
    Test v;
    Penalty sigma(10.0);
    Normal n;

    Field f([](const point_t &x) {
      return -2 * (1 + x[0]) * x[1] * x[2] * std::exp(x[0]);
    });
    Field g(
        [](const point_t &x) { return x[0] * x[1] * x[2] * std::exp(x[0]); });

    M.int3d<Q3>(u * v, sym);
    K.int3d<Q3>(dot(grad(u), grad(v)), sym);
    K.int2d<Q2>(-dot(avg(grad(u)), jump(v)) - dot(avg(grad(v)), jump(u)) +
                    sigma * dot(jump(u), jump(v)),
                {0, 1, 2, 3, 4, 5, 6}, sym);

    F.int3d<Q3>(f * v);
    G.int2d<Q2>(g * (-dot(grad(v), n) + sigma * v), {1, 2, 3, 4, 5, 6});

    K.finalize();
    M.finalize();

    // functions time dependent part and initial conditions
    auto f_t = [](scalar_t t) { return std::sin(t); };
    auto g_t = [](scalar_t t) { return std::sin(t); };

    dense_t U0(Vh.ndof());
    U0.setZero();

    Field uet_t0(
        [](const point_t &x) { return x[0] * x[1] * x[2] * std::exp(x[0]); });
    LinearForm U1(Vh);
    U1.int3d<Q3>(uet_t0 * v);

    // time integration
    TimeIntegrator<solver_t> integrator(options);

    int block_size = (p + 3) * (p + 2) * (p + 1) / 6;
    integrator.solver().init(Th.nelements(), block_size);

    integrator.compute(M.matrix().selfadjointView<Eigen::Upper>());

    dense_t uh = integrator.solve(M.matrix().selfadjointView<Eigen::Upper>(),
                                  K.matrix().selfadjointView<Eigen::Upper>(),
                                  U0, U1.vector(), std::tie(F.vector(), f_t),
                                  std::tie(G.vector(), g_t));

    // exact solution at final time
    scalar_t t1 = 1.0;

    Field ue_t1([&](const point_t &x) {
      return std::sin(t1) * x[0] * x[1] * x[2] * std::exp(x[0]);
    });

    Field due_t1(
        [&](const point_t &x) {
          return std::sin(t1) * (1 + x[0]) * x[1] * x[2] * std::exp(x[0]);
        },
        [&](const point_t &x) {
          return std::sin(t1) * x[0] * x[2] * std::exp(x[0]);
        },
        [&](const point_t &x) {
          return std::sin(t1) * x[0] * x[1] * std::exp(x[0]);
        });

    // error at final time
    auto error = Vh.error<Q3>(uh, ue_t1, due_t1);

    // print info
    timer.stop();
    std::cout << "meshfile = " << meshfile << "\n";
    std::cout << "polynomial degree = " << p << "\n";
    std::cout << "number of dof = " << Vh.ndof() << "\n";
    std::cout << "mesh diameter = " << h << "\n";
    std::cout << "L2 error = " << error.first << "\n";
    std::cout << "H10 error = " << error.second << "\n";
    std::cout << std::endl;
  }

  return 0;
}