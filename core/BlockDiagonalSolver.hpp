#pragma once

// #include <omp.h> (no omp functions are used)
#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>
#include "Traits.hpp"

// disable eigen automatic parallalezation to use openmp on groups of blocks
// instead than in each single block
#define EIGEN_DONT_PARALLELIZE

namespace wave {
/**
  \brief A BlockDiagonalSolver is a solver that exploit the structure of a
  matrix

  This solver is meant to be applied to a sparse matrix that has a block
  diagonal structure. \tparam Solver The Eigen solver to be used to solve the
  block systems (for example Eigen::LLT<MatrixXd, Eigen::Upper> for symmetric
  positive definite block systems stored in the upper triangular part of a
  matrix)
*/
template <typename Solver>
class BlockDiagonalSolver;

/**
  \brief Construct the solver by taking a concrete Eigen solver to solve the
  block systems \tparam Solver The concrete Eigen solver to be used on blocks
  \tparam Matrix This is deduced from the Solver.
  \tparam OPTIONS These are deduces from the Solver.
*/
template <template <typename, int...> typename Solver, typename Matrix,
          int... OPTIONS>
class BlockDiagonalSolver<Solver<Matrix, OPTIONS...> > {
 public:
  using scalar_t = Traits::scalar_t;
  using dense_t = Traits::dense_t;
  using sparse_t = Traits::sparse_t;
  using matrix_t = Matrix;
  using block_solver_t = Solver<Matrix, OPTIONS...>;

  /// \brief Default constructor
  BlockDiagonalSolver() = default;
  /// \brief Construct solver from the sizes of the blocks.
  BlockDiagonalSolver(const std::vector<int> &block_sizes) {
    init(block_sizes);
  }
  /// \brief Construct solver from the sizes of the blocks (by move)
  BlockDiagonalSolver(std::vector<int> &&block_sizes) {
    init(std::move(block_sizes));
  }
  /// \brief Construct solver giving the number of blocks and the number of rows
  /// of a block (assuming all blocks are homogeneous)
  BlockDiagonalSolver(int nblocks, int size)
      : BlockDiagonalSolver(std::vector<int>(nblocks, size)) {}

  /// \brief Initialize an uninitialized solver form the block sizes
  void init(const std::vector<int> &block_sizes) {
    int size = _block_sizes.size();

    _nblocks = size;
    _block_sizes = block_sizes;
    _offsets.resize(size);
    _offsets[0] = 0;
    std::partial_sum(_block_sizes.begin(), _block_sizes.end() - 1,
                     _offsets.begin() + 1);
    _solvers.resize(size);
  }
  /// \brief Initialize an uninitialized solver form the block sizes (by move)
  void init(std::vector<int> &&block_sizes) {
    int size = block_sizes.size();

    _nblocks = size;
    _block_sizes = std::move(block_sizes);
    _offsets.resize(size);
    _offsets[0] = 0;
    std::partial_sum(_block_sizes.begin(), _block_sizes.end() - 1,
                     _offsets.begin() + 1);
    _solvers.resize(size);
  }
  /// \brief Initialize an uninitialized solver form the number of blocks and an
  /// homogeneous size.
  void init(int nblocks, int size) {
    this->init(std::vector<int>(nblocks, size));
  }

  /// \brief Compute the factorizations of the blocks of a block diagonal sparse
  /// matrix.
  auto &compute(const sparse_t &matrix) {
#pragma omp parallel for default(none) \
    shared(_nblocks, _block_sizes, _offsets, _solvers, matrix)
    for (int i = 0; i < _nblocks; ++i) {
      int size = _block_sizes[i];
      int offset = _offsets[i];
      _solvers[i].compute(matrix_t(matrix.block(offset, offset, size, size)));
    }

    return *this;
  }
  /// \brief Solve the system using the precomputed factorizations
  /// \param vector The vector that represents the right hand side in the linear
  /// system.
  dense_t solve(const dense_t &vector) const {
    dense_t solution(vector.size());

#pragma omp parallel for default(none) \
    shared(_nblocks, _block_sizes, _offsets, _solvers, solution, vector)
    for (int i = 0; i < _nblocks; ++i) {
      int size = _block_sizes[i];
      int offset = _offsets[i];
      solution.segment(offset, size) =
          _solvers[i].solve(vector.segment(offset, size));
    }

    return solution;
  }

 private:
  /// The total number of blocks
  int _nblocks{0};
  /// The vector that contains the size of each block
  std::vector<int> _block_sizes;
  /// The vector that contains the location of the left upper corner of each
  /// block
  std::vector<int> _offsets;
  /// The vector that contains the factorizations of each block and block
  /// solvers.
  std::vector<block_solver_t> _solvers;
};

}  // namespace wave

#undef EIGEN_DONT_PARALLELIZE