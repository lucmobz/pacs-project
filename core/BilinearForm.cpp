#include "BilinearForm.hpp"
#include <iterator>

namespace wave {

void BilinearForm::finalize() {
  std::vector<triplet_t> triplets;

  int ntriplets = 0;
  for (const auto &ts : _triplets) {
    ntriplets += ts.second.size() * (1 + ts.first);
  }
  triplets.reserve(ntriplets);

  for (auto &ts : _triplets) {
    if (!is_symmetric() && ts.first) {
      for (const auto &t : ts.second) {
        if (t.row() != t.col())
          triplets.emplace_back(t.col(), t.row(), t.value());
      }
    }
    triplets.insert(triplets.end(), std::make_move_iterator(ts.second.begin()),
                    std::make_move_iterator(ts.second.end()));
  }
  _triplets.clear();

  _matrix.setFromTriplets(triplets.begin(), triplets.end());
  _matrix.prune(_matrix.coeff(0, 0));
}

}  // namespace wave