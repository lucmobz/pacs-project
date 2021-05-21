#pragma once

#include <array>
#include "Traits.hpp"
#include "utilities_include.hpp"

namespace wave {
/// \brief Sample 2D quadrature class to be used for testing. Any quadrature
/// with the 4 static methods explained in QuadratureGaussLegendre is
/// admissible.
// see John Burkardt home page
class QuadratureDunavantTriangle25 {
 private:
  // this class is a stateless singleton
  QuadratureDunavantTriangle25() = delete;
  static constexpr int DIM = 2;
  static constexpr int SIZE = 25;
  static constexpr int EXACTNESS = 10;

 public:
  using scalar_t = Traits::scalar_t;
  using qpoint_t = Traits::qpoint_t<DIM>;
  using weights_t = std::array<scalar_t, SIZE>;
  using nodes_t = std::array<qpoint_t, SIZE>;

  /// \brief Get the number of quadrature nodes and weights.
  static constexpr auto size() -> int { return SIZE; }
  /// \brief Get the order of exactness of the quadrature rule.
  static constexpr auto exactness() -> int { return EXACTNESS; }
  /// \brief Get the measure of the reference simplex.
  static constexpr auto measure() -> scalar_t {
    return 1 / static_cast<scalar_t>(utl::factorial(DIM));
  }
  /// \brief Get the quadrature weights.
  static constexpr auto weights() -> weights_t {
    return weights_t{0.045408995191377,  0.0183629788782335, 0.0183629788782335,
                     0.0183629788782335, 0.022660529717764,  0.022660529717764,
                     0.022660529717764,  0.036378958422710,  0.036378958422710,
                     0.036378958422710,  0.036378958422710,  0.036378958422710,
                     0.036378958422710,  0.0141636212655285, 0.0141636212655285,
                     0.0141636212655285, 0.0141636212655285, 0.0141636212655285,
                     0.0141636212655285, 0.0047108334818665, 0.0047108334818665,
                     0.0047108334818665, 0.0047108334818665, 0.0047108334818665,
                     0.0047108334818665};
  }
  /// \brief Get the quadrature nodes.
  static auto nodes() -> nodes_t {
    return nodes_t{qpoint_t{1. / 3., 1. / 3.},
                   qpoint_t{0.028844733232685, 0.485577633383657},
                   qpoint_t{0.485577633383657, 0.028844733232685},
                   qpoint_t{0.485577633383657, 0.485577633383657},
                   qpoint_t{0.781036849029926, 0.109481575485037},
                   qpoint_t{0.109481575485037, 0.781036849029926},
                   qpoint_t{0.109481575485037, 0.109481575485037},
                   qpoint_t{0.141707219414880, 0.307939838764121},
                   qpoint_t{0.141707219414880, 0.550352941820999},
                   qpoint_t{0.307939838764121, 0.550352941820999},
                   qpoint_t{0.307939838764121, 0.141707219414880},
                   qpoint_t{0.550352941820999, 0.307939838764121},
                   qpoint_t{0.550352941820999, 0.141707219414880},
                   qpoint_t{0.025003534762686, 0.246672560639903},
                   qpoint_t{0.025003534762686, 0.728323904597411},
                   qpoint_t{0.246672560639903, 0.728323904597411},
                   qpoint_t{0.246672560639903, 0.025003534762686},
                   qpoint_t{0.728323904597411, 0.246672560639903},
                   qpoint_t{0.728323904597411, 0.025003534762686},
                   qpoint_t{0.009540815400299, 0.066803251012200},
                   qpoint_t{0.009540815400299, 0.923655933587500},
                   qpoint_t{0.066803251012200, 0.009540815400299},
                   qpoint_t{0.066803251012200, 0.923655933587500},
                   qpoint_t{0.923655933587500, 0.066803251012200},
                   qpoint_t{0.923655933587500, 0.009540815400299}};
  }
};
}  // namespace wave