#pragma once

#include <muParser.h>
#include <array>
#include <string>

/// muParser function (not used)
template <typename Scalar, char... VARIABLES>
class muFunction {
 public:
  using scalar_t = Scalar;
  static constexpr int NVARIABLES = sizeof...(VARIABLES);

  muFunction() { _init(); }
  muFunction(const muFunction<scalar_t, VARIABLES...>& other)
      : _parser{other._parser} {
    _init();
  }
  muFunction(const std::string& expr) {
    _init();
    set(expr);
  }

  muFunction& operator=(const muFunction& other) {
    if (this == &other) return *this;
    _parser = other._parser;
    _init();
    return *this;
  }

  void set(const std::string& expr) { _parser.SetExpr(expr); }

  template <typename Point>
  scalar_t operator()(const Point& point) const {
    for (int i = 0; i < NVARIABLES; ++i) _values[i] = point[i];
    return _parser.Eval();
  };

 protected:
  mu::Parser _parser{};
  mutable std::array<scalar_t, NVARIABLES> _values{};

  void _init() {
    int i = 0;
    ((_parser.DefineVar(std::string(1, VARIABLES), &_values[i]), ++i), ...);
  }
};