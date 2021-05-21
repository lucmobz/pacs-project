#pragma once

namespace utl {
constexpr long unsigned factorial(long unsigned number) {
  if (number == 0lu)
    return 1lu;
  else if (number == 1lu)
    return 1lu;
  else
    return number * factorial(number - 1lu);
}
}  // namespace utl