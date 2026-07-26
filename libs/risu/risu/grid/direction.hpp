#pragma once

#include <risu/prelude.hpp>

template <typename T> constexpr auto sin90(T theta90) -> int {
  if (theta90 % 2 == 0) return 0;
  return (theta90 % 4 + 4) % 4 < 2 ? 1 : -1;
}
template <typename T> constexpr auto cos90(T theta90) -> int {
  return sin90(theta90 + 1);
}
template <typename T> constexpr auto sin45(T theta45) -> int {
  if (theta45 % 4 == 0) return 0;
  return (theta45 % 8 + 8) % 8 < 4 ? 1 : -1;
}
template <typename T> constexpr auto cos45(T theta45) -> int {
  return sin45(theta45 + 2);
}

