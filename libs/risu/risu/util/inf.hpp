#pragma once

#include <risu/prelude.hpp>

#include <limits>

template <typename T> constexpr auto inf() -> T {
  if constexpr (std::numeric_limits<T>::has_infinity) return std::numeric_limits<T>::infinity();
  return std::numeric_limits<T>::max() / 2.125L;
}

