#pragma once

#include <risu/prelude.hpp>

template <typename T, typename U> inline auto chmin(T &&a, const U b) -> bool {
  const auto compare = a > b;
  if (compare) a = b;
  return compare;
}
template <typename T, typename U> inline auto chmax(T &&a, const U b) -> bool {
  const auto compare = a < b;
  if (compare) a = b;
  return compare;
}

