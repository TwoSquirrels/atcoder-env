#pragma once

#include <risu/prelude.hpp>

#include <algorithm>
#include <ranges>
#include <type_traits>

// [l, r) as a range, empty when r < l
template <typename T, typename U> inline auto iota_range(T l, U r) {
  using V = std::common_type_t<T, U>;
  return std::views::iota(V(l), std::max(V(l), V(r)));
}
