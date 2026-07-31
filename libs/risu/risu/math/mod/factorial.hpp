#pragma once

#include <risu/prelude.hpp>
#include <risu/math/mod/mint_inv.hpp>

#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>
#if __has_include(<atcoder/modint>)
#  include <atcoder/modint>
#endif

template <typename T> auto fact(int n, bool inv = false) -> T {
  assert(n >= 0);
  static auto factorials = std::vector<std::pair<T, T>>{ { 1, 1 } };
  for (auto i = std::ssize(factorials); i <= n; ++i) factorials.emplace_back(i * factorials[i - 1].first, 0);
  if (inv && factorials[n].second == 0) {
    if constexpr (std::is_integral_v<T>) factorials[n].second = factorials[n].first <= 1 ? 1 : 0;
#if __has_include(<atcoder/modint>)
    else if constexpr (atcoder::internal::is_modint<T>::value) factorials[n].second = mint_inv(factorials[n].first);
#endif // <atcoder/modint>
    else factorials[n].second = 1 / factorials[n].first;
  }
  return inv ? factorials[n].second : factorials[n].first;
}

