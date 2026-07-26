#pragma once

#include <risu/prelude.hpp>

#include <array>

#if __has_include(<atcoder/modint>)
#  include <atcoder/modint>

template <typename T = atcoder::modint>
auto mint_inv(T x) -> T {
  constexpr int memo_limit = 1 << 24;
  static std::array<T, memo_limit> memo;
  if (x.val() >= memo_limit) return x.inv();
  if (memo[x.val()] == 0) memo[x.val()] = x.inv();
  return memo[x.val()];
}
#endif // <atcoder/modint>
