#pragma once

#include <risu/prelude.hpp>
#include <risu/math/mod/factorial.hpp>

template <typename T> auto perm(long long n, int k) -> T {
  if (n <= 5000000) return fact<T>(n) * fact<T>(n - k, true);
  auto res = T(1);
  for (auto i = n; i > n - k; --i) res *= i;
  return res;
}
template <typename T> auto comb(long long n, int k) -> T {
  if (n - k < k) k = n - k;
  return perm<T>(n, k) * fact<T>(k, true);
}

